#! /usr/bin/perl -w

use strict;
sub  error_state_check (@ ); # defined at the bottom
sub make_fake_top (@);

defined $ARGV[0] ||
    die "Usage: gromacs.pl <pdb name root>/no_protein [<list of ligand names> -remd|-min].\n";

$| = 1; # turn on autoflush;

# remd flag must be the last, for now
my $remd = 0;
($ARGV[$#ARGV] eq "-remd" ) && ($remd = 1);
my $minimization = 0;
($ARGV[$#ARGV] eq "-min" ) && ($minimization = 1);

##############################################################################
##############################################################################
my $np = 0;

my $cmd_log = 1;
my $pdb2gro;
my $gromacs_path = "/usr/local/bin";
my $mpirun = "/usr/local/bin/mpirun";
my ($mdrun, $mdrun_mpi, $grompp);
my $gromacs_run;
my $groc;
my $grompp_root;
my $protein = 1;

my $hostname = `hostname`;
chomp $hostname;

   
$gromacs_run = "$gromacs_path/mdrun -nt 1";
$grompp      = "$gromacs_path/grompp";
$groc        = "/home/ivanam/perlscr/gromacs/gro_concat.pl";
$pdb2gro     = "/home/ivanam/perlscr/translation/pdb2gro.pl";


if ($hostname eq "donkey" || $hostname eq "bii-sharonc") {

    $gromacs_path = "/usr/bin";
    $grompp       = "/usr/bin/grompp";
    $gromacs_run  = "/usr/bin/mdrun -nt 1";
    $mdrun_mpi    = "/usr/bin/mdrun_mpi";
    $mpirun       = "/usr/bin/mpirun";
}

##############################################################################
##############################################################################

my $name        = $ARGV[0]; # expect ok pdb file called $name.pdb present
my $forcefield  = "amber99sb";
#my $forcefield = "oplsaa";
#my $forcefield = "gmx";
#my $box_type   = "cubic";
my $box_type    = "triclinic";
my $box_edge    =  0.7; # distance of the box edge from the molecule in nm 
my $neg_ion     = "Cl"; # theses names depend on the choice of the forcefield (check "ions.itp")
my $pos_ion     = "Na";
my $water      = "tip3p";
my $genion_solvent_code = 12;

if (  $forcefield  eq  "oplsaa" || $forcefield  eq "amber99sb" ) {
    $neg_ion  = "CL";
    $pos_ion  = "NA";
}

##############################################################################
##############################################################################
## in  ~/moldyn/gromacs/tutorials/generic_mdps
foreach ("em.mdp", "pr.mdp", "md.mdp" ) {
    ( -e $_ ) || die "\n$_ not found.\n\n";
}

if ( $name eq "no_protein") {
    $protein = 0;
    (  defined $ARGV[1] ) ||
	die "No protein, no small molecule ... what are we doing here?\n";
    $box_edge    =  1.5;
} else {
    ( -e "$name.pdb") || die "$name.pdb not found\n";
}

##############################################################################
##############################################################################
$cmd_log &&  ( open (CMD_LOG, ">cmd.log") || die "Cno cmd.log: $!.\n");


##############################################################################
##############################################################################
## check if itp and gro given for all ligands
my @ligand_names = ();
my %multiplicity = (); # how many time each ligand appears
my $mult;
my $ligand_name;
my $file;
# this nomenclature corresponds to Gromacs forcefield
my %ion = ( "K", 1, "Na",  1, "Ca", 2,  "Mg", 2, "Cl", -1, "Zn", 2,
	    "NA",  1, "CA", 2,  "MG", 2, "CL", -1, "ZN", 2, 
	    "na",  1, "ca", 2,  "mg", 2, "cl", -1, "zn", 2 );
# is this correct?
my %opls_ion_names = ( "mg", "MG2+", "ca", "CA2+",   "li", "LI+",   
		       "na", "NA+",   "k", "K+",   "rb", "Rb+",   
		       "cs", "Cs+",   "f", "F-",   "cl", "CL-",   
		       "br", "BR-",   "i", "I-");

if ( defined $ARGV[1] ) {
    my @aux;
    open ( LF, "<$ARGV[1]") ||
	die "Cno $ARGV[1]: $!.\n";
    while ( <LF> ) {
	chomp;
	@aux = split;
	push @ligand_names, $aux[0];
	if ( defined  $aux[1] ) {
	    $multiplicity{$aux[0]} = $aux[1];
	} else {
	    die "Number of instances of molecule $aux[0] not defined in $ARGV[1].\n";
	}
       
    }
    $protein || ($name = join "_", @ligand_names);

    close LF;
    foreach $ligand_name ( @ligand_names ) {
	
	if ( $ligand_name  eq "water" ) {
            #all waters are understood to be in the same file; 
            #itp is standard for water
	    $file = "water.gro"; 
	    ( -e $file) || die "$file not found in ".`pwd`;

	} else {
	    for ($mult=1; $mult <= $multiplicity{$ligand_name}; $mult ++ ) {
		my $ln = (lc $ligand_name);
		my $name_root = $ln;
		my $cmd;
		

		if ( $multiplicity{$ligand_name} > 1 ) {
		    $name_root = $ln.$mult;
		}

		$file = $name_root.".gro";
		next if ( -e $file);
		( -e "$name_root.pdb" ) || die "Cno $name_root.pdb: $! .\n";

		# handle ions
		if ( defined $ion{$ligand_name} ) {
		    if (  $forcefield  eq  "oplsaa" ) {
			my $tmp = "tmp.ion";
			my ($ion_type, $new_ion_type, $new_line);
			my ($serial, $res_seq, $x, $y, $z);
			open (ION_IF, "<$name_root.pdb");
			open (ION_OF, ">$name_root.gro");
			while ( <ION_IF> ) {

			    $ion_type = substr $_, 12, 4;  $ion_type =~ s/\s//g;
			    $serial   = substr $_,  6, 5;  $serial =~ s/\s//g;
			    $res_seq  = substr $_, 22, 4;  $res_seq=~ s/\s//g;
			    $x = substr $_, 30, 8; $x =~ s/\s//g;
			    $y = substr $_, 38, 8; $y =~ s/\s//g;
			    $z = substr $_, 46, 8; $z =~ s/\s//g;
 
			    (defined $opls_ion_names{lc $ion_type}) || 
				die "opls name for  $ion_type not defined.".
				"\n(at least not in $0).\n";
			   
			    $new_ion_type = sprintf "%-5s",  $opls_ion_names{lc $ion_type};
			   
			    printf ION_OF  "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
			    $res_seq, $new_ion_type,  $ion_type, $serial, $x/10, $y/10, $z/10;
			}
			close ION_IF;
			close ION_OF;

		    } else {
			$cmd = "$pdb2gro <  $name_root.pdb >  $name_root.gro ";
			(system $cmd) && die "Error running $cmd\n";
			( $cmd_log ) &&  print CMD_LOG "$cmd\n";
		    }

		# handle other types of ligands
		} else {

		    ( -e "$name_root.gro") || die "$name_root.gro must be prepared in advance\n";
		    ( -e "$ln.itp") || die "$ln.itp must be prepared in advance\n";
		}
	    }
	    if ( ! defined $ion{$ligand_name} ) { # ion "topology" is included by default
		# all other molecules have their own gro and itp files
		$file = lc $ligand_name.".itp"; 
		( -e $file) || die "$file not found in ".`pwd`;
	    }
	}
    } 
} 

############################################################################## 
##############################################################################
my $command;
my $ret;
my $program;
my @aux;
my $charge;
my $input_system_for_md;
my $log;

###############
# process pdb into gro and topology files
###############


if ( ! -e "$name.top" || ! -e "$name.gro" )  {

    if ( $protein ) {
	$program = "$gromacs_path/pdb2gmx";
	$log = "pdb2gmx.log";
	print "\t runnning $program \n";
	# -ignh instructs pdb2gmx to ingore H and place its own
	$command  = "$program  -ignh -ff $forcefield -water $water -f $name.pdb -o $name.gro -p $name.top ";
	( -e "pdb2gmx_in") && ( `rm pdb2gmx_in`); 



	if ( -e "ssbridges" ) {
	    # ssbridges has one "y\n" for each ssbridge (sequentailly
	    # by the first cysteine) that we want to maintain
	    # I guess "n\n" for the ones we do not want too
	    `cat ssbridges >> pdb2gmx_in`;
	    $command  .= " -ss";
	} 

	if ( -e "termini" ) {
	    # termini has the following format: "2\n2\n" for each chain 
	    # (4 x 2 for two chains and so on) 
	    `cat termini >> pdb2gmx_in`;
	    $command  .= " -ter";
	}
	
	(-e "pdb2gmx_in") && ($command  .= " < pdb2gmx_in");
	$command  .= " > $log 2>&1";
	($cmd_log)  &&  print CMD_LOG "$command\n";

	(system $command)
	    && die "Error:\n$command\nerror"; # looks like it exits with something on success
	error_state_check ( $log, ("masses will be determined based on residue and atom names"));
    }

    if ( @ligand_names ) {
	my $written;
	####################################
	# concatenate gro files
	if ( ! $protein && @ligand_names == 1 ) {
	    ( -e "$name.gro" ) || die "$name.gro not found\n";
	} else {
	    if ( $protein ) {
		`cp $name.gro $name.orig.gro`;
		$command = "$groc  $name.gro ";
	    } else  {
		$command = "$groc ";
	    }
	    foreach $ligand_name ( @ligand_names ) {
		if ($ligand_name eq "water") {
		    $command .=  $ligand_name.".gro ";
		    next;
		}
		if ( $multiplicity{$ligand_name} == 1 ) {
		    $command .=  (lc $ligand_name).".gro ";
		} else {
		    for ($mult=1; $mult <= $multiplicity{$ligand_name}; $mult ++ ) {
			$command .=  (lc $ligand_name).$mult.".gro ";
		    }
		}
	    }
	    $command .= " > tmp";
	    my $ret;
	    $ret = system $command ;
	    $ret	&& die "Error:\n$command\n$ret"; 
	    ( $cmd_log ) &&  print CMD_LOG "$command\n";
	    `mv tmp $name.gro`;
	}

	####################################
	# change the topology (*.top) file
	if ( ! $protein ) {
	    make_fake_top ($name);
	}
	open ( TOP, "<$name.top") ||
	    die "Error opening $name.top: $!.\n";

	open ( OF, ">tmp") ||
	    die "Error opening tmp output file: $!.\n";
	$written = 0;
	while ( <TOP> ) {
	    print OF;
	    if ( !$written && /^\#include/ ) {
		foreach $ligand_name ( @ligand_names ) {
		    next if  ($ligand_name eq "water");
		    next if ( defined $ion{$ligand_name} );
		    print OF "#include \"". lc $ligand_name.".itp\"\n";
		}
		$written = 1;
	    } elsif ( /^\[ molecules \]/ ) {
		while ( <TOP> ) {
		    if ( /\S/ ) {
			print OF;
		    }
		}
		foreach $ligand_name ( @ligand_names ) { 
		    if  ($ligand_name eq "water"){
			 print OF "SOL     $multiplicity{$ligand_name}\n";
		    } else { 
			if ( $forcefield eq "oplsaa" && defined $opls_ion_names{lc $ligand_name} ) {
			    print OF  $opls_ion_names{lc $ligand_name}. "  $multiplicity{$ligand_name}\n";
			} else {
			    print OF  $ligand_name."     $multiplicity{$ligand_name}\n";
			}
		    }
		}
		print OF "\n";
	    }
	}
	close OF;
	close TOP;
	`mv tmp $name.top`;

	
    }

} else {
    print "\t gro and top files found\n";
}


###############
# place the sytem in a box
###############
$file = "$name.box.gro";
if (!  -e  $file) {
    $program = "$gromacs_path/editconf";
    $log = "editconf.log";
    print "\t runnning $program \n";
    $command = "$program -f $name.gro -o $file  -bt $box_type -d $box_edge -c > $log  2>&1";
    # -c is the centering command
    (system $command)
	&& die "Error running\n$command\n"; 
    ( $cmd_log ) &&  print CMD_LOG "$command\n";
    # some bug in editconf: if the box triclinic it reports "no boxtype specified,
    # but constructs something which has all three unit cell vectors perpendicular 
    # and of different lengths
    error_state_check ( $log, ("masses will be determined based on residue and atom names",
				    "No boxtype specified"));
} else {
    print "\t $file found\n"; 
}


###############
# add water
###############
$file = "$name.water.gro";
if (!  -e  $file) {
    $program = "$gromacs_path/genbox";
    $log = "genbox.log";
    print "\t runnning $program \n";
    $command = "$program -cp $name.box.gro -o $file  -cs spc216.gro -p $name.top > $log  2>&1";
    # -c is the centering command
   ( system $command) 
	&& die "Error:\n$command\nerror"; 
    ( $cmd_log ) &&  print CMD_LOG "$command\n";
    error_state_check ( $log, ("masses will be determined based on residue and atom names",
			"radii will be determined"));


    # I need to fix the number of waters now:
    if ( defined $multiplicity{"water"} ) {
	my ( $line, $totsol, @aux2);
	$ret = `grep SOL $name.top`;


	@aux = split '\n', $ret;
	$totsol = 0;
	foreach $line ( @aux ) {
	    @aux2 = split " ", $line;
	    $totsol += $aux2[1];
	} 
	`grep -v SOL  $name.top > tmp`;
	`mv tmp $name.top`; 
	`echo \"SOL    $totsol\" >> $name.top`
    }
   
} else { 
    print "\t $file found\n";  
}




###############
# geometry "minimization"
###############
# preprocessing
$file = "$name.em_input.tpr";
if (!  -e  $file) {

    my $ion_name;
    $program = "$grompp";
    $log     = "grompp.log";
    print "\t runnning $program \n";
    $command = "$program -f em.mdp -c $name.water.gro -p $name.top -o $file> $log  2>&1";
   
    # grompping
    # -c is the centering command
    ( $cmd_log ) &&  print CMD_LOG "$command\n";
    system $command 
	|| die "Error:\n$command\nerror"; 
    $charge = 0;
    $ret = `grep \'System has non-zero total charge\' $log`;
    if ( $ret ) {
	@aux = split " ", $ret;
	$charge = int ( sprintf "%5.0f", pop @aux ) ;
	print "\t charge $charge\n"; 
       ( $charge ) && print "\t system has nonzero charge -  adding counterions\n";
    }    
    if ( $charge ) {
	error_state_check ( "$log", ("maxwarn", "WARNING 1", 
				    "There was 1 warning", 
				    "defaults to zero instead of generating an error"));
    } else {
	error_state_check ( "$log", ("maxwarn","defaults to zero instead of generating an error"));
    }

    if ( $charge ) { # adding counterions

	$log = "grompp.counterion.log";
	my  $ion_request;

	if ( $np) {
	    # redo grompp to make it serial
	    $program = "$grompp_root";
	    print "\t re-runnning $program \n";
	    $command = "$program -f em.mdp -c $name.water.gro -p $name.top -o $file > $log  2>&1";
	    # grompping
	    ( $cmd_log ) &&  print CMD_LOG "$command\n";
	    (system $command)
		&& die "Error:\n$command\nerror"; 
	}

	#genion
	$program = "$gromacs_path/genion";
	print "\t\t runnning $program \n";
	if ($charge > 0) {
	    $ion_name = $neg_ion;
	    $ion_request = " -nname $ion_name -nn $charge ";
	} else {
	    $charge = -$charge; # notice we turn the $charge variable into its absolute value
	    $ion_name = $pos_ion;
	    $ion_request =  " -pname $ion_name -np $charge ";
	}

	#   figure out what is the numer of "Water" group
	$ret = `echo 0 | $gromacs_path/genion -s $name.em_input.tpr -o $name.ion.gro -pname Na -np 13   2> groups`;
	$ret = `grep Water groups | grep -v non`;

	@aux = split  " ", $ret;
	my $water_group_number = $aux[1];

	$command = "echo $water_group_number | $program -s $file -o $name.ion.gro $ion_request> $log  2>&1";
	( $cmd_log ) &&  print CMD_LOG "$command\n";
	(system $command)  
	    && die "Error:\n$command\nerror";  

	@aux = split " ", `grep Group $log | grep SOL`; 
	if ( $aux[1] != $genion_solvent_code ) {
	    $genion_solvent_code = $aux[1];
	    $command = "echo $genion_solvent_code | $program -s $file -o $name.ion.gro $ion_request> $log  2>&1";
	    ( $cmd_log ) &&  print CMD_LOG "$command\n";
	    (system $command)  
		&& die "Error:\n$command\nerror";  
	}
	( -e  "$name.ion.gro" ) || die "Error:\n$command\nerror";  
	error_state_check ( "$log", ("turning of free energy, will use lambda=0"));
	print "\t\t added $charge $ion_name \n";

	# "fix" the topology file
	if ( 1 ) {
	    my $number_of_waters;
	    my $new_solvent_line;

	    $file = "$name.top";
	    open ( OLD_FILE, "<$file") ||
		die "error opening $file: $!.";
	    open ( NEW_FILE, ">$file.new") ||
		die "error opening $file.new: $!.";
	    while ( <OLD_FILE> ){
		if ( /^SOL/ ) {
		    @aux = split " ", $_;
		    $number_of_waters = pop @aux;
		    push @aux, $number_of_waters - $charge;
		    $new_solvent_line = join "  ", @aux;
		    print NEW_FILE $new_solvent_line."\n";
		    print NEW_FILE "$ion_name     $charge\n";
		} else {
		    print NEW_FILE;
		}
	    }
	    close NEW_FILE;
	    close OLD_FILE;
	    `mv $file $file.old`;
	    `mv $file.new $file`;
	}


	$input_system_for_md = "$name.ion.gro";

	#grompp again
	$program = "$grompp";
	print "\t runnning $program \n";
	$command = "$program -f em.mdp -c $input_system_for_md  -p $name.top -o $name.em_input.tpr > $log  2>&1";
	(system $command) 
	    && die "Error:\n$command\nerror"; 
	( $cmd_log ) &&  print CMD_LOG "$command\n";
	error_state_check ( "$log", ("maxwarn", 
				    "defaults to zero instead of generating an error"));

    } else {
	$input_system_for_md = "$name.water.gro";
    }
    
} else { 
    print "\t $file found\n"; 
    if ( -e  "$name.ion.gro" ) {
	$input_system_for_md = "$name.ion.gro";
    } elsif  ( -e  "$name.water.gro" ) {
	$input_system_for_md = "$name.water.gro";
    } else {
	die "Error: neither $name.ion.gro nor $name.water.gro found.\n";
    }
} 




# minimization
$file = "$name.em.edr";
if (!  -e  $file) {
    $program = "$gromacs_run";
    $log = "energy_minimization.log";
    print "\t runnning $program -- energy minimization \n";
    $command = "$program -s $name.em_input.tpr  -o $name.em_output.trr  -c $input_system_for_md  -e $file> $log  2>&1";
    # -c is the centering command
    (system $command)
	&& die "Error:\n$command\nerror"; 
    ( $cmd_log ) &&  print CMD_LOG "$command\n";
    error_state_check ( "$log", ("masses will be determined based on residue and atom names"));
    $ret = `grep \'Steepest Descents converged\' $log`;
    if ( ! $ret ) {
	print "$program did not converge. Please check the file called $log.\n";
	#exit (1);
    }
} else {
    print "\t $file found\n"; 
}

$minimization && exit;


###############
# position-restrained md
###############
# preprocessing
$file = "$name.pr.tpr";
if (!  -e  $file) {
    $program = "$grompp";
    $log  = "grompp_before_pos_restr_md.log";
    print "\t runnning $program \n";
    $command = "$program -f pr.mdp -c $input_system_for_md  -r $input_system_for_md  -p $name.top -o $file -maxwarn 2  > $log  2>&1";
    ( $cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command )
	&& die "Error:\n$command\n"; 
    error_state_check ( "$log", ("maxwarn", 
				    "defaults to zero instead of generating an error"));
} else {
    print "\t $file found\n"; 
}


# md 
$file = "$name.pr.edr";
if (!  -e  $file) {
    $program = "$gromacs_run";
    $log     = "pos_restr_md.log";
    print "\t runnning $program -- position restrained md \n";
    $command = "$program -s $name.pr.tpr -o $name.pr.trr -c $input_system_for_md  -e $file > $log  2>&1";
    # 
    ( $cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command)
	&& die "Error:\n$command\nerror"; 
    error_state_check ( "$log", ("maxwarn"));
} else {
    print "\t $file found\n"; 
}



###############
# ... and finally ... the MD SIMULATION!
###############
# preprocessing
###############
$file = "$name.md.tpr";
if (!  -e  $file) {
    $program = "$grompp";
    $log  = "grompp_before_production_run.log";
    print "\t runnning $program \n";
    $command = "$program -f md.mdp -c $input_system_for_md  -r $input_system_for_md  -p $name.top -o $file -maxwarn 2 > $log  2>&1";
    # 
    ( $cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command) 
	&& die "Error:\n$command\n"; 
    error_state_check ( "$log", ("maxwarn", 
				    "defaults to zero instead of generating an error"));
} else {
    print "\t $file found\n"; 
}


# if we are turning this into a replica exchange sim,
# we exit here - and let some other script take it from here

($remd) && exit;


# md 
$file = "$name.md.edr";
if (!  -e  $file) {
    $program = "$gromacs_run";
    $log  = "production_run.log";
    print "\t runnning $program --  md  simulation proper\n";
    $command = "$program -s $name.md.tpr -o $name.md.trr -c $input_system_for_md ".
	"  -e $file > $log  2>&1";
    #
    ( $cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command) 
	&& die "Error:\n$command\n"; 
   error_state_check ( "$log", ("maxwarn"));
} else {
    print "\t $file found\n"; 
}


###############
#  compress the trajectory, strip hydrogens and gzip
###############

=pod
$command = "$gromacs_path/trjconv -f  $name.md.trr -o $name.md.xtc";
system $command 
	|| die "Error:\n$command\nerror";
`rm $name.md.trr`;


$command = "echo 0 | $gromacs_path/trjconv -s $name.md.tpr -f $name.md.xtc -o tmp.pdb ";
system $command 
	|| die "Error:\n$command\nerror";

$command = "  awk -F\ '\' \'\$14 != \"H\"\' tmp.pdb >  trajectory.pdb";
system $command 
	|| die "Error:\n$command\nerror";

if ($np) {
    $command = " gzip  trajectory.pdb";
    system $command 
	|| die "Error:\n$command\nerror";
}
`rm tmp.pdb`;
=cut



#################################################################################

###############
# extending the run
###############

#tpbconv -f old_name.trr -s old_name.tpr -e old_name.edr -o new_name.tpr -until <time in ps>
#mdrun -v -deffnm new_name -s new_name.tpr


#################################################################################
#################################################################################
sub  error_state_check (@ ) {

    my $logfile = shift @_;
    my @tolerated_warnings = @_;
    my $ret;
    my @lines;
    my $serious;
    my $line;

    
    $ret = `grep -i error $logfile | grep -v gcq`; # the retarded quote may have the word error in it
    if ( $ret ) {
	$serious = "";
	@lines = split "\n", $ret;
	foreach $line ( @lines ) {
	    if ( grep  {$line =~ $_ } @tolerated_warnings) {
	    } else {
		$serious = $line;
		last;
	    }
	}
	if ( $serious ) {
 	    print "$program ended in error state. Please check the file called $logfile.\n";
	    print "$serious\n\n";
	    exit (1);
	}
    }

    $ret = `grep -i warning $logfile`;
 
    if ( $ret ) {
	$serious = "";
	@lines = split "\n", $ret;
	foreach $line ( @lines ) {
	    if ( grep  {$line =~ $_ } @tolerated_warnings) {
	    } else {
		$serious = $line;
		last;
	    }
	}
	if ( $serious ) {
 	    #print "$program issued a warning. Please check the file called $logfile.\n";
	    #print "$serious\n\n";
	    #exit (1);
	}
   }
 


}



sub make_fake_top (@){

    my $name = shift @_;

    open (FAKE, ">$name.top") ||
	die "Cno $name.top:$!\n";
    
    print FAKE "; topology wrapper written by $0.\n";
    print FAKE "; Include forcefield parameters\n";
    print FAKE "#include  \"$forcefield.ff/forcefield.itp\n";

    print FAKE "; Include Position restraint file\n";
    print FAKE "#ifdef POSRES\n";
    print FAKE "#include \"posre.itp\"\n";
    print FAKE "#endif\n";

    print FAKE "; Include water topology\n";
    print FAKE "#include \"$forcefield.ff/spc.itp\"\n";

    print FAKE "#ifdef POSRES_WATER\n";
    print FAKE "; Position restraint for each water oxygen\n";
    print FAKE "[ position_restraints ]\n";
    print FAKE ";  i funct       fcx        fcy        fcz\n";
    print FAKE "   1    1       1000       1000       1000\n";
    print FAKE "#endif\n";

    print FAKE "; Include topology for ions\n";
    print FAKE "#include \"$forcefield.ff/ions.itp\"\n";

    print FAKE "[ system ]\n";
    print FAKE "; Name\n";
    print FAKE "Protein in water\n";

    print FAKE "[ molecules ]\n";
    print FAKE "; Compound        #mols\n";

    close FAKE;

}
