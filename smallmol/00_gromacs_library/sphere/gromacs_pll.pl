#! /usr/bin/perl -w

use strict;
sub  error_state_check (@ ); # defined at the bottom

defined $ARGV[0] ||
    die "Usage: gromacs.pl <pdb name root> [<list of ligand names>].\n";

#`unlimit`;

##############################################################################
##############################################################################
my $sphere = 1;
my $np = 4;

my $pdb2gro;
my $prodrg = "/home/i/imihalek/downloads/prodrg/prodrg";
my $prodrg_param = "/home/i/imihalek/downloads/prodrg/prodrg.param";
my $gromacs_path = "/home/i/imihalek/downloads/gromacs_old/i686-pc-linux-gnu/bin";
#my $gromacs_path = "/home/i/imihalek/downloads/gromacs/bin";
my $mpirun = "/usr/local/bin/mpirun";
my ($mdrun, $grompp);
my $gromacs_run;
my $groc;
my $grompp_root;
if ( $sphere ) {
    $gromacs_path =  "/usr/local/bin/";
    $mdrun = "/usr/local/bin/mdrun";
    if ( $np > 1 ) {
	$gromacs_run = "$mpirun  -machinefile \$PBS_NODEFILE -np $np $mdrun -np $np";
    } else {
	$gromacs_run = $mdrun;
    }
    $grompp_root = "/home/imihalek/gromacs/gromacs-3.2.1/src/kernel/grompp";
    if ( $np > 1 ) {
	$grompp = "$grompp_root  -np $np";
    } else {
	$grompp = $grompp_root;
    }
    $groc = "/home/imihalek/perlscr/gro_concat.pl"


} else {
    $gromacs_run = "mdrun";
    $grompp = "grompp";
    $groc = "/home/i/imihalek/perlscr/gro_concat.pl";

}


##############################################################################
##############################################################################

my $name        = $ARGV[0]; # expect ok pdb file called $name.pdb present
my $forcefield  = "gmx";
my $box_type    = "cubic";
my $box_edge    =  0.7; # distance of the box edge from the molecule in nm 
my $neg_ion     = "Cl"; # theses names depend on the choice of the forcefield (check "ions.itp")
my $pos_ion     = "Na";
my $genion_solvent_code = 12;

##############################################################################
##############################################################################
## in  ~/moldyn/gromacs/tutorials/generic_mdps
foreach ( "em.mdp", "pr.mdp", "md.mdp" ) {
    ( -e $_ ) || die "\n$_ not found.\n\n";
}


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
    close LF;
    foreach $ligand_name ( @ligand_names ) {
	print $ligand_name, "\n";
	if ( $ligand_name  eq "water" ) {
	    $file = "water.gro"; #all waters are understood to be in the same file; itp is standard for water
	    ( -e $file) || die "$file not found in ".`pwd`;
	} else {
	    for ($mult=1; $mult <= $multiplicity{$ligand_name}; $mult ++ ) {
		my $ln = (lc $ligand_name);
		my $name_root = $ln.$mult;
		my $cmd;
		$file = $name_root.".gro";
		next if ( -e $file);
		( -e "$name_root.pdb" ) || die "Cno $name_root.pdb: $! .\n";
		if ( defined $ion{$ligand_name} ) {
		    $cmd = "$pdb2gro <  $name_root.pdb >  $name_root.gro ";
		    (system $cmd) && die "Error running $cmd\n";
		} else {
		    `mkdir tmp_prodrg`;
		    chdir "tmp_prodrg";
		    $cmd = "$prodrg  ../$name_root.pdb $prodrg_param CGRP ";
		    (system $cmd) && die "Error running $cmd\n";
		    `mv DRGFIN.GRO ../$name_root.gro`;
		    `mv DRGGMX.ITP ../$ln.itp`; # no mult!
		    chdir "..";
		    `rm -rf tmp_prodrg`;
		}
	    }
	    if ( ! defined $ion{$ligand_name} ) { # ion "topology" is included by default
		$file = lc $ligand_name.".itp"; # all other molecules have their own gro and itp files
		( -e $file) || die "$file not found in ".`pwd`;
	    }
	}
    } 
} 

#exit;
############################################################################## 
##############################################################################
my $command;
my $ret;
my $program;
my @aux;
my $charge;
my $input_system_for_md;

###############
# process pdb into gro and topology files
###############
if ( ! -e "$name.top" || ! -e "$name.gro" )  {

    $program = "$gromacs_path/pdb2gmx";
    print "\t runnning $program \n";
    $command = "$program  -ignh -ff $forcefield -f $name.pdb -o $name.gro -p $name.top >& log";
    # -ignh instructs pdb2gmx to ingore H and place its own
    system $command 
	|| die "Error:\n$command\nerror"; # looks like it exits with something on success
    error_state_check ( "log", ());

    if ( @ligand_names ) {
	my $written;

	# concatenate gro files
	$command = "$groc  $name.gro ";
	foreach $ligand_name ( @ligand_names ) {
	    if ($ligand_name eq "water") {
		$command .=  $ligand_name.".gro ";
		next;
	    }
	    for ($mult=1; $mult <= $multiplicity{$ligand_name}; $mult ++ ) {
		$command .=  (lc $ligand_name).$mult.".gro ";
	    }
	}
	$command .= " > tmp";
	my $ret;
	$ret = system $command ;
	$ret	&& die "Error:\n$command\n$ret"; 
	`mv tmp $name.gro`;
	# change the topology (*.top) file
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
			print OF  $ligand_name."     $multiplicity{$ligand_name}\n";
		    }
		}
		print OF "\n";
	    }
	}
	close OF;
	close TOP;
	`mv tmp $name.top`;
	# fix the pr.mdp and md.mdp.files
	foreach $file ( "pr.mdp", "md.mdp" ) {

	     # temp coupling not defined for ions
	    open ( OLD_FILE, "<$file") ||
		die "error opening $file: $!.";
	    open ( NEW_FILE, ">$file.new") ||
		die "error opening $file.new: $!.";
	    while ( <OLD_FILE> ){
		if ( /^tc_grps/ ) {
		    @aux = split " ", $_;
		    foreach $ligand_name ( @ligand_names ) {
			next if  ($ligand_name eq "water");
			push @aux, $ligand_name; 
		    }
		    print NEW_FILE (join "  ", @aux), "\n";
		    
		} elsif ( /^tau_t/ || /^ref_t/) {
		    my $tau_t;
		    @aux = split " ", $_;
		    $tau_t = pop @aux;
		    push @aux,  $tau_t; 
		    foreach $ligand_name ( @ligand_names ) {
			next if  ($ligand_name eq "water");
			push @aux,  $tau_t; 
		    }
		    print NEW_FILE (join "  ", @aux), "\n";
		} else {
		    print NEW_FILE;
		}
	    }
	   
	    close NEW_FILE;
	    close OLD_FILE;
	    `mv $file $file.old`;
	    `mv $file.new $file`;
	    
	}
	
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
    print "\t runnning $program \n";
    $command = "$program -f $name.gro -o $file  -bt $box_type -d $box_edge -c >& log";
    # -c is the centering command
    system $command 
	|| die "Error:\n$command\nerror"; 
    error_state_check ( "log", ("masses will be determined based on residue and atom names"));
} else {
    print "\t $file found\n"; 
}


###############
# add water
###############
$file = "$name.water.gro";
if (!  -e  $file) {
    $program = "$gromacs_path/genbox";
    print "\t runnning $program \n";
    $command = "$program -cp $name.box.gro -o $file  -cs spc216.gro -p $name.top >& log";
    # -c is the centering command
    system $command 
	|| die "Error:\n$command\nerror"; 
    error_state_check ( "log", ("masses will be determined based on residue and atom names"));

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
    print "\t runnning $program \n";
    $command = "$program -f em.mdp -c $name.water.gro -p $name.top -o $file >& log";
    print "$command\n", 
    # grompping
    # -c is the centering command
    system $command 
	|| die "Error:\n$command\nerror"; 
    $charge = 0;
    $ret = `grep \'System has non-zero total charge\' log`;
    if ( $ret ) {
	@aux = split " ", $ret;
	$charge = int ( sprintf "%5.0f", pop @aux ) ;
	print "\t charge $charge\n"; 
       ( $charge ) && print "\t system has nonzero charge -  adding counterions\n";
    }    
    if ( $charge ) {
	error_state_check ( "log", ("maxwarn", "WARNING 1", "There was 1 warning"));
    } else {
	error_state_check ( "log", ("maxwarn"));
    }
    if ( $charge ) { # adding counterions

	my  $ion_request;

	if ( $np > 1) {
	    # redo grompp to make it serial
	    $program = "$grompp_root";
	    print "\t re-runnning $program \n";
	    $command = "$program -f em.mdp -c $name.water.gro -p $name.top -o $file  >& log";
	    # grompping
	    # -c is the centering command
	    system $command 
		|| die "Error:\n$command\nerror"; 
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
	$genion_solvent_code += @ligand_names;
	$command = "echo $genion_solvent_code | $program -s $file -o $name.ion.gro $ion_request >& log";
	system $command  
	    || die "Error:\n$command\nerror";  
	@aux = split " ", `grep Group log | grep SOL`; 

	print"********* $command\n";
	print"********* @aux\n";
	if ( $aux[1] != $genion_solvent_code ) {
	    $genion_solvent_code = $aux[1];
	    $command = "echo $genion_solvent_code | $program -s $file -o $name.ion.gro $ion_request >& log";
	    system $command  
		|| die "Error:\n$command\nerror";  
	}
	( -e  "$name.ion.gro" ) || die "Error:\n$command\nerror";  
	error_state_check ( "log", ("turning of free energy, will use lambda=0"));
	print "\t\t added $charge $ion_name \n";

	# "fix" the toplogy file
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

	# fix the pr.mdp and md.mdp.files
	foreach $file ( "pr.mdp", "md.mdp" ) {

	    $ret = `grep tc_grps $file | grep $ion_name `;
	    next if ( $ret);
	     # temp coupling not defined for ions
	    open ( OLD_FILE, "<$file") ||
		die "error opening $file: $!.";
	    open ( NEW_FILE, ">$file.new") ||
		die "error opening $file.new: $!.";
	    while ( <OLD_FILE> ){
		if ( /^tc_grps/ ) {
		    @aux = split " ", $_;
		    push @aux,  $ion_name; 
		    print NEW_FILE (join "  ", @aux), "\n";
		    
		} elsif ( /^tau_t/ || /^ref_t/) {
		    my $tau_t;
		    @aux = split " ", $_;
		    $tau_t = pop @aux;
		    push @aux,  $tau_t; 
		    push @aux,  $tau_t; 
		    print NEW_FILE (join "  ", @aux), "\n";
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
	$command = "$program -f em.mdp -c $input_system_for_md  -p $name.top -o $name.em_input.tpr >& log";
	system $command 
	    || die "Error:\n$command\nerror"; 
	error_state_check ( "log", ("maxwarn"));

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
    print "\t runnning $program -- energy minimization \n";
    $command = "$program -s $name.em_input.tpr  -o $name.em_output.trr  -c $input_system_for_md  -e $file >& log";
    # -c is the centering command
    system $command 
	|| die "Error:\n$command\nerror"; 
    error_state_check ( "log", ("masses will be determined based on residue and atom names"));
    $ret = `grep \'Steepest Descents converged\' log`;
    if ( ! $ret ) {
	print "$program did not converge. Please check the file called \"log\".\n";
	#exit (1);
    }
} else {
    print "\t $file found\n"; 
}

###############
# position-restrained md
###############
# preprocessing
$file = "$name.pr.tpr";
if (!  -e  $file) {
    $program = "$grompp";
    print "\t runnning $program \n";
    $command = "$program -f pr.mdp -c $input_system_for_md  -r $input_system_for_md  -p $name.top -o $file  >& log";
    # -c is the centering command
    system $command 
	|| die "Error:\n$command\nerror"; 
    error_state_check ( "log", ("maxwarn"));
} else {
    print "\t $file found\n"; 
}
# md 
$file = "$name.pr.edr";
if (!  -e  $file) {
    $program = "$gromacs_run";
    print "\t runnning $program -- position restrained md \n";
    $command = "$program -s $name.pr.tpr -o $name.pr.trr -c $input_system_for_md  -e $file  >& log";
    # -c is the centering command
    system $command 
	|| die "Error:\n$command\nerror"; 
    error_state_check ( "log", ("maxwarn"));
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
    print "\t runnning $program \n";
    $command = "$program -f md.mdp -c $input_system_for_md  -r $input_system_for_md  -p $name.top -o $file  >& log";
    # -c is the centering command
    system $command 
	|| die "Error:\n$command\nerror"; 
    error_state_check ( "log", ("maxwarn"));
} else {
    print "\t $file found\n"; 
}
# md 
$file = "$name.md.edr";
if (!  -e  $file) {
    $program = "$gromacs_run";
    print "\t runnning $program --  md  simulation proper\n";
    $command = "$program -s $name.md.tpr -o $name.md.trr -c $input_system_for_md  -e $file  >& log";
    # -c is the centering command
    system $command 
	|| die "Error:\n$command\nerror"; 
    error_state_check ( "log", ("maxwarn"));
} else {
    print "\t $file found\n"; 
}


###############
#  compress the trajectory, strip hydrogens and gzip
###############

$command = "$gromacs_path/trjconv -f  $name.md.trr -o $name.md.xtc";
system $command 
	|| die "Error:\n$command\nerror";

=pod
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

`rm $name.md.trr`;


#################################################################################

###############
# extending the run
###############

#tpbconv -f old_name.trr -s old_name.tpr -e old_name.edr -o new_name.tpr -until <time in ps>
#mdrun -v -deffnmnew_name -s new_name.tpr


#################################################################################
#################################################################################
sub  error_state_check (@ ) {

    my $logfile = shift @_;
    my @tolerated_warnings = @_;
    my $ret;

    
    $ret = `grep -i error $logfile | grep -v gcq`; # the retarded quote may have the word error in it
    if ( $ret ) {
	print "$program ended in error state. Please check the file called \"log\".\n";
	print "$ret\n\n";
	exit (1);
    }
    $ret = `grep -i warning $logfile`;
    if ( $ret ) {
	if ( grep  {$ret =~ $_ } @tolerated_warnings) {
	} else {
	    print "$program issued a warning. Please check the file called \"log\".\n";
	    print "$ret\n\n";
	    exit (1);
	}
    }

}
