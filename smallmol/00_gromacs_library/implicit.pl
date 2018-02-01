#! /usr/bin/perl -w

# NOTE: water is hardcoded as tip3p here

use strict;
sub  error_state_check (@ ); # defined at the bottom
sub make_fake_top (@);
sub fix_impl_params (@);

defined $ARGV[0] ||
    die "Usage: gromacs.pl <pdb name root>/no_protein [<list of ligand names> -postp|-min].\n";

$| = 1; # turn on autoflush;

# remd flag must be the last, for now
my $postp = 0;
($ARGV[$#ARGV] eq "-postp") && ($postp = 1);
my $minimization = 0;
($ARGV[$#ARGV] eq "-min" ) && ($minimization = 1);

##############################################################################
##############################################################################
my $np = 0;

my $write_cmd_log = 1;
my $pdb2gro;
my $gromacs_path = "/usr/local/bin";
my $perl_path    = "/home/ivanam/perlscr/gromacs";
#my $perl_path   = "/Users/ivana/perlscr/gromacs";
my $mpirun       = "/usr/local/bin/mpirun";
my ($mdrun, $mdrun_mpi, $grompp);
my $gromacs_run;
my $groc;
my $genrestr; 
my $trjconv;
my $g_energy;
my $babel;
my $gmx_wrapper;
my $grompp_root;
my $protein = 0;
my $home = `pwd`; chomp $home;

my $hostname = `hostname`;
chomp $hostname;

  
$gromacs_run = "$gromacs_path/mdrun";
$grompp      = "$gromacs_path/grompp";
$groc        = "$perl_path/gro_concat.pl";
$pdb2gro     = "$perl_path/pdb2gro.pl";
$gmx_wrapper = "$perl_path/gromacs.pl";
$genrestr    = "$gromacs_path/genrestr";
$trjconv     = "$gromacs_path/trjconv";
$g_energy     = "$gromacs_path/g_energy";



if ($hostname eq "donkey") {

    $gromacs_path = "/usr/bin";
    $grompp       = "/usr/bin/grompp";
    $gromacs_run  = "/usr/bin/mdrun";
    $mdrun_mpi    = "/usr/bin/mdrun_mpi";
    $mpirun       = "/usr/bin/mpirun";
    $genrestr    = "/usr/bin/genrestr";
    $trjconv     = "/usr/bin/trjconv";

} else {
    $mdrun_mpi = "";
}

foreach ($gromacs_path,  $perl_path, $gromacs_run, $grompp, $groc, $pdb2gro, $genrestr) {
    ( -e $_ ) || die "\n$_ not found.\n\n";
}

$gromacs_run .= " -nt 1";


if ( $postp) {
    foreach ($trjconv) {
	( -e $_ ) || die "\n$_ not found.\n\n";
    }
}

if ( $np && ! -e $mpirun) {
    die "\n$mpirun not found.\n\n";
}


##############################################################################
##############################################################################

my $name        = $ARGV[0]; # expect ok pdb file called $name.pdb present
my $forcefield  = "amber99sb";
my $water = "tip3p";
#my $forcefield = "oplsaa";
#my $forcefield = "gmx";
#my $box_type   = "cubic";
my $box_type    = "triclinic";
my $box_edge    =  0.7; # distance of the box edge from the molecule in nm 
my $neg_ion     = "Cl"; # theses names depend on the choice of the forcefield (check "ions.itp")
my $pos_ion     = "Na";
my $genion_solvent_code = 12;

if (  $forcefield  eq  "oplsaa" || $forcefield  eq "amber99sb" ) {
    $neg_ion  = "CL";
    $pos_ion  = "NA";
}

##############################################################################
##############################################################################
my ($in_dir, $top_dir, $em1_dir, $em2_dir,  $production, $postpcss) = 
    ("00_input", "01_topology", "02_em_steepest", "03_em_lbfgs", "04_production", "05_postpcss");

foreach ( $in_dir, "$in_dir/em_steep.mdp", "$in_dir/em_lbfgs.mdp",  "$in_dir/md.mdp" ) {
    ( -e $_ ) || die "\n$_ not found.\n\n";
}

if ( $postp) {
    foreach ("$in_dir/pr.mdp") {
	( -e $_ ) || die "\n$_ not found.\n\n";
    }
}

if ( $name eq "no_protein") {
    $protein = 0;
    (  defined $ARGV[1] ) ||
	die "No protein, no small molecule ... what are we doing here?\n";
    $box_edge    =  2.0;
} else {
    ( -e "$in_dir/$name.pdb") || die "$in_dir/$name.pdb not found\n";
    $protein = 1;
}

#foreach ( $top_dir, $em1_dir, $em2_dir, $pr1_dir, $pr2_dir, $production) {
#    (-e $_) || `mkdir $_`;
#}

##############################################################################
##############################################################################
$write_cmd_log &&  ( open (CMD_LOG, ">cmd.log") || die "Cno cmd.log: $!.\n");


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

(-e  $top_dir) || `mkdir $top_dir`;

if ( defined $ARGV[1] ) {
    my @aux;
    open ( LF, "<$in_dir/$ARGV[1]") ||
	die "Cno $in_dir/$ARGV[1]: $!.\n";
    while ( <LF> ) {
	next if (! /\S/);
	chomp;
	@aux = split;
	push @ligand_names, $aux[0];

	if ( $aux[0] ne "water" ) {
	    if ( defined  $aux[1] ) {
		$multiplicity{$aux[0]} = $aux[1];
	    } else {
		die "Number of instances of molecule $aux[0] not defined in $ARGV[1].\n";
	    }
	}
       
    }
    $protein || ($name = join "_", @ligand_names);

    close LF;
    foreach $ligand_name ( @ligand_names ) {

	for ($mult=1; $mult <= $multiplicity{$ligand_name}; $mult ++ ) {
	    my $ln = (lc $ligand_name);
	    my $name_root = $ln;
	    my $cmd;
		

	    if ( $multiplicity{$ligand_name} > 1 ) {
		$name_root = $ln.$mult;
	    }

	    $file = $name_root.".gro";
	    next if ( -e "$in_dir/$file");
	    ( -e "$in_dir/$name_root.pdb" ) || die "Cno $in_dir/$name_root.pdb: $!.\n";

	    # handle ions
	    if ( defined $ion{$ligand_name} ) {
		if (  $forcefield  eq  "oplsaa" ) {
		    my $tmp = "tmp.ion";
		    my ($ion_type, $new_ion_type, $new_line);
		    my ($serial, $res_seq, $x, $y, $z);
		    open (ION_IF, "<$in_dir/$name_root.pdb");
		    open (ION_OF, ">$top_dir/$name_root.gro");
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
		    $cmd = "$pdb2gro <  $in_dir/$name_root.pdb >  $top_dir/$name_root.gro ";
		    (system $cmd) && die "Error running $cmd\n";
		    ( $write_cmd_log ) &&  print CMD_LOG "$cmd\n";
		}

		# handle other types of ligands
	    } else {

		( -e "$in_dir/$name_root.gro") || die "$in_dir/$name_root.gro must be prepared in advance\n";
		( -e "$in_dir/$ln.itp") || die "$in_dir/$ln.itp must be prepared in advance\n";
		# check whether we have the parameters for GB, if not xi the file
		my $ret = "" || `grep implicit_genborn_params $in_dir/$ln.itp`;
		$ret || fix_impl_params ("$in_dir/$ln.itp");
	    }
	}

	if ( ! defined $ion{$ligand_name} ) { # ion "topology" is included by default
	    # all other molecules have their own gro and itp files
	    $file = lc $ligand_name.".itp"; 
	    ( -e "$in_dir/$file") || die "$in_dir/$file not found in ".`pwd`;

	    # check whether we have the parameters for GB, if not xi the file
	    my $ret = "" || `grep implicit_genborn_params $in_dir/$file`;
	    $ret || fix_impl_params ("$in_dir/$file");
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

chdir $home;
(-e  $top_dir) || `mkdir $top_dir`;
chdir $top_dir;

if ( ! -e "$name.top" || ! -e "$name.gro" )  {

    if ( $protein ) {
	$program = "$gromacs_path/pdb2gmx";
	$log = "pdb2gmx.log";
	print "\t runnning $program \n";
	# -ignh instructs pdb2gmx to ingore H and place its own
	$command  = "$program  -ignh -ff $forcefield -water none  -f ../$in_dir/$name.pdb -o $name.gro -p $name.top ";
	( -e "pdb2gmx_in") && ( `rm pdb2gmx_in`); 
	if ( -e "ssbridges" ) {
	    # ssbridges has one "y\n" for each ssbridge (sequentailly
	    # by the first cysteine) that we want to maintain
	    # I guess "n\n" for the ones we do not want too
	    `cat ssbridges > pdb2gmx_in`;
	    $command  .= " -ss";
	} 

	if ( -e "termini" ) {
	    # termini has the following format: "2\n2\n" for each chain 
	    # (4 x 2 for two chains and so on) 
	    `cat termini >> pdb2gmx_in`;
	    $command  .= " -ter";
	}
	
	( -e "pdb2gmx_in") && ($command  .= " < pdb2gmx_in");
	$command  .= " > $log 2>&1";
	($write_cmd_log)  &&  print CMD_LOG "$command\n";

	(system $command)
	    && die "Error:\n$command\nerror"; # looks like it exits with something on success
	error_state_check ( $log, ("masses will be determined based on residue and atom names"));

	( -e "pdb2gmx_in") && `rm pdb2gmx_in`;
    }

    if ( @ligand_names ) {
	my $written;
	`cp ../$in_dir/*gro ../$in_dir/*.itp .`; # too messy otherwise
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
	    ( $write_cmd_log ) &&  print CMD_LOG "$command\n";
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
		    print OF "#include \"$home/$top_dir/". lc $ligand_name.".itp\"\n";
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
# place the sytem in a box - we don't need it but I'm too lazy to change the script
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
    ( $write_cmd_log ) &&  print CMD_LOG "$command\n";
    # some bug in editconf: if the box triclinic it reports "no boxtype spcified,
    # but constructs something which has all three unit cell vectors perpendicular 
    # and of different lengths
    error_state_check ( $log, ("masses will be determined based on residue and atom names",
				    "No boxtype specified"));
} else {
    print "\t $file found\n"; 
}

############################################################################## 
##############################################################################


###############
# geometry "minimization"
###############
# preprocessing  = making of  tpr file

chdir $home;
(-e $em1_dir) || `mkdir $em1_dir`;
chdir $em1_dir;

$file = "$name.em_input.tpr";
if (!  -e  $file) {

    # grompping
    # -c is the centering command
    my $ion_name;
    $program = "$grompp";
    $log     = "grompp.log";
    print "\t runnning $program \n";
    $command = "$program -f ../$in_dir/em_steep.mdp -c ../$top_dir/$name.box.gro ".
	"-p ../$top_dir/$name.top -o $file ";
    ( -e "../$in_dir/groups.ndx" ) && ($command .=  " -n ../$in_dir/groups.ndx ");
    $command 	.= "> $log  2>&1";
    ($write_cmd_log ) &&  print CMD_LOG "$command\n";
    system $command 
	|| die "Error:\n$command\nerror";  # if no charge, we are done


    $charge = 0;
    $ret = `grep \'System has non-zero total charge\' $log`;
    if ( $ret ) {
	@aux = split " ", $ret;
	$charge = int ( sprintf "%5.0f", pop @aux ) ;
	print "\t charge $charge\n"; 
       ( $charge ) && print "\t system has nonzero charge - should be taken care of by the impl solvent\n";
    }    
 
    $input_system_for_md = "$name.box.gro";
    
    
} else { 
    print "\t $file found\n"; 
} 


# minimization - round 1
$file = "$name.em_out.gro";
if (!  -e  $file) {
    #
    $program = "$gromacs_run";
    $log = "energy_minimization.log";
    print "\t runnning $program -- steepest descent energy minimization \n";
    $command = "$program -s $name.em_input.tpr -c $name.em_out.gro -o $name.em_out.trr > $log  2>&1";
    ( $write_cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command)
	&& die "Error:\n$command\nerror"; 
    error_state_check ( "$log", ("masses will be determined based on residue and atom names"));
    $ret = `grep \'Steepest Descents converged\' $log`;
    if ( ! $ret ) {
	print "$program did not converge. Please check the file called $log.\n";
	#exit (1);
    }
    #`rm  $name.em_out.trr`;
  
} else {
    print "\t $file found\n"; 
}



# minimization - round 2
# will have to dom some re-grompping, though
chdir $home;
(-e $em2_dir) || `mkdir $em2_dir`;
chdir $em2_dir;

$file = "$name.em_input.tpr";
if (! -e  $file) {

    #grompp again
    $log = "grompp.log";
    $program = "$grompp";
    print "\t runnning $program \n";
    $command = "$program -f ../$in_dir/em_lbfgs.mdp  -c ../$em1_dir/$name.em_out.gro  ".
	"  -p ../$top_dir/$name.top -o $name.em_input.tpr  -maxwarn 1 ";
    ( -e "../$in_dir/groups.ndx" ) && ($command .=  " -n ../$in_dir/groups.ndx ");
    $command 	.= "> $log  2>&1";
   ( $write_cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command) 
	&& die "Error:\n$command\nerror"; 
     error_state_check ("$log", 
			("maxwarn", "defaults to zero instead of generating an error"));
} else {
    print "\t $file found\n"; 
}

$file = "$name.em_out.gro";
if (!  -e  $file) {
    #
    $program = "$gromacs_run";
    $log = "energy_minimization.log";
    print "\t runnning $program -- lbfgs energy minimization \n";
    $command = "$program -s $name.em_input.tpr -c $name.em_out.gro -o $name.em_out.trr > $log  2>&1";
    ( $write_cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command)
	&& die "Error:\n$command\nerror"; 
    error_state_check ( "$log", ("masses will be determined based on residue and atom names"));
    $ret = `grep \'Low-Memory BFGS Minimizer converged\' $log`;
    if ( ! $ret ) {
	print "$program did not converge. Please check the file called $log.\n";
	#exit (1);
    }
    #`rm  $name.em_out.trr`;
} else {
    print "\t $file found\n"; 
}


$minimization && exit;


############################################################################## 
##############################################################################


###############
# ... and finally ... the production run!
###############

chdir $home;
(-e $production) || `mkdir $production`;
chdir $production;


# preprocessing
$file = "$name.md.tpr";
if (!  -e  $file) {


    if ( $protein ) {
	# there are some horror stories with the protein unfolding under implicit waters
	# so restrain the the backbone (posres.itp file):
	$program = "$genrestr";
	print "\t runnning $program \n";
	$command = "echo 4 | $program  -f ../$em2_dir/$name.em_out.gro > /dev/null 2>&1 ";
	( $write_cmd_log ) &&  print CMD_LOG "$command\n";
	(system $command) 
	    && die "Error:\n$command\n"; 
    }


    # grompp
    $program = "$grompp";
    $log  = "grompp_before_production_run.log";
    print "\t runnning $program \n";
    # -DPOSRES in mdp file will  instruct grompp to look for posres.itp file
    $command = "$program  -f  ../00_input/md.mdp -c ../$em2_dir/$name.em_out.gro  ".
	" -p ../$top_dir/$name.top -o $file -maxwarn 1 ";
    ( -e "../$in_dir/groups.ndx" ) && ($command .=  " -n ../$in_dir/groups.ndx ");
    $command 	.= "> $log  2>&1";
    ( $write_cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command) 
	&& die "Error:\n$command\n"; 
    error_state_check ( "$log", ("maxwarn", 
				    "defaults to zero instead of generating an error"));
} else {
    print "\t $file found\n"; 
}


# md 
$file = "$name.md.edr";
if (!  -e  $file) {
    $program = "$gromacs_run";
    $log  = "production_run.log";
    print "\t runnning $program --  md  simulation proper\n";
    $command = "$program -s $name.md.tpr -o $name.md.trr -c $name.md.gro ".
	"  -e $file > $log  2>&1";
    #
    ( $write_cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command) 
	&& die "Error:\n$command\n"; 
   error_state_check ( "$log", ("maxwarn"));
} else {
    print "\t $file found\n"; 
}

############################################################################## 
##############################################################################

$postp || exit;


###############
# postprocessing - specific for the implicit system
###############
sub postprocess_ligand (@_);
sub postprocess_protein ();
sub run_pos_restr ();
sub energies_from_posrestr();

chdir $home;
(-e $postpcss) || `mkdir $postpcss`;
chdir $postpcss;

#
# make pdb file
#
$file = "$name.md.pdb";
if (!  -e  $file) {
    $program = "$trjconv";
    $log     = "trjconv.log";
    print "\t runnning $program\n";
    $command = "(echo 0; echo 0) |  $program -s ../$production/$name.md.tpr ".
	"-f ../$production/$name.md.trr -o $name.md.pdb  -fit progressive  -skip 10  > $log  2>&1";
    #
    ( $write_cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command) 
	&& die "Error:\n$command\n"; 
   error_state_check ( "$log", ("maxwarn"));
} else {
    print "\t $file found\n"; 
}


#
# process individual frames
#
my $protein_pdb = "";
my $ligand_pdb  = "";
my @LIGAND_NAMES = map { uc $_} @ligand_names;
my $posres_in_dir = "00_input";
my $res_name;
my $current_ligand_name = "";
my $frame_ct = 0;
my $begin_time = 0;

(-e $posres_in_dir) && `rm -rf $posres_in_dir`;
`mkdir $posres_in_dir`;
`cp ../$in_dir/resolvation_input/*  $posres_in_dir`;

open (CUM, ">cumulative.table") || 
    die "Cno cumulative.table: $!.\n";
my   $avg_coul = 0;
my    $avg_lj  = 0;

open (IF, "<$file") || 
    die "Cno $file: $!.\n";

while ( <IF> ) {
 
    if  ($_=~/^ENDMDL/ || $_=~/^TER/)   {

	$frame_ct++;
	$begin_time = time;

	`rm -rf 0[1-4]*`;
	$ligand_pdb  && postprocess_ligand ($current_ligand_name);
	$protein_pdb && postprocess_protein();
	run_pos_restr ();
	energies_from_posrestr();

	printf "processed frame %d (%ds)\n", $frame_ct, time-$begin_time;



	$protein_pdb = "";
	$ligand_pdb  = "";
	$current_ligand_name = "";
	    
    } elsif ($_=~/^ATOM/) {

	$res_name = substr $_,  17, 3; $res_name=~ s/\s//g;

	if ( grep { /$res_name/} (@LIGAND_NAMES, "LIG") ) {
	    $current_ligand_name || ($current_ligand_name = lc $res_name);
	    $ligand_pdb .= $_;
	} else {
	    $protein_pdb .= $_;
	}

    } 
}
close IF;
cllose CUM;


sub energies_from_posrestr(){
    
    chdir "04_nvt_eq";

    my $log = "energy.log";
    my $ret;
    my $cmd = "(echo 44; echo 45) | $g_energy -f ener.edr  | tail -n50 | tee tmp";
    (system $cmd) && die "Error running $cmd.\n";

    my ($field_name, $avg, $err, $rmsd, $drift) = ();
    ($field_name, $avg, $err, $rmsd, $drift) = split " ", `grep \'LJ-SR\' tmp`;

    my $LJ = $avg;

    ($field_name, $avg, $err, $rmsd, $drift) = ();
    ($field_name, $avg, $err, $rmsd, $drift) = split " ", `grep \'Coul-SR\' tmp`;

    my $Coul = $avg;

    $avg_coul += $Coul;
    $avg_lj   += $LJ;

    printf CUM " %4d   %8.3f    %8.3f      %8.3f    %8.3f  \n", 
    $frame_ct, $LJ, $Coul,  $avg_coul/$frame_ct, $avg_lj/$frame_ct;
    
    chdir "..";

}


sub run_pos_restr () {

    my $log = "posres.log";

 
    $program = $gmx_wrapper;
    $command = "$program no_protein ligands -posres > $log 2>&1";
    ( $write_cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command) 
	&& die "Error:\n$command\n";    
 
}

sub postprocess_ligand (@_ ) {

    my $log = "lig_postp.log";
    my $name = $_[0];
    my $number_of_atoms = 0;
    my @ret;
    open ( PDB, ">ligand.pdb" ) ||
	die "Cno ligand.pdb: $!.\n";
    print PDB $ligand_pdb;
    close PDB;


    # make gro file
    $file = "$name.gro";
    $program = $pdb2gro ;
    $command = "$program  < ligand.pdb > tmpfile";
    ( $write_cmd_log ) &&  print CMD_LOG "$command\n";
    (system $command) 
	&& die "Error:\n$command\n"; 

    `echo ligand.pdb > $file`;
    @ret = split " ", `wc -l tmpfile`;
    `echo $ret[0] >> $file`;
    `cat tmpfile >> $file`;
    `echo '  0.000   0.000   0.000 ' >> $file`;
    `rm tmpfile`;



    # make restraint file - doesn't know about coords, so we can do it only once
    $file = "posre_ligand.itp";
    if ( ! -e $file ) {
	$program = $genrestr;
	$command = "echo 0 | $program -f $name.gro -o $file  -fc 1000 1000 1000 > $log  2>&1";
	($write_cmd_log ) &&  print CMD_LOG "$command\n";
	(system $command) 
	    && die "Error:\n$command\n"; 
    }
    `cp $file $posres_in_dir`;

    `mv $name.gro $posres_in_dir`;

    return;
}


sub postprocess_protein () {

    open ( PDB, ">protein.pdb" ) ||
	die "Cno protein.pdb: $!.\n";
    print PDB $protein_pdb;
    close PDB;
    `mv protein.pdb  $posres_in_dir`;
  
   
    $file = "posre.itp";

    if ( ! -e $file ) {

	(-e "../04_production/posre.itp") || 
	    die "../04_production/posre.itp not found\n";
	`ln -s ../04_production/posre.itp .`;

    }
   `cp $file $posres_in_dir`;

    return;
}


############################################################################## 
##############################################################################



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

    print FAKE "; Ligand position restraints\n";
    print FAKE "#ifdef POSRES_LIGAND\n";
    print FAKE "#include \"posre_ligand.itp\"\n";
    print FAKE "#endif\n";

    print FAKE "; Include water topology\n";
    print FAKE "#include \"$forcefield.ff/$water.itp\"\n";

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



############################################################################## 
##############################################################################
sub fix_impl_params (@){


    ($forcefield  eq "amber99sb") 
	|| die "GB params assumes  amber99sb forcefield ...\n";
 
    my $itpfile = $_[0];
    my %params = ();

    # sp2, all-atom, aromatic
    $params{"C"}  = "         0.172    1      1.554    0.1875    0.72 ; C";
    $params{"C*"} = "         0.172    0.012  1.554    0.1875    0.72 ; C";
    $params{"CA"} = "            0.18     1      1.037    0.1875    0.72 ; C";
    $params{"CB"} = "            0.172    0.012  1.554    0.1875    0.72 ; C";
    $params{"CC"} = "           0.172    1      1.554    0.1875    0.72 ; C";

    $params{"CN"} = "           0.172    0.012  1.554    0.1875    0.72 ; C";
    $params{"CR"} = "           0.18     1      1.073    0.1875    0.72 ; C";
    $params{"CV"} = "           0.18     1      1.073    0.1875    0.72 ; C";
    $params{"CW"} = "           0.18     1      1.073    0.1875    0.72 ; C";

    # sp3 all-atom 
    $params{"CT"} = "           0.18     1      1.276    0.190     0.72 ; C";

 
    $params{"H"} =  "           0.1      1      1        0.115     0.85 ; H";
    $params{"HC"} = "           0.1      1      1        0.125     0.85 ; H";
    $params{"H1"} = "           0.1      1      1        0.125     0.85 ; H";
    $params{"HA"} = "           0.1      1      1        0.125     0.85 ; H";
    $params{"H4"} = "           0.1      1      1        0.115     0.85 ; H";
    $params{"H5"} = "           0.1      1      1        0.125     0.85 ; H";
    $params{"HO"} = "           0.1      1      1        0.105     0.85 ; H";
    $params{"HS"} = "           0.1      1      1        0.125     0.85 ; H";
    $params{"HP"} = "           0.1      1      1        0.125     0.85 ; H";

    $params{"N"} =  "           0.155    1      1.028    0.17063   0.79 ; N";
    $params{"NA"} = "           0.155    1      1.028    0.17063   0.79 ; N";
    $params{"NB"} = "           0.155    1      1.215    0.17063   0.79 ; N";
    $params{"N2"} = "           0.16     1      1.215    0.17063   0.79 ; N";
    $params{"N3"} = "           0.16     1      1.215    0.1625    0.79 ; N";

    $params{"O"} = "            0.15     1      0.926    0.148     0.85 ; O";
    $params{"OH"} = "           0.152    1      1.080    0.1535    0.85 ; O";
    $params{"O2"} = "           0.17     1      0.922    0.148     0.85 ; O";
    $params{"S"} = "            0.18     1      1.121    0.1775    0.96 ; S";
    $params{"SH"} = "           0.18     1      1.121    0.1775    0.96 ; S";
    $params{"BR"}= "            0.1      1      1        0.125     0.85 ; BR";
    $params{"F"}= "             0.1      1      1        0.156     0.85 ; F";
    $params{"CL"}= "            0.1      1      1        0.70      0.85 ; CL";
    $params{"I"}= "             0.1      1      1        0.206     0.85 ; I";
    $params{"P5"}= "            0.1      1      1        0.190     0.85 ; P5";
         
    # take as identical (Ivana);
    # these are "united" how did they end up in the same list as all-atom?
    # sp3
    $params{"C2"} =  $params{"CT"}; 
    $params{"C3"} =  $params{"CT"}; 
    $params{"CH"} =  $params{"CT"}; 
    $params{"CS"} =  $params{"CT"}; 
    $params{"CD"} =  $params{"CA"}; # should united atome appear here at all ...?
    $params{"CP"} =  $params{"CA"};
    $params{"CX"} =  $params{"CT"}; # CX - tip of that funny three-membered ring
    $params{"CE"} =  $params{"C*"};
    $params{"CG"} =  $params{"C*"};
    $params{"H2"} =  $params{"H"};

    $params{"NT"} =  $params{"N3"}; # sp3 nitrogen with 4 substituents
    $params{"NH"} =  $params{"N2"}; # sp2 nitrogen in base NH2 group or arginine NH2
    $params{"N*"} =  $params{"N2"}; # sp2 nitrogen in base NH2 group or arginine NH2
    $params{"ND"} =  $params{"N"};  # sp2 nitrogen in amide
    $params{"HN"} =  $params{"H"};  # amide or imino hydrogen


    $params{"OS"} =  $params{"OH"}; # sther or esther O params  replaced by alcohol
    $params{"SS"} =  $params{"S"};    
    $params{"N1"} =  $params{"N2"}; #triple bond in CN?
    $params{"NC"} =  $params{"NB"};
    $params{"NO"} =  $params{"NB"}; # nitrobenzyl?
    $params{"SY"} =  $params{"S"}; #sulfur in 5 ring member
    $params{"SS"} =  $params{"S"}; # sulfur dioxide
    $params{"C1"} =  $params{"C"}; #sp carbon

    `cp $itpfile $itpfile.orig`;
    my @lines = split "\n", `cat $itpfile.orig`;
    my $new_itp = "";
    my $new_field  =   "[ implicit_genborn_params ]\n".
	";name    sar      st     pi       gbr       hct\n";
    my $reading = 0;
    my $name;
   
    foreach my $line (@lines) {

	if ( $line =~ /atomtype/ ) {
	    $reading = 1;

	} elsif ( $reading && $line =~ /\[/ ) {
	    
	    # add the hacked gb field;
	    $new_itp .= $new_field."\n";
	    $reading = 0;

	} elsif  ($reading  && $line =~ /\S/ && $line !~ /name/) {  # not empty or a header a header line
		($name) = split " ", $line;
		$new_field .= " $name ";
		if ( defined $params{$name} ) {
		    $new_field .= $params{$name};

		} elsif ( defined $params{uc $name} ) {
		    $new_field .= $params{uc $name};

		#} elsif (  defined $params{uc substr $name, 0, 1}  ) {
		#    $new_field .= defined $params{uc substr $name, 0, 1} ;

		} else {
		    die "GB params for atom type $name not found.\n";
		}
		$new_field .= "\n";

		
	}

	$new_itp .= $line."\n";

	

    }

    open (NEW_ITP, ">$itpfile") || die "Cno $itpfile: $!\n";
    print NEW_ITP $new_itp;
    close NEW_ITP;

 
    return;
}


=pod

Amber Forcefield Atom List

Atom	Description
C	sp2 carbonyl all-atom carbon and aromatic carbon with hydroxyl substituient in tyrosine
C*	sp2 aromatic all-atom carbon in 5-membered ring with 1 substituent
C2	sp3 united carbon with 2 hydrogens
C3	sp3 united carbon with 3 hydrogens
CA	sp2 aromatic all-atom carbon in 6-membered ring with 1 substituent
CB	sp2 aromatic all-atom carbon in junction between 5- and 6-membered rings
CC	sp2 aromatic all-atom carbon in 5-membered ring with 1 substituent and next to a nitrogen
CD	sp2 aromatic united carbon in 6-membered ring with 1 hydrogen
CE	sp2 aromatic united carbon in 5-membered ring between 2 nitrogens and with 1 hydrogen
CF	sp2 aromatic united carbon in 5-membered ring next to a nitrogen without a hydrogen
CG	sp2 aromatic united carbon in 5-membered ring next to a N-H
CH	sp3 united carbon with 1 hydrogen
CI	sp2 united carbon in 6-membered ring between 2 NC nitrogens
CJ	sp2 united carbon in pyrimidine at positions 5 or 6
CK	sp2 aromatic carbon in 5-membered ring between 2 nitrogens and with 1 hydrogen
CM	sp2 all-atom carbon in pyrimidine at positions 5 or 6
CN	sp2 aromatic junction all-atom carbon in between 5- and 6-membered rings
CP	sp2 aromatic united carbon in 5-membered ring between 2 nitrogens and with 1 hydrogen (in HIS)
CQ	sp2 all-atom carbon in 6-membered ring of purine between two NC nitrogens and with 1 hydrogen
CR	sp2 aromatic all-atom carbon in 5-membered ring between 2 nitrogens and with 1 hydrogen (in HIS)
CT	sp3 all-atom carbon with 4 explicit substituents
CV	sp2 aromatic all-atom carbon in 5-membered ring bonded to 1 N and 1 H
CW	sp2 aromatic all-atom carbon in 5-membered ring bonded to N-H and 1 H
H	amide or imino hydrogen
H2	amino hydrogen in NH2
H3	hydrogen of lysine or arginine (positively charged)
HC	explicit hydrogen attached to carbon
HO	hydrogen in hydroxyl group
HS	hydrogen attached to sulfur
HW	hydrogen in water
LP	lone pair
N	sp2 nitrogen in amide
N*	sp2 nitrogen in purine or pyrimidine wiht alkyl group
N2	sp2 nitrogen in base NH2 group or arginine NH2
N3	sp3 nitrogen with 4 substituents
NA	sp2 nitrogen in 5-membered ring with hydrogen attached
NB	sp2 nitrogen in 5-membered ring with lone pairs
NC	sp2 nitrogen in 6-membered ring with lone pairs
NP	
NT	sp3 nitrogen with 3 substituents
O	carbonyl oxygen
O2	carboxyl or phosphate non-bonded oxygen
OH	alcohol oxygen
OS	ether or ester oxygen
OW	water oxygen
P	phosphorus in phosphate group
S	sulfur in disulfide linkage or methionine
SH	sulfur in thiol
C0	calcium ion (+2)
IM	chlorine ion (-1)
CU	copper ion (+2)
I	Iodine ion (-1)
MG	Magnesium ion (+2)
QC	cesium ion (+1)
QK	potassium ion (+1)
QL	lithium ion (+1)
QN	sodium ion (+1)
QR	rubidium ion (+1)
CS	sp3 carbon
AC	alpha-anomeric carbon
BC	beta-anomeric carbon
HT	sp3 hydrogen
AH	alpha-anomeric hydrogen
BH	beta-anomeric hydrogen
HY	hydroxyl hydrogen
OT	hydroxyl oxygen
OA	alpha-anomeric oxygen
OB	beta-anomeric oxygen
OE	ring oxygen
h$	atom for automatic parameter assignment
c$	atom for automatic parameter assignment
n$	atom for automatic parameter assignment
o$	atom for automatic parameter assignment
s$	atom for automatic parameter assignment
p$	atom for automatic parameter assignment
ospc	in SPC water molecule - used for rattle routine
otip	in TIP3P water molecule - used for rattle routine
=cut
