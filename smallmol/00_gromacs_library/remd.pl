#! /usr/bin/perl -w

use strict;
sub  error_state_check (@ ); # defined at the bottom

defined $ARGV[0] ||
    die "Usage: $0 <pdb name root> [<list of ligand names>].\n";

$| = 1; # turn on autoflush;

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
my $gromacs_perl = "/home/ivanam/perlscr/gromacs/gromacs.pl";


my $hostname = `hostname`;

chomp $hostname;


if ( $np ) {
    $mdrun = "/home/imihalek/gromacs/gromacs-3.3.1/src/kernel/mdrun";
    $gromacs_run = "$mpirun  -machinefile \$PBS_NODEFILE -np $np $mdrun -np $np";
    $grompp_root = "/home/imihalek/gromacs/gromacs-3.3.1/src/kernel/grompp";
    $grompp  = "$grompp_root  -np $np";
    $groc    = "/home/imihalek/perlscr/gromacs/gro_concat.pl";
    $pdb2gro = "/home/imihalek/perlscr/pdb2gro.pl";
} else {
   
    $gromacs_run = "$gromacs_path/mdrun";
    $grompp      = "$gromacs_path/grompp";
    $groc        = "/home/ivanam/perlscr/gromacs/gro_concat.pl";
    $pdb2gro     = "/home/ivanam/perlscr/translation/pdb2gro.pl";
}

if ($hostname eq "donkey") {
    $mdrun       = "/usr/bin/mdrun";
    $mdrun_mpi   = "/usr/bin/mdrun_mpi";
    $mpirun      = "/usr/bin/mpirun";
} else {
    $mdrun_mpi = "";
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
my $genion_solvent_code = 12;

if (  $forcefield  eq  "oplsaa" || $forcefield  eq "amber99sb" ) {
    $neg_ion  = "CL";
    $pos_ion  = "NA";
}

##############################################################################
##############################################################################
## in  ~/gromacs/generic_mdps
foreach ( $gromacs_perl, "start", "start/$name.pdb", "start/em.mdp", "start/pr.mdp", 
	  "start/md.mdp", "start/temp_distr" ) {
    ( -e $_ ) || die "\n$_ not found.\n\n";
}
##############################################################################
##############################################################################
my $home = `pwd`;
chomp $home;
my $cmd;

##############################################################################
##############################################################################
# read in temp distr file and create equilibratiion directories
my ($replica_num);
my $temp;
my @replica_temp;

open (TDF, "<start/temp_distr")
    || die "CNo tempstart/_distr: $!\n";

@replica_temp  = ();

while (<TDF>) {
    next if ( !/\S/);
    next if ( /^#/);
    ($replica_num, $temp) = split " ", $_;
    push  @replica_temp, $temp;
}
close TDF;


##################################################
#
$np = 8;
($hostname eq "donkey") && ($np = 16);

##################################################
#   energy minimization - the same for all
my $mindir = "minimization";
(-e $mindir) || `mkdir $mindir`;
chdir $mindir;
print "\n**************************************\n";
print "minimization\n\n";
if ( ! -e "$name.em_output.trr" ) {	
`cp ../start/* .`;
$cmd = "$gromacs_perl @ARGV -min";
(system $cmd) && 
    die "Error running $cmd in\n";
} else {
 print "$name.em_output.trr found\n";

}
chdir $home;



##################################################
# create one directory for each replica
my $eqdir = "replica_equilibration";
( -e $eqdir) || `mkdir $eqdir`;
my $repdir;



foreach  $replica_num (0 .. $np-1) {

     chdir $home;
     $repdir = "$eqdir/replica_$replica_num";
     ( -e $repdir) || `mkdir $repdir`;
     chdir $repdir;
     
     (-e "$ARGV[0].md.tpr") && next;


     `ln -s $home/$mindir/* .`;

     $temp = $replica_temp[$replica_num];

     foreach my $mdpfile ( "pr.mdp", "md.mdp") {
	 my @lines = split "\n", `cat $mdpfile`;
	 open ( OF, ">$home/$mindir/$mdpfile" ) || die "Cno $home/$mindir/$mdpfile: $!\n";
	 foreach my $line (@lines) {
	     if ( $line =~ "ref_t") {
		 print OF "ref_t               =   $temp    $temp \n";
	     } elsif  ( $line =~ "gen_temp") {
		 print OF "gen_temp            =   $temp    \n";
	     } else {
		 print OF "$line\n";
	     }
	 }
	 close OF;
     }

     print "\n**************************************\n";
     print " running equilibration for replica $replica_num ($temp K)\n";
     $cmd = "$gromacs_perl @ARGV -remd";
     (system $cmd) && 
	 die "Error running $cmd in $repdir\n";
     

}



#############################################################
# collect all replicas and run the actual REMD sim
chdir $home;
my $prod = "production_run";
(-e $prod) || `mkdir $prod`;


foreach  $replica_num (0 .. $np-1)          {

    $repdir = "$eqdir/replica_$replica_num";
    chdir $home;
    chdir $repdir;
    
    `mv $ARGV[0].md.tpr  $home/$prod/$ARGV[0]\_$replica_num.tpr`;

}


chdir $home;
chdir $prod;

$cmd = "$mpirun -np $np $mdrun_mpi -np $np -s $ARGV[0]\_.tpr -multi $np -replex 1000 -o $ARGV[0].trr";
(system $cmd) && 
    die "Error running $cmd.\n";


exit;
