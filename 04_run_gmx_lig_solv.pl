#! /usr/bin/perl -w
(@ARGV ==2) ||
    die "Usage:  $0  <full path to acpype dir>  <root name (single)> \n";

($acpype_path, $root)  = @ARGV;


$home = `pwd`;
chomp $home;

$acpype_dir = "$acpype_path/$root.acpype";
$generic_mdps_dir = "/home/ivanam/perlscr/gromacs/generic_input/leapfrog";

foreach ($acpype_path, $acpype_dir, $generic_mdps_dir) {
    (-e $_) || die "$_ not found\n"; 
}


(-e $root) || `mkdir $root`;
(-e "$root/start") || `mkdir $root/start`;
chdir "$root/start";

$grofile = "$root\_GMX.gro";
(-e "$acpype_dir/$grofile") || die "$grofile not found in $acpype_dir\n";
(-e "$root.gro") || `cp $acpype_dir/$grofile $root.gro`;

$itpfile = "$root\_GMX.itp";
(-e "$acpype_dir/$itpfile") || die "$itpfile not found in $acpype_dir\n";
(-e "$root.itp") || `cp $acpype_dir/$itpfile $root.itp`;

foreach ("em.mdp", "pr.mdp", "md.mdp") {

    (-e "$generic_mdps_dir/$_") || die "$_ not found in $generic_mdps_dir\n";

    open (IF, "<$generic_mdps_dir/$_") || 
	die "Cno $generic_mdps_dir/$_: $!.\n";
    open (OF, ">$_") || 
	die "Cno $_: $!.\n";
    
    foreach $line ( <IF> ) {
	if ( $line =~ /tau_t/) {
	    print OF "tau_t                 =     0.1 \n";
	} elsif ( $line =~ /tc_grps/) {
	    print OF "tc_grps               =     system \n";
	    
	} elsif ( $line =~ /ref_t/) {
	    print OF "ref_t                 =     300 \n";
	} else {
	    print OF $line;
	}


	
    }
    ($_ =~ "md.mdp") || next;
    print OF "energygrps          = LIG SOL\n";
    close OF;

    close IF;
}


`echo $root 1 > ligands`;

`cp * ../`;

print " cd $root && /home/ivanam/perlscr/gromacs/sloppy_gmx.pl  no_protein ligands &\n";
