#! /usr/bin/perl -w
(@ARGV ==3) ||
    die "Usage:  $0  <full path to acpype dir>  <root name (single)>  <receptor pdb>\n";

($acpype_path, $root, $rec)  = @ARGV;


$home = `pwd`;
chomp $home;

$acpype_dir = "$acpype_path/$root.acpype";
$generic_mdps_dir = "/home/ivanam/perlscr/gromacs/generic_input/leapfrog";

foreach ($acpype_path, $acpype_dir, $generic_mdps_dir, $rec) {
    (-e $_) || die "$_ not found\n"; 
}


(-e $root) || `mkdir $root`;
(-e "$root/start") || `mkdir $root/start`;
chdir "$root/start";

$grofile = "$root\_GMX.gro";
(-e "$acpype_dir/$grofile") || die "$grofile not found in $acpype_dir\n";
 n(-e "$root.gro") || `cp $acpype_dir/$grofile $root.gro`;

$ret = `head -n3 $root.gro | tail -n 1`;

@aux = split " ", $ret;
$ligname = $aux[1];

$itpfile = "$root\_GMX.itp";
(-e "$acpype_dir/$itpfile") || die "$itpfile not found in $acpype_dir\n";
(-e "$root.itp") || `cp $acpype_dir/$itpfile $root.itp`;

foreach ("em.mdp", "pr.mdp", "md.mdp") {
    (-e "$generic_mdps_dir/$_") || die "$_ not found in $generic_mdps_dir\n";
    `grep -v energygrps $generic_mdps_dir/$_ > $_`;
    ($_ =~ "md.mdp") || next;
    `echo energygrps          = $ligname  SOL  protein >> $_`;
     
}

@aux = split "\/", $rec;

$rec_pdb = pop @aux;
print "copying  $rec to ".`pwd`;
(-e "$rec_pdb") || `cp $rec $rec_pdb`;

`echo $root 1 > ligands`;

`cp * ../`;

print " cd $root && /home/ivanam/perlscr/gromacs/sloppy_gmx.pl  2xdlA  ligands &\n";
