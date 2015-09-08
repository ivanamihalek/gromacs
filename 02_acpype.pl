#! /usr/bin/perl -w
# the pdbfiles here are supposed to be already optimized in gamess

@ARGV ||
    die "Usage:  $0  <full path to gamess dir>  <root name> \n";


($pdb_path, $root)  = @ARGV;
$acpype = "/home/ivanam/docking/02_ligand_prep/acpype/acpype.py";



foreach ($pdb_path, $acpype) {
    (-e $_) || die "$_ not found\n"; 
}

$home = `pwd`;
chomp $home;

chdir  $pdb_path;
@pdbids = split "\n", `ls -d $root\*`;
chdir  $home;

@pdbids || die "no $root\*.pdb files found in $pdb_path\n";

#print join "\n", @pdbids;
#print "\n";
#exit;


foreach $pdbid (@pdbids) {

    printf "\n $pdbid:\n";
    $gmslog = "$pdb_path/$pdbid/$pdbid.gms_log";
    if ( ! -e $gmslog) {
	print "$gmslog not found\n";
	next;
    }

    $cmd = "tail $gmslog | ".
	" grep \'EXECUTION OF GAMESS TERMINATED NORMALLY\'";
    $ret = "" || `$cmd`;
    if ( ! $ret) {
	print "apparently, gamess failed for $pdbid\n";
	next;
    }
    print "$ret\n";


    $pdbfile = "$pdb_path/$pdbid/$pdbid\_charges.pdb";
    if ( ! -e $pdbfile) {
	printf "$pdbfile not found\n";
	next;
    }

    chdir  $home;
    (-e  "$pdbid.pdb") || `ln -s $pdbfile $pdbid.pdb`;

    if ( ! -e "$pdbid.acpype" ) {
	$cmd = "$acpype -i $pdbid.pdb";
	(system $cmd) && die "error running $cmd\n";
    }


    #exit;

   
}

