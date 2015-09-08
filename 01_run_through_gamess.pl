#! /usr/bin/perl -w

@ARGV ||
    die "Usage:  $0  <full path to pdb files>  <root name>  <gamess infile header>\n";


($pdb_path, $root, $hdr)  = @ARGV;
$babel  = "/usr/local/bin/babel";
$rungms = "/home/ivanam/chem/gamess/rungms";
$extract_geometry = "/home/ivanam/perlscr/translation/gamess_stdout_2_pdb.pl";

foreach ($pdb_path, $hdr, $babel, $rungms, $extract_geometry) {
    (-e $_) || die "$_ not found\n"; 
}

$home = `pwd`;
chomp $home;

chdir  $pdb_path;
@pdbfiles = split "\n", `ls $root\*.pdb`;
chdir  $home;


@pdbfiles || die "no $root\*.pdb files found in $pdb_path\n";


foreach $pdbfile (@pdbfiles) {

    
    chdir  $home;

    $pdbid = $pdbfile;
    $pdbid =~ s/\.pdb$//;
    (-e $pdbid) || `mkdir $pdbid`;

    chdir $pdbid;
    
    $gamin = "$pdbid.gamin";
    if ( ! -e $gamin  || -z $gamin) {
	$cmd = "$babel  $pdb_path/$pdbfile  $gamin";
	(system $cmd) && die "error running $cmd\n";
	`sed \'s/COORD=CART/COORD=UNIQUE/g\' $pdbid.gamin -i`
    }
    
    $inp = "$pdbid.inp";
    if ( ! -e $inp || -z $inp) {
	`cat ../$hdr $gamin > $inp`;
    }

    $gmslog = "$pdbid.gms_log";
    if ( ! -e $gmslog || -z $gmslog) {
	$cmd = "$rungms  $inp > $gmslog ";
	(system $cmd); # && die "error running $cmd\n";; some crap on exit
    }

    $new_pdb = "$pdbid\_charges.pdb";
    if ( ! -e $new_pdb || -z $new_pdb) {
	$cmd = "$extract_geometry  $pdbid.gms_log $pdbid";
	(system $cmd) && die "error running $cmd\n";
    }
    
}
