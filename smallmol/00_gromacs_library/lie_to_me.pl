#! /usr/bin/perl -w



$gromacs_path = "/usr/local/bin";
$make_ndx     = "$gromacs_path/make_ndx";
$tpbconv      = "$gromacs_path/tpbconv";
$mdrun        = "$gromacs_path/mdrun";
$g_energy     = "$gromacs_path/g_energy";
$g_lie        = "$gromacs_path/g_lie";
$perl_path    = "/home/ivanam/perlscr/gromacs";

foreach ($gromacs_path, $perl_path, $make_ndx, $tpbconv, $mdrun, $g_energy, $g_lie ) {

    (-e $_) || die "$_ not found\n";
}



(  -e "ligand" &&  -e "protein") || 
    die  "both protein and ligand dirs ".
    "must be present to do LIE calc.\n";


$begin_water = "";
$end_water   = "";
$begin       = "";
$end         = "";

if (@ARGV > 0 ) {
    
    if ( @ARGV < 4 ) {
	die "If you choose to give arguments to $0, you have to give all four:\n".
	    "$0  <begin time in water>   <end time in water>  ".
	    "<begin time in protein>  <end time in protein>.\n".
	    "(all in picoseconds).\n";
    }
    
    ($begin_water, $end_water, $begin, $end) = @ARGV;
}


$home = `pwd`; chomp $home;


$outfile = "energy_calc.log";
open (OF, ">$outfile") ||
    die "Cno $outfile: $!.\n";


##########################################
#
# process the ligand in water simulation
#

$simdir  = "ligand";
chdir "$home/$simdir";
$trr = "";
$tpr = "";

if ( -e "07_extension" ) {

    chdir "07_extension";

    #########################################################
    $trr = "";
    @trrs = split "\n", `ls *trr`;
    ( @trrs == 0  ) && 
	die "no trr file found in $simdir/07_extension\n";

    if (@trrs > 1) {
	$trr = "";
 	foreach  (@trrs){
	    if (/concat/ ) {
		$trr = $_;
		last;
	    } 
	}
	if ( ! $trr ) {
	    foreach  (@trrs){
		if (/ext/ ) {
		    $trr = $_;
		    last;
		} 
	    }
	}
	$trr || 
	    die "multiple trr files found in $simdir/07_extension,"
	    ." neither one having concat or ext in the name:"
	    ."@trrs\n";

	$trr = "07_extension/$trr";
    } else {

	$trr = "07_extension/$trrs[0]";

    }
    #########################################################
    $tpr = "";
    @tprs = split "\n", `ls *tpr`;
    ( @tprs == 0  ) && 
	die "no tpr file found in $simdir/07_extension\n";

    if (@tprs > 1) {
	$tpr = "";
 	foreach  (@tprs){
	    if (/concat/ ) {
		$tpr = $_;
		last;
	    } 
	}
	if ( ! $tpr ) {
	    foreach  (@tprs){
		if (/ext/ ) {
		    $tpr = $_;
		    last;
		} 
	    }
	}
	$tpr || 
	    die "multiple tpr files found in $simdir/07_extension,"
	    ." neither one having concat or ext in the name:"
	    ."@tprs\n";

	$tpr = "07_extension/$tpr";
    } else {

	$tpr = "07_extension/$tprs[0]";

    }
    #########################################################
    $gro = "";
    @gros = split "\n", `ls *gro`;
    ( @gros == 0  ) && 
	die "no gro file found in $simdir/07_extension\n";

    if (@gros > 1) {
	$gro = "";
 	foreach  (@gros){
	    if (/concat/ ) {
		$gro = $_;
		last;
	    } 
	}
	if ( ! $gro ) {
	    foreach  (@gros){
		if (/ext/ ) {
		    $gro = $_;
		    last;
		} 
	    }
	}
	$gro || 
	    die "multiple gro files found in $simdir/07_extension,"
	    ." neither one having concat or ext in the name:"
	    ."@gros\n";

	$gro = "07_extension/$gro";
    } else {

	$gro = "07_extension/$gros[0]";

    }
    #########################################################
    $edr = "";
    @edrs = split "\n", `ls *edr`;
    ( @edrs == 0  ) && 
	die "no edr file found in $simdir/07_extension\n";

    if (@edrs > 1) {
	$edr = "";
 	foreach  (@edrs){
	    if (/concat/ ) {
		$edr = $_;
		last;
	    } 
	}
	if ( ! $edr ) {
	    foreach  (@edrs){
		if (/ext/ ) {
		    $edr = $_;
		    last;
		} 
	    }
	}
	$edr || 
	    die "multiple edr files found in $simdir/07_extension,"
	    ." neither one having concat or ext in the name:"
	    ."@edrs\n";

	$edr = "07_extension/$edr";
    } else {

	$edr = "07_extension/$edrs[0]";

    }
    #########################################################
 

} elsif ( -e "06_production" ){

    chdir "06_production";

    @trrs = split "\n", `ls *trr | grep -v zeroq`;
    (@trrs == 0) && 
	die "no trr file found in $simdir/06_production\n";

    (@trrs > 1) && 
	die "multiple trr files found in $simdir/06_production\n:"
	."@trrs\n";;
    $trr = "06_production/$trrs[0]";


    @tprs = split "\n", `ls *tpr | grep -v zeroq`;
    (@tprs == 0) && 
	die "no tpr file found in $simdir/06_production\n";

    (@tprs > 1) && 
	die "multiple tpr files found in $simdir/06_production\n:"
	."@tprs\n";;
    $tpr = "06_production/$tprs[0]";


    @gros = split "\n", `ls *gro | grep -v zeroq`;
    (@gros == 0) && 
	die "no gro file found in $simdir/06_production\n";

    (@gros > 1) && 
	die "multiple gro files found in $simdir/06_production\n:"
	."@gros\n";;
    $gro = "06_production/$gros[0]";


    @edrs = split "\n", `ls *edr | grep -v zeroq`;
    (@edrs == 0) && 
	die "no edr file found in $simdir/06_production\n";

    (@edrs > 1) && 
	die "multiple edr files found in $simdir/06_production\n:"
	."@edrs\n";;
    $edr = "06_production/$edrs[0]";


} else {
    
    die "neither 06_production nor 07_extension dirs found for $simdir.\n";

}

printf "$simdir:    $trr   $tpr   $gro  $edr \n"; 
chdir "$home/$simdir";

$outdir = "10_analysis";
(-e $outdir) || `mkdir $outdir`;

chdir $outdir;

if ( ! -e "index.ndx" ) {
    $cmd = "(echo 2; echo keep 2;  echo q) | $make_ndx -f  ../$gro";
    (system $cmd) && die "Error running $cmd.\n";
}


if ( ! -e "zeroq.tpb" ) {
    $cmd = "$tpbconv -s ../$tpr -f ../$trr -n index.ndx -zeroq -o zeroq.tpb";
    (system $cmd) && die "Error running $cmd.\n";
}

if ( ! -e "zeroq.edr" ) {
    $cmd = "$mdrun -s zeroq.tpb -rerun ../$trr -deffnm zeroq";
    (system $cmd) && die "Error running $cmd.\n";
}


$time_span = "";
if ( $begin_water &&  $end_water) {
    $time_span = "  -b $begin_water  -e $end_water ";
}

$cmd = "(echo 55) | $g_energy -f zeroq.edr $time_span | tail -n 50 | tee tmp";
(system $cmd) && die "Error running $cmd.\n";




($field_name, $avg, $err, $rmsd, $drift) = ();
($field_name, $avg, $err, $rmsd, $drift) = split " ", `grep \'LJ-SR\' tmp`;
$clean_lj_in_water = $avg;
$clean_lj_in_water_err = $err;

$time_span = "";
if ( $begin_water &&  $end_water) {
    $time_span = "  -b $begin_water  -e $end_water ";
}

$cmd = "(echo 54; echo 55) | $g_energy -f ../$edr $time_span | tail -n 50 | tee tmp";
(system $cmd) && die "Error running $cmd.\n";

($field_name, $avg, $err, $rmsd, $drift) = ();
($field_name, $avg, $err, $rmsd, $drift) = split " ", `grep \'LJ-SR\' tmp`;
$compounded_lj_in_water = $avg;
$compounded_lj_in_water_err = $err;

($field_name, $avg, $err, $rmsd, $drift) = ();
($field_name, $avg, $err, $rmsd, $drift) = split " ", `grep \'Coul-SR\' tmp`;
$compounded_coul_in_water = $avg;
$compounded_coul_in_water_err = $err;

$sum     = $compounded_lj_in_water     + $compounded_coul_in_water;
$sum_err = $compounded_lj_in_water_err + $compounded_coul_in_water_err;

$clean_coul_in_water =  sprintf "  %8.3f ", $sum-$clean_lj_in_water;


printf OF " clean_lj_in_water   =  %8.3f   pm  %8.3f  kJ/mol \n", 
    $clean_lj_in_water, $clean_lj_in_water_err;
printf OF "              sum    =  %8.3f   pm  %8.3f  kJ/mol\n", 
    $sum, $sum_err;
printf OF " clean_coul_in_water =  %8.3f   pm  %8.3f  kJ/mol \n", 
    $sum - $clean_lj_in_water, $sum_err+$clean_lj_in_water_err;

`rm -f *.[1-9]*  tmp `;





##########################################
#
# now move to protein
#

$simdir  = "protein";
chdir "$home/$simdir";
$trr = "";
$tpr = "";

$edr = "";


if ( -e "07_extension" ) {

    chdir "07_extension";

   
    @edrs = split "\n", `ls *edr`;

    (@edrs == 0) && 
	die "no edr file found in $simdir/07_extension\n";

 
    if (@edrs > 1) {
	$edr = "";
	foreach  (@edrs){
	    if (/concat/ ) {
		$edr = $_;
		last;
	    } 
	}
	if ( ! $edr ) {
	    foreach  (@edrs){
		if (/ext/ ) {
		    $edr = $_;
		    last;
		} 
	    }
	}

	$edr || 
	    die "multiple edr files found in $simdir/07_extension,"
	    ." neither one having concat or ext in the name:"
	    ."@edrs\n";
    }
    $edr = "07_extension/$edr";
  

} elsif ( -e "06_production" ){

    chdir "06_production";


    @edrs = split "\n", `ls *edr`;
    (@edrs == 0) && 
	die "no edr file found in $simdir/06_production\n";

    (@edrs > 1) && 
	die "multiple edr files found in $simdir/06_production\n:"
	."@edrs\n";;
    $edr = "06_production/$edrs[0]";


} else {
    
    die "neither 06_production nor 07_extension dirs found for $simdir.\n";

}

printf "$simdir:   $edr \n"; 

chdir "$home/$simdir";

#figure out the ligand name
$ligand = "";

@itps = split "\n", `ls 01_topology/*.itp | grep -v posre.itp`;
(@itps == 0 ) && die "no itp found in $home/$simdir/01_topology";
(@itps >  1 ) && die "more than one itp found in $home/$simdir/01_topology";



$ret = `grep moleculetype $itps[0] -A 2 | tail -n 1`;


@aux = split " ", $ret;

$ligand = $aux[0];
$ligand = uc $ligand;



$outdir = "10_analysis";
(-e $outdir) || `mkdir $outdir`;
chdir $outdir;


$time_span = "";
if ( $begin &&  $end) {
    $time_span = "  -b $begin  -e $end ";
}
 
$cmd = "$g_lie -f ../$edr -Elj $clean_lj_in_water  ".
    " -Eqq $clean_coul_in_water  -ligand $ligand  $time_span | tail -n 50 | tee tmp";

(system $cmd) && die "Error running $cmd.\n";

$ret = `grep DGbind tmp`;
chomp $ret;
$ret =~ s/[\=\(\)]//g;
($field_name, $dgbind, $err) = split " ", $ret;

printf OF "\n";
printf OF "         DG binding =  %8.3f   pm  %8.3f  kJ/mol\n",    $dgbind, $err;
printf OF "         DG binding =  %8.3f   pm  %8.3f  kcal/mol\n",  $dgbind/4.184, $err/4.184;
    
printf OF "\n";

$expfactor =  $dgbind/4.184*1.678;
$K_avg     =  exp($expfactor)*1.e6;
$expfactor =  ($dgbind+$err)/4.184*1.678;
$lower_K   =  exp($expfactor)*1.e6;
$expfactor =  ($dgbind-$err)/4.184*1.678;
$upper_K   =  exp($expfactor)*1.e6;

printf OF "         Ka at 300K =   %8.3e uM\n", $K_avg ;  

printf OF "         Ka range   =  [ %8.3e .. %8.3e ]  uM\n", $lower_K , $upper_K;  

`rm -f *.[1-9]*  tmp `;


close OF;
