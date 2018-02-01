#! /usr/bin/perl -w

sub round (@);

(@ARGV >= 1 ) || 
    die "Usage: current_raptor.pl <dir/name file>.\n"; 

$infile = shift @ARGV;
$home = `pwd`; chomp $home;
$trjconv = "/usr/local/gromacs/bin/trjconv";
$g_rms = "/usr/local/gromacs/bin/g_rms";

open (IF, "<$infile") ||
    die "Cno $infile:$!.\n";
($stet, $time, $lambda) = ();

$curr_dir = "$infile.trj";

( -e $curr_dir) && `rm -r $curr_dir`;
`mkdir $curr_dir`;
`touch $curr_dir/current_runs`;

$ctr = 0;
$time = "0";
while (<IF>) { 
    ($dir, $name) = split; 
    chdir "$home/$dir";
    print "\n$dir: $name\n";
    print "\t", `ls $name.md.trr`;
    @timelines = split "\n",  `grep -A1 Time md0.log | tail -n10 `;
    
    for ( $t_ctr = $#timelines; $t_ctr >  0; $t_ctr--) {
	next if ($timelines[$t_ctr] !~ /\S/ );
	$timeline2 = $timelines[$t_ctr];
	$timeline1 = $timelines[$t_ctr-1];
	if ( $timeline1 =~ /Time/  && (@aux = split " ", $timeline2) == 3 ) {
	    $time = $aux[1];
	    $time = round ( $time/1000);
	    $time .= "ns";
	    print "\ttime in ns: $time\n";
	    last;
	}
    }
    $ctr ++;


    `echo 2 > twotwo`;
    `echo 2 >> twotwo`;

    $cmd = "echo 0 | $trjconv -s $name.md.tpr -f $name.md.trr -pbc nojump  -o $name.$time.xtc  >& /dev/null";
    system ($cmd) && die "Failure running $cmd.\n";
    $cmd = "$trjconv -s $name.md.tpr -f $name.$time.xtc -fit progressive -o $name.$time.pdb < twotwo >& /dev/null";
    system ($cmd) && die "Failure running $cmd.\n";

    print `ls  $name.$time.pdb`;
    $mangle = "$name.$time.pdb_$ctr";
    `mv  $name.$time.pdb $home/$curr_dir/$mangle`;
    `echo $dir  $name.$time.pdb $mangle >> $home/$curr_dir/current_runs`;
    
=pod
    $cmd = "$g_rms -s $name.md.tpr -f $name.md.trr -o $name.$time.rms.xvg < twotwo >& /dev/null";
    system ($cmd) && die "Failure running $cmd.\n";
    print `ls  $name.$time.rms.xvg`;
    $mangle = "$name.$time.rms_xvg_$ctr";
    `mv  $name.$time.rms.xvg $home/$curr_dir/$mangle`;
    `echo $dir  $name.$time.rms.xvg $mangle >> $home/$curr_dir/current_runs`;
=cut
   

}


close IF;

chdir $home;

`tar -cvf  $curr_dir.tar  $curr_dir`;
( -e "$curr_dir.tar.gz") && `rm  $curr_dir.tar.gz`;
`gzip $curr_dir.tar `;



sub round (@){
    my($number) = shift @_;
    return int($number + .5);
}
