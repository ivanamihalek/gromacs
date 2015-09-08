#! /usr/bin/perl -w


(@ARGV >= 3) || 
    die "Usage: $0 <l start>  <l step size>  <# of l steps>  [no_protein|<protein name>]\n";

($l_start, $l_step_size, $no_of_steps, $mode) = @ARGV;

if ( defined $mode ) {
    $gmx = "/home/ivanam/perlscr/gromacs/gromacs.pl";
    foreach ( "start", $gmx) {    
	(-e $_) || die "$_ not found\n";
    }
}

$home = `pwd`;
chomp $home;

foreach $i ( 0 .. $no_of_steps-1 ) {

    $lambda = $l_start+$i*$l_step_size;
     
    $label = int (100*$lambda);
    if ( $label < 10 ) {
	$label = "00$label";
    } elsif (  $label < 100 ) {
	$label = "0$label";

    }elsif (  $label > 100 ) {
	die "$lambda > 1 (?)\n";
    }

    $dir = "$label\_lambda";

    print "$dir\n";

    chdir $home;
    ( -e $dir) || `mkdir $dir`;

    chdir $dir;
    ( -e "start") || `cp -r ../start 00_input`;
    
    chdir  "00_input";

    @mdps = split "\n", `ls *.mdp`;
    foreach $mdp ( @mdps ) {
	`cp $mdp tmp`;
	@lines = split "\n", `cat tmp`;

	open (OF, ">$mdp") ||
	    die "Cno $mdp: $!";

	foreach $line (@lines) {
	    if ($line =~ /init_lambda/) {
		printf OF "init_lambda              =%5.2f\n", $lambda;
	    } elsif ($line =~ /foreign_lambda/) {

		printf OF "foreign_lambda           =";
		if ( $lambda-$l_step_size >= 0.0 ) {
		    printf OF " %5.2f ", $lambda-$l_step_size;
		} 
		if ( $lambda+$l_step_size <= 1.0 ) {
		    printf OF " %5.2f ", $lambda+$l_step_size;
		} 
		print  OF "\n";
		
	    } else {
		print OF "$line\n";
	    }
	}
	close OF;

    }
    `rm tmp`;
    if ( defined $mode) {
	chdir "..";
	$cmd = "$gmx  $mode ligands";
	$pid = fork();
	$pid || (exec $cmd);
    }

}
