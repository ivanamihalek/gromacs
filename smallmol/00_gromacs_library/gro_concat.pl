#! /usr/bin/perl -w
use IO::Handle;         #autoflush
# FH -> autoflush(1);

(defined $ARGV[1]) ||
    die "Usage: gro_concat.pl <gro1> <gro2> ...\n";

@grofiles = @ARGV;

$output = "";
$res_ctr = 0;
$atom_ctr = 0; 
$box_vectors = "";
foreach $grofile ( @grofiles ) {

    $resno_old   = -1;
    $resname_old = "";

    open (GF, "<$grofile") ||
	die "Error: opening $grofile: $!.\n";
    while ( <GF> ) {
	next if ( ! /\S/);
	# check  whether this is the coordinate line (it should end with a float number)
	@aux = split;
	$test = pop @aux;
	next if ( $test !~ /\d\.\d/);

	$resno   = substr $_, 0, 5;
	$resname = substr $_, 5, 5;
	#print "$resname\n";
	if ( @aux == 2 ) { # I did pop above
	    $box_vectors || ($box_vectors = $_); #assume these are ok
	                       # basically, I am assuming that I am concating one big
                               # and several smaller moelcules
	    next;
	}
	$atom_ctr ++;
	if ( $resno ne $resno_old || $resname ne $resname_old) {
	    $res_ctr ++;
	    $resno_old  = $resno;
	    $resname_old  = $resname;
	}
	$aux = sprintf "%5d",  $res_ctr;
	(substr $_, 0, 5) = $aux;
	$aux = sprintf "%5d",  $atom_ctr;
	(substr $_, 15, 5) = $aux;

	$resname =~ s/\s//g; 
	$aux = sprintf "%-5s", $resname;
	( substr $_, 5, 5) = $aux;

	$output .= $_;
   }
    close GF;
}


print "some name\n";
printf "%5d\n", $atom_ctr;
print $output;
print $box_vectors;




#12345678901234567890123456789012345678901234567890
#    1CIT  O1       1  -2.124   5.795  10.179
#  754PRO     O2 7400  -0.680   4.604  12.266
#C format
#    "%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" 
