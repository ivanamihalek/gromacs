#! /usr/bin/perl -w

foreach $type ( "charge", "vdw") {
 
    foreach $i ( 0,1,2,4,6,8,9,10 ) {

	if ( $i < 10 ) {
	    $dir = "0$i\_lambda";
	    $of = "$type.0$i.xvg";
	} else {
	    $dir = "$i\_lambda";
	    $of = "$type.$i.xvg";
	}
	$filename = "$type/$dir/06_production/dhdl.xvg";

	print "$filename\n";
	print "$of\n";
	(-e $of) && next;

	open (IF, "<$filename" ) 
	    || die "Cno $filename: $!.\n";

	open (OF, ">$of" ) 
	    || die "Cno $of: $!.\n";

	while ( <IF> ) {
	    if ( /^\@/ ) {
		( /s[12] / ) ||  print OF $_;
		next;
	    }
	    if ( /^#/ ) {
		print OF $_;
		next;
	    }
	    if ( !/\S/ ) {
		print OF $_;
		next;
	    }
    
	    @aux = split;
    
	    printf OF " %10.4f   %10.4f  \n",  @aux[0 .. 1];

	}
	close IF;
	close OF;
    }

}
