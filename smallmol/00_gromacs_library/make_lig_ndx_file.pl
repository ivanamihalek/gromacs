#! /usr/bin/perl -w

(@ARGV >= 2) ||
    die "Usage:  $0  <gro file>  <lig name>\n";

($filename, $lig) = @ARGV;
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";

open (OF, ">lig.ndx") 
    || die "Cno lig.ndx.\n";

print "energygrps  ";
while ( <IF> ) {
    /^Protein/ && next;
    $resn = substr $_, 5, 3;
    ($resn ne $lig) && next;


    $atomtype = substr $_, 11, 4; $atomtype =~ s/\s//g;
    $num      = substr $_, 15, 5; $num =~ s/\s//g;

    

    print OF "[ L$atomtype ]\n  $num\n" ;

    print "  L$atomtype ";

}
print "\n";
close IF;
close OF

#123456789012345678901234567890123456789012345678901234567890123456789
# 224LIG     H8 3301   2.988   3.519   2.907 -1.0517  0.6115 -0.0309
 
