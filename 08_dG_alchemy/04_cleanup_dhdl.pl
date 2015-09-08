#! /usr/bin/perl -w

=pod
This source code is part of smallmol pipeline for estimate of dG 
upon small modification of a protein molecule
Written by Ivana Mihalek. opyright (C) 2011-2015 Ivana Mihalek.

Gromacs @CCopyright 2015, GROMACS development team. 
Acpype  @CCopyright 2015 SOUSA DA SILVA, A. W. & VRANKEN, W. F.
Gamess  @Copyright 2015m ISUQCG and and contributors to the GAMESS package

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version. This program is distributed in 
the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program. If not, see<http://www.gnu.org/licenses/>.
Contact: ivana.mihalek@gmail.com.
=cut



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
