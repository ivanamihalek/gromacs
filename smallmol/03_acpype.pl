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

# the pdbfiles here are supposed to be already optimized in gamess

@ARGV ||
    die "Usage:  $0  <full path to gamess dir>  <root name> \n";


($pdb_path, $root)  = @ARGV;
$acpype = "../00_gromacs_library/acpype.py";



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

