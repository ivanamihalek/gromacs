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


die "Set up the babel and GAMESS paths below, then comment out this line.\n";

@ARGV ||
    die "Usage:  $0  <full path to pdb files>  <root name>  <gamess infile header>\n";


($pdb_path, $root, $hdr)  = @ARGV;
$babel  = "/usr/local/bin/babel";
$rungms = "/home/ivanam/chem/gamess/rungms";
$extract_geometry = "../00_gromacs_library/gamess_stdout_2_pdb.pl";

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
