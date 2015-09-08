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

(@ARGV ==3) ||
    die "Usage:  $0  <full path to acpype dir>  <root name (single)>  <receptor pdb>\n";

($acpype_path, $root, $rec)  = @ARGV;


$home = `pwd`;
chomp $home;

$acpype_dir = "$acpype_path/$root.acpype";
$generic_mdps_dir = "../00_gromacs_library/generic_input/leapfrog";

foreach ($acpype_path, $acpype_dir, $generic_mdps_dir, $rec) {
    (-e $_) || die "$_ not found\n"; 
}


(-e $root) || `mkdir $root`;
(-e "$root/start") || `mkdir $root/start`;
chdir "$root/start";

$grofile = "$root\_GMX.gro";
(-e "$acpype_dir/$grofile") || die "$grofile not found in $acpype_dir\n";
 n(-e "$root.gro") || `cp $acpype_dir/$grofile $root.gro`;

$ret = `head -n3 $root.gro | tail -n 1`;

@aux = split " ", $ret;
$ligname = $aux[1];

$itpfile = "$root\_GMX.itp";
(-e "$acpype_dir/$itpfile") || die "$itpfile not found in $acpype_dir\n";
(-e "$root.itp") || `cp $acpype_dir/$itpfile $root.itp`;

foreach ("em.mdp", "pr.mdp", "md.mdp") {
    (-e "$generic_mdps_dir/$_") || die "$_ not found in $generic_mdps_dir\n";
    `grep -v energygrps $generic_mdps_dir/$_ > $_`;
    ($_ =~ "md.mdp") || next;
    `echo energygrps          = $ligname  SOL  protein >> $_`;
     
}

@aux = split "\/", $rec;

$rec_pdb = pop @aux;
print "copying  $rec to ".`pwd`;
(-e "$rec_pdb") || `cp $rec $rec_pdb`;

`echo $root 1 > ligands`;

`cp * ../`;

print " cd $root && /home/ivanam/perlscr/gromacs/sloppy_gmx.pl  2xdlA  ligands &\n";
