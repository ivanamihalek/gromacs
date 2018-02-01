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

(@ARGV ==2) ||
    die "Usage:  $0  <full path to acpype dir>  <root name (single)> \n";

($acpype_path, $root)  = @ARGV;


$home = `pwd`;
chomp $home;

$acpype_dir = "$acpype_path/$root.acpype";
$generic_mdps_dir = "/home/ivanam/perlscr/gromacs/generic_input/leapfrog";

foreach ($acpype_path, $acpype_dir, $generic_mdps_dir) {
    (-e $_) || die "$_ not found\n"; 
}


(-e $root) || `mkdir $root`;
(-e "$root/start") || `mkdir $root/start`;
chdir "$root/start";

$grofile = "$root\_GMX.gro";
(-e "$acpype_dir/$grofile") || die "$grofile not found in $acpype_dir\n";
(-e "$root.gro") || `cp $acpype_dir/$grofile $root.gro`;

$itpfile = "$root\_GMX.itp";
(-e "$acpype_dir/$itpfile") || die "$itpfile not found in $acpype_dir\n";
$fixing = 0;
if ( ! -e "$root.itp") {
    open (IF, "< $acpype_dir/$itpfile") || 
	die "Cno $:  $acpype_dir/$itpfile$!.\n";
    open (OF, ">$root.itp") || 
	die "Cno: $root.itp $!.\n";

    foreach $line ( <IF> ) {
	if ( $line =~ /atoms/) {
	    $fixing = 1;
	    print  OF  $line;
	    next;
	} elsif ( $line =~ /bonds/){
	    $fixing = 0;
	}

	if ( $fixing ) {
	    if ( $line =~/^;/ || $line !~ /\S/) {
		print OF  $line;
	    }  else {
		chomp $line;
		@aux = split " ", $line;
		# charg is in the 7th column
		$aux[6] = " 0.000 ";
		print OF  join "   " , @aux;
		print OF "\n";
	    }

	} else {
	    print  OF  $line;
	}
    }
    close OF;
    close IF;
}

foreach ("em.mdp", "pr.mdp", "md.mdp") {

    (-e "$generic_mdps_dir/$_") || die "$_ not found in $generic_mdps_dir\n";

    open (IF, "<$generic_mdps_dir/$_") || 
	die "Cno $generic_mdps_dir/$_: $!.\n";
    open (OF, ">$_") || 
	die "Cno $_: $!.\n";
    
    foreach $line ( <IF> ) {
	if ( $line =~ /tau_t/) {
	    print OF "tau_t                 =     0.1 \n";
	} elsif ( $line =~ /tc_grps/) {
	    print OF "tc_grps               =     system \n";
	    
	} elsif ( $line =~ /ref_t/) {
	    print OF "ref_t                 =     300 \n";
	} else {
	    print OF $line;
	}


	
    }
    ($_ =~ "md.mdp") || next;
    print OF "energygrps          = LIG SOL\n";
    close OF;

    close IF;
}


`echo $root 1 > ligands`;

`cp * ../`;

print " cd $root && ../00_gromacs_library/gromacs/sloppy_gmx.pl  no_protein ligands &\n";
