#!/usr/bin/perl -w

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




sub q_run (@);
$max_jobs = 7;

###########################################
###########################################
###########################################
#
#  check for all dependencises 
#
###########################################

# the number of steps for each type of run
# no protein
$no_steps[0]{"em_steep"}  = 10000;
$no_steps[0]{"em_lbfgs"}  =  5000;
$no_steps[0]{"pr_nvt"}    = int (200.0/2e-3); # 200ps
$no_steps[0]{"pr_npt"}    = int (200.0/2e-3); # 200ps 
$no_steps[0]{"md"}        = int ( 10.0/2e-6); # 10ns

# with protein
$no_steps[1]{"em_steep"}  = 10000;
$no_steps[1]{"em_lbfgs"}  = 5000;
$no_steps[1]{"pr_nvt"}    = int (100.0/2e-3); # 100ps
$no_steps[1]{"pr_npt"}    = int (100.0/2e-3); # 100ps
$no_steps[1]{"md"}        = int (  3.0/2e-6); # 3ns;


# lambdas are hardcoded for now:
($l_start, $l_step_size, $no_of_steps) = (0.0, 0.1, 11);



$gmx = "/home/ivanam/perlscr/gromacs/gromacs.pl";
( -e "$gmx" ) || die "Main ingredient missing: $gmx.\n";

( -e "gmx_input" ) || die "Main ingredient missing: gmx_input directory.\n";

foreach ( "itps", "gros", "mdps") {



sub q_run (@);
$max_jobs = 7;

###########################################
###########################################
###########################################
#
#  check for all dependencises 
#
###########################################

# the number of steps for each type of run
# no protein
$no_steps[0]{"em_steep"}  = 10000;
$no_steps[0]{"em_lbfgs"}  =  5000;
$no_steps[0]{"pr_nvt"}    = int (200.0/2e-3); # 200ps
$no_steps[0]{"pr_npt"}    = int (200.0/2e-3); # 200ps 
$no_steps[0]{"md"}        = int ( 10.0/2e-6); # 10ns

# with protein
$no_steps[1]{"em_steep"}  = 10000;
$no_steps[1]{"em_lbfgs"}  = 5000;
$no_steps[1]{"pr_nvt"}    = int (100.0/2e-3); # 100ps
$no_steps[1]{"pr_npt"}    = int (100.0/2e-3); # 100ps
$no_steps[1]{"md"}        = int (  3.0/2e-6); # 3ns;


# lambdas are hardcoded for now:
($l_start, $l_step_size, $no_of_steps) = (0.0, 0.1, 11);



$gmx = "/home/ivanam/perlscr/gromacs/gromacs.pl";
( -e "$gmx" ) || die "Main ingredient missing: $gmx.\n";

( -e "gmx_input" ) || die "Main ingredient missing: gmx_input directory.\n";

foreach ( "itps", "gros", "mdps") {
   ( -e "gmx_input/$_" ) || die "Main ingredient missing: $_ subdirectory.\n";
}

# do we have the protein pdb file?
@gmx_input_contents = split "\n", `ls gmx_input`;

$pdb_id = "";
foreach $name ( @gmx_input_contents ) {
    ( -d $name) && next; # this is directory
    ($name =~ /pdb$/) || next;

    $pdb_id = $name;
    $pdb_id =~ s/\.pdb//;
}
$pdb_id || die "no pdb file found in gmx_input\n";
print "\n***********************\n $pdb_id.pdb  found in gmx_input - I'm taking it is the protein file.\n";


############################################
$home = `pwd`; chomp $home;
($nr, $name_root) = ();

foreach $system_dir ("no_protein", "with_protein" ) {
    
    foreach $subdir ( "04_from", "05_to") {
	($nr, $name_root) = split "_", $subdir;
	foreach $input_type ( "itp", "gro" ) {
	    $input_file = "$home/gmx_input/$input_type"."s/$name_root.$input_type";
	    ( -e $input_file) || die "$input_file not found\n";
	}
    }

    foreach $subdir ( "01_chg_off", "02_vdw", "03_chg_on") {
	@aux  = split "_", $subdir;
	shift @aux;
	$name_root  = join "_", @aux;
	foreach $input_type ( "itp", "gro" ) {
	    $input_file = "$home/gmx_input/$input_type"."s/mutant.$name_root.$input_type";
	    ( -e $input_file) || die "$input_file not found\n";
	}
    }

}



###########################################
###########################################
###########################################
#
#  run - make directory structure as we go
#
###########################################

@running_jobs = ();
($nr, $name_root) = ();

foreach $system_dir ("no_protein", "with_protein" ) {
    
    chdir $home;
    ( -e $system_dir) || `mkdir  $system_dir`;
    chdir $system_dir;

    

    foreach $subdir ( "04_from", "05_to") {

	chdir $home;
	chdir $system_dir;
	( -e $subdir) || `mkdir $subdir`;
	chdir $subdir;

	($nr, $name_root) = split "_", $subdir;

	print "\n***********************\n $system_dir    $subdir   \n";

	$inp_dir = "00_input";
	( -e $inp_dir) || `mkdir $inp_dir`;

	foreach $input_type ( "itp", "gro" ) {
	    $input_file = "$home/gmx_input/$input_type"."s/$name_root.$input_type";
	    `cp  $input_file  $inp_dir`;
	}	

	`cp $home/gmx_input/mdps/* $inp_dir`;
	($system_dir eq "with_protein") && `cp $home/gmx_input/$pdb_id.pdb $inp_dir`;
	`echo $name_root 1 > $inp_dir/ligands`;
	# get rid of free energy, but add energygrps to be used by LIE later
	chdir $inp_dir;
	# find out the name of the ligand in the gro file
	# for that I  need
	# gro line format "%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
	$ret = `head -n3 $name_root.gro | tail -n1`;
	$ligand_name = substr $ret, 5, 5;
	$ligand_name =~ s/\s//g;
	foreach $mdp (split "\n", `ls *mdp`) {
	   

	    $runtype = $mdp;
	    $runtype =~ s/\.mdp//;
	    @lines = split "\n", `cat $mdp`;
	    # we'll rewrite the original mdp
	    open (OF, ">$mdp") ||
		die "Cno $mdp: $!";

	    foreach $line (@lines) {
		($line =~ /free_energy/)  && next;
		($line =~ /lambda/)  && next;
		if ( $line =~ /nsteps/ ) {
		    # how many steps do we want?
		    $with_protein = ($system_dir eq "with_protein") ? 1 : 0;
		    $nst = $no_steps[$with_protein]{$runtype};
		    print OF "nsteps  = $nst \n";
		    next;
		}
		print OF "$line\n";
	    }

	    if ( $system_dir eq "with_protein") {
		print OF "energygrps = $ligand_name SOL protein \n";
	    } else {
		print OF "energygrps = $ligand_name SOL  \n";	
	    }
	    close OF;

	}
	chdir "..";
	
	
	if ( $system_dir eq "with_protein") {
	    $cmd = "$gmx  $pdb_id   ligands > /dev/null";
	} else {
	    $cmd = "$gmx  no_protein  ligands > /dev/null";
	}
	
	(q_run $cmd) && die "error running $cmd ";
    }


    foreach $subdir ( "01_chg_off", "02_vdw", "03_chg_on") {

	chdir $home;
	chdir $system_dir;
	( -e $subdir) || `mkdir $subdir`;
	chdir $subdir;

	@aux  = split "_", $subdir;
	shift @aux;
	$name_root  = join "_", @aux;



	# here we need to further chop the thing into lambdas
	foreach $i ( 0 .. $no_of_steps-1 ) {

	    chdir "$home/$system_dir/$subdir";

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
	    (-e $dir) || `mkdir  $dir`;
	    chdir $dir;

	    print "\n***********************\n $system_dir    $subdir   $dir\n";

	    $inp_dir = "00_input";
	    ( -e $inp_dir) || `mkdir $inp_dir`;
	    foreach $input_type ( "itp", "gro" ) {
		$input_file = "$home/gmx_input/$input_type"."s/mutant.$name_root.$input_type";
		`cp  $input_file  $inp_dir/mutant.$input_type`;
	    }


	    `cp $home/gmx_input/mdps/* $inp_dir`;
	    ($system_dir eq "with_protein") && `cp $home/gmx_input/$pdb_id.pdb $inp_dir`;
	    `echo mutant 1 > $inp_dir/ligands`;
	
	    chdir $inp_dir;
	   
	    @mdps = split "\n", `ls *.mdp`;
	    foreach $mdp ( @mdps ) {

		$runtype = $mdp;
		$runtype =~ s/\.mdp//;

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
		
		    } elsif ($line =~ /nsteps/) {
			# how many steps do we want?
			$with_protein = ($system_dir eq "with_protein") ? 1 : 0;
			$nst = $no_steps[$with_protein]{$runtype};
			print OF "nsteps  = $nst \n";
		    } else {
			print OF "$line\n";
		    }
		}
		close OF;
	    }
	    `rm tmp`;
	    chdir "..";


	
	    if ( $system_dir eq "with_protein") {
		$cmd = "$gmx  $pdb_id   ligands > /dev/null";
	    } else {
		$cmd = "$gmx  no_protein  ligands > /dev/null";
	    }
	
	    (q_run $cmd) &&  die "Error running $cmd ";
	}
    }

}


##########################################
##########################################

sub q_run (@) {

    my $cmd = shift @_;
    my $pid;

    while ( @running_jobs >= $max_jobs ) {
	# how many jobs
	print "\nrunning jobs greater than max jobs ($max_jobs)\n";
	print "@running_jobs\n";
	# check if all jobs alive:
	sleep 600; # check every 10mins
	@prev_jobs    =  @running_jobs;
	@running_jobs = ();
	foreach $pid (@prev_jobs) {
	    my $alive = "";
	    $alive = `ps -u ivanam | grep -v defunct |  awk \'\$1==$pid\'`;
	    if ( $alive ) {
		push @running_jobs, $pid;
	    }
	}
    }

    #$cmd = "sleep 10";
    $pid = fork();

    if ($pid) { # the parent process

	push @running_jobs, $pid;

    } else {

	#print "$pid: $cmd\n";
	exec $cmd;
	exit;
	#print "$pid: done\n";
    }
    return 0;
}
