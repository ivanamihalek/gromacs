#! /usr/bin/perl -w

# 49  LJ-14:LIG-LIG                       50  Coul-SR:LIG-SOL                   
# 51  LJ-SR:LIG-SOL                       52  LJ-LR:LIG-SOL                     
# 53  Coul-14:LIG-SOL                     54  LJ-14:LIG-SOL                     

# for mdrun -rerun
# 49  Coul-SR:MOL-SOL                     50  LJ-SR:MOL-SOL                     
# 51  Coul-LR:MOL-SOL                     52  Coul-14:MOL-SOL                   
# 53  LJ-14:MOL-SOL                       54  Coul-SR:SOL-SOL   
        
@cmpds =  ( "orig_cmpd", "cmpd11", "cmpd10", "cmpd9", "cmpd8", 
	    "cmpd7", "cmpd6", "cmpd5", "cmpd4", "cmpd3", "cmpd2", "cmpd1");

@runs = ("00_full_system",  "01_lig+solv",  
	 "02_lig+solv_uncharged");

foreach $cmpd ( @cmpds) {
    foreach $run (@runs) {
	( -e "$run/$cmpd") || die "$run/$cmpd not found\n";
    }
}


$home = `pwd`; chomp $home;
($name, $avg, $err, $crap, $crap) = ();

foreach $cmpd ( @cmpds) {

    foreach $run (@runs[1..2]) {
	chdir $home;
	chdir  "$run/$cmpd";


	$edr = "$cmpd.md.edr";
	$cmd = "(echo 50 && echo 53 && echo 51 && echo 52  && echo 54) "
	    ."| g_energy -f $edr | grep \'kJ/mol\'";
	$ret = `$cmd`;

	$ret_orig{$run} = $ret;


	# go for the rerun thing

=pod
  #should I really do that - the correction
  #from  Could LR is in the 4th place,
  #I'm losing info about LJ-lR  in taht case
  #and the samppling is only on the handful
  #of trajectory snapshots ...
  
	$edr =  "rerun.edr";

	if ( ! -e $edr) {
	    # create rerun mdp
	    open (IF, "<md.mdp") || 
		die "Cno md.mdp: $!.\n";
	    open (OF, ">rerun.mdp") || 
		die "Cno rerun.mdp: $!.\n";
	
	    foreach $line ( <IF> ) {
		if ( $line =~ /coulombtype/) {
		    print OF "\n";
		} elsif ( $line =~ /coulombtype        =   Cut-off/) {
		    print OF "\n";
	    
		} elsif ( $line =~ /rcoulomb/) {
		    print OF "rcoulomb            =  1.7\n";

		} elsif ( $line =~ /rlist/) {
		    print OF "rlist               =  1.6\n";

		} elsif ( $line =~ /rvdw/) {
		    print OF "rvdw                =  1.6\n";
		
		} elsif ( $line =~ /nstlist/) {
		    print OF "nstlist             =  1\n";

		} else {
		    print OF $line;
		}

	    }
	    close OF;
	    close IF;
	
	    # grompp
	    $cmd = "grompp -f rerun.mdp -c $cmpd.water.gro -p $cmpd.top -o rerun.tpr";
	    (system $cmd) && die "Error running $cmd.\n";

	    # mdrun -rerun
	    $cmd = "mdrun -rerun  $cmpd.md.trr -s rerun.tpr -e $edr ";
	    (system $cmd) && die "Error running $cmd.\n";
	}
=cut
	# g_eneregy
	$cmd = "(echo 49 && echo 50 && echo 51 && echo 52  && echo 53 ) "
	    ."| g_energy -f $edr | grep \'kJ/mol\'";
	$ret = `$cmd`;

	$ret{$run} = $ret;


	@lines = split "\n", $ret;
	foreach $line (@lines) {
	    ($name, $avg, $err, $crap, $crap) = split " ", $line;
	    ($name, $lig) = split "\:", $name;
	    $term{$run}{$name} = $avg;
	}
    }

    foreach $run (@runs[1..2]) {
	print "$run:\n";
	foreach $name ( keys %{$term{$run}} ) {
	    print "$name  $term{$run}{$name}\n";
	}
    }
    $Eqq{$cmpd} = $term{"01_lig+solv"}{"Coul-SR"}+
	$term{"01_lig+solv"}{"Coul-14"};
 
    (defined $term{"01_lig+solv"}{"Coul-LR"}) &&  ($Eqq{$cmpd}+=$term{"01_lig+solv"}{"Coul-LR"});

    $Elj{$cmpd} = $term{"01_lig+solv"}{"LJ-SR"}+$term{"01_lig+solv"}{"LJ-LR"}+
	$term{"01_lig+solv"}{"LJ-14"};

    $correction = $term{"02_lig+solv_uncharged"}{"LJ-SR"}+
	$term{"02_lig+solv_uncharged"}{"LJ-LR"}+
	$term{"02_lig+solv_uncharged"}{"LJ-14"};


    $Eqq{$cmpd} += $Elj{$cmpd} - $correction;



    $run = "00_full_system";
    chdir $home;
    chdir "$run/$cmpd";
    $edr = "2xdlA.md.edr";

    
    $lig =~ s/\-SOL//;
    $cmd = "g_lie -f $edr -Eqq $Eqq{$cmpd} -Elj $Elj{$cmpd} ".
	"-ligand $lig | grep DGbind";
    $ret = `$cmd`;
    ($name, $crap, $dgb, $crap) = split " ", $ret;
    $dgbind{$cmpd} = $dgb;
    #print"$cmd\n";
    #last if ($cmpd eq "cmpd10");
}


foreach $cmpd (@cmpds) {
    printf "%10s   %8.3f    %8.3f    %8.3f  \n", $cmpd,
    $Eqq{$cmpd}, $Elj{$cmpd}, $dgbind{$cmpd};
    
}


=pod
foreach $run (@runs[1..2]) {
    print "$run:\n";
    print "$ret_orig{$run}\n";
    print "----------------------\n";
    print "$ret{$run}\n\n";
}
=cut
