#! /usr/bin/perl -w

use Switch;

sub output_atomtypes();
sub output_moleculetype;
sub output_atoms ();

# take that geometries are the same, apart from the new piece 
# or mutated atom type (otherwise the whole thing
# makes no sense, does it

# problems" box vectors are crap - this
# will work if this goes through gro_concat, which
# will discard them and take box vectors from the protein

@ARGV ||
    die "Usage:  <gro from>  <itp from>  <gro to>  <itp to>  <output (root)> \n";

($gro_from, $itp_from, $gro_to,  $itp_to, $output) = ();
($gro_from, $itp_from, $gro_to,  $itp_to, $output) = @ARGV;

# gro line format "%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"




#######################################
# "from" coords
$filename = $gro_from;
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";

<IF>; # skip header
$line_ctr = 1;
$no_atoms  = 100000;
@atom_names = ();

while ( <IF> ) {
    $line_ctr ++;
    if ($line_ctr > $no_atoms+2){
	#$box_line = $_;
	last;
    }

    if ( $line_ctr == 2) {
	chomp;
	$no_atoms = $_;
	$no_atoms =~ s/\s//g;
	
    } else {

	$atom_name = substr $_, 10, 5;
	$atom_name =~ s/\s//g;
	push @atom_names, $atom_name;

	$x{$atom_name} = substr $_, 20, 8;
	$x{$atom_name} =~ s/\s//g;

	$y{$atom_name} = substr $_, 28, 8;
	$y{$atom_name} =~ s/\s//g;

	$z{$atom_name} = substr $_, 36, 8;
	$z{$atom_name} =~ s/\s//g;

	$atom_nr{$atom_name}  =  substr $_, 15, 5;
	$atom_nr{$atom_name}  =~ s/\s//g;
    }
}

close IF;


#######################################
# to coords
$filename = $gro_to;
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";

<IF>; # skip header
$line_ctr = 1;
$no_atoms  = 100000;
@new_atom_names = ();

while ( <IF> ) {
    $line_ctr ++;
    last if ($line_ctr > $no_atoms+2);

    if ( $line_ctr == 2) {
	chomp;
	$no_atoms = $_;
	$no_atoms =~ s/\s//g;
	
    } else {
	$atom_name = substr $_, 10, 5;
	$atom_name =~ s/\s//g;
	push @new_atom_names, $atom_name;

	$new_x{$atom_name} = substr $_, 20, 8;
	$new_x{$atom_name} =~ s/\s//g;

	$new_y{$atom_name} = substr $_, 28, 8;
	$new_y{$atom_name} =~ s/\s//g;

	$new_z{$atom_name} = substr $_, 36, 8;
	$new_z{$atom_name} =~ s/\s//g;

	$new_atom_nr{$atom_name}  =  substr $_, 15, 5;
	$new_atom_nr{$atom_name}  =~ s/\s//g;

    }
}

close IF;



#######################################
# sanity checking

@internal_rep = ();
$ctr = 0;
foreach $atom_name ( @atom_names) {
    # find the atom in the same position in the "to" file
    $maps_to = "";
    foreach $new_atom_name ( @new_atom_names) {
	$d = 0;
	$aux = $x{$atom_name} - $new_x{$new_atom_name};
	$d += $aux*$aux;
	$aux = $y{$atom_name} - $new_y{$new_atom_name};
	$d += $aux*$aux;
	$aux = $z{$atom_name} - $new_z{$new_atom_name};
	$d += $aux*$aux;
	$d = sqrt ($d);
	if ( $d < 0.01) {
	    $maps_to = $new_atom_name;
	    last;
	}

    }

    # outcomes:
    if ( $maps_to ) {

	#$old_map{$atom_name} = $maps_to;
	$new_map{$maps_to}   = $atom_name;

	push @internal_rep, "$atom_name  $atom_nr{$atom_name}   $maps_to  $new_atom_nr{$maps_to}";
	$name2internal{$atom_name} = $ctr;
	$number2internal{$atom_nr{$atom_name}} = $ctr;

	$newname2internal{$maps_to} = $ctr;
	$newnumber2internal{$new_atom_nr{$maps_to}} = $ctr;



	$ctr++;
    }  else {

	push @internal_rep, "$atom_name  $atom_nr{$atom_name} -  - ";
	$name2internal{$atom_name} = $ctr;
	$number2internal{$atom_nr{$atom_name}} = $ctr;

	$ctr++;

    }
}


foreach $new_atom_name ( @new_atom_names) {
    (defined $new_map{$new_atom_name}) && next;
    
    push @internal_rep, " -  -   $new_atom_name  $new_atom_nr{$new_atom_name}";

    $newname2internal{$new_atom_name} = $ctr;
    $newnumber2internal{$new_atom_nr{$new_atom_name}} = $ctr;

    $ctr++;

}


################################
# assign new names to the internal rep
($oldname, $oldnum, $newname, $newnum_i) = ();
foreach $intern ( 0 .. $#internal_rep) {
    ($oldname, $oldnum, $newname, $newnum_i) = split " ", $internal_rep[$intern];
    
    if ( $oldname eq "-" ) {
	$newname =~ /(\D+)(\d*)/;
	$new_type = $1;
	$type = $new_type;
    } elsif ( $newname eq "-" ) {
	$oldname  =~ /(\D+)(\d*)/;
	$old_type = $1;
	$type     = $old_type;
    } else {
	$newname =~ /(\D+)(\d*)/;
	$new_type = $1;
	$oldname  =~ /(\D+)(\d*)/;
	$old_type = $1;
	
	if ( $new_type eq $old_type ) {
	    $type =  $new_type;
	} else  {
	    $type = $old_type.$new_type;
	}

    }
    (defined $count{$type}) || ($count{$type} = 0);
	
    $count{$type}++;
    $intern_name[$intern] = $type.$count{$type};

}

################################
foreach $intern ( 0 .. $#internal_rep) {

    ($oldname, $oldnum, $newname, $newnum_i) = split " ", $internal_rep[$intern];
    
    printf " %4s   %4s  %4s   %4s  %4s  \n", 
    $intern_name[$intern], $oldname, $oldnum, $newname, $newnum_i;
}
exit;

################################


open ( VDW_ITP, ">$output.vdw.itp") ||
    die "Cno $output.itp: $!\n";

open ( CHG_ITP, ">$output.chg.itp") ||
    die "Cno $output.itp: $!\n";


################################
################################
################################
################################
# parse both input itp files

foreach $input_itp ( $itp_from, $itp_to ) {
 
    $filename = $input_itp;
    open (ITP, "<$filename" ) 
	|| die "Cno $filename: $!.\n";

    @{$sections{$input_itp}} = ();
    $section_name = "";
    $dihedrals    = 0;

    while ( <ITP> )  {


	if ( /\[(.+?)\]/ ) {
	    $section_name = $1; 
	    $section_name =~ s/\s//g;
	    if (  $section_name eq "dihedrals" ) {
		$dihedrals++;
		if ( $dihedrals == 2 ) {
		    $section_name = "impropers";
		}
	    }
	    
	} elsif ( /^\s*\;/ && $section_name) {
	    $comment{$section_name} = $_;

	} elsif  ( !/\S/ )  {
	    $section_name = "";

	} elsif ($section_name)  {

	    if ( ! defined $section{$input_itp}{$section_name})  {
		push @{$sections{$input_itp}}, $section_name;
		$section{$input_itp}{$section_name} = "";
	    }
	    $section{$input_itp}{$section_name} .= $_;
	}
	
    }
    close ITP;
}




################################
################################
################################
#  write ITP files

open ( VDW_ITP, ">$output.vdw.itp") ||
    die "Cno $output.itp: $!\n";

open ( CHG_ITP, ">$output.chg.itp") ||
    die "Cno $output.itp: $!\n";


$date = `date`; chomp $date;

print VDW_ITP "; $output.vdw.itp created by $0 on $date\n";
print CHG_ITP "; $output.chg.itp created by $0 on $date\n";


foreach $section_name ( "atomtypes", "moleculetype", "atoms" ) {

    (defined $section{$itp_from}{$section_name}) ||
	(defined $section{$itp_to}{$section_name}) || next;

    print VDW_ITP "\n[ $section_name ]\n";
    print CHG_ITP "\n[ $section_name ]\n";
    if  ( defined  $comment{$section_name} ) {
	print VDW_ITP  $comment{$section_name};
	print CHG_ITP  $comment{$section_name};
    }


    switch ($section_name) {
	case "atomtypes"     {output_atomtypes()}
	case "moleculetype"  {output_moleculetype()}
	case "atoms"         {output_atoms()}
    }
}

print VDW_ITP  "\n";
print CHG_ITP  "\n";

close VDW_ITP;
close CHG_ITP;


######################################################################
######################################################################
######################################################################
######################################################################

sub output_atomtypes() {
    @atom_types = ();

    foreach $input_itp ( $itp_from, $itp_to ) {
	( defined $section{$input_itp}{$section_name} ) || next;
	foreach ( split "\n", $section{$input_itp}{$section_name} ) {
	    chomp;
	    @aux = split;
	    $name_bond = "$aux[0] $aux[1]";
	    if ( ! defined $atomtypes{$name_bond} ) {
		$atomtypes{$name_bond} = $_;
		push @atom_types, $name_bond;
	    }
	}
    }
 
    foreach $name_bond ( @atom_types ) {
	print VDW_ITP  "$atomtypes{$name_bond}\n";
	print CHG_ITP  "$atomtypes{$name_bond}\n";
    }
    # dummy atom type
    print VDW_ITP  "dum      dum          0.00000  0.00000   A     0.00000e+00   0.00000e+00 ; 0.00  0.0000\n";
}

################################
sub output_moleculetype() {

    foreach $input_itp ( $itp_from, $itp_to ) {
	( defined $section{$input_itp}{$section_name} ) || next;
	foreach ( split "\n", $section{$input_itp}{$section_name} ) {
	    chomp;
	    @aux = split;
	    $nrexcl{$input_itp} = pop @aux;
	}
    }
    if ($nrexcl{$itp_from} != $nrexcl{$itp_to})  {
	die "nrexcl in $itp_from different from the one in $itp_to. (?!)\n";
    }
    print VDW_ITP  " $output  $nrexcl{$itp_to} \n";
    print CHG_ITP  " $output  $nrexcl{$itp_to} \n";
}

################################

sub output_atoms () {
    ($oldname, $oldnum, $newname, $newnum)  = ();

    if ( defined $section{$itp_to}{$section_name} ) {
	foreach ( split "\n", $section{$itp_to}{$section_name} ) {
	    
	    $atom_name = substr $_, 25, 5;
	    $atom_name =~ s/\s//g;

	    $new_atom_line{$atom_name} = $_;
	    
	}
    }

    if ( defined $section{$itp_from}{$section_name} ) {
	foreach $old_atom_line ( split "\n", $section{$itp_from}{$section_name} ) {

	    $atom_number = substr $old_atom_line, 0, 6;
	    $atom_number =~ s/\s//g;
	
	    $intern = $number2internal{$atom_number};
	    
	    $atom_name = substr $old_atom_line, 25, 5;
	    $atom_name =~ s/\s//g;

	    $res_name = substr $old_atom_line, 20, 3;
	    $res_name =~ s/\s//g;
	    

	    ($oldname, $oldnum, $newname, $newnum) = split " ", $internal_rep[$intern];

	    $atom_type_old = $atom_name;
	    $atom_type_old =~ s/\d//g;
    
	    $type_old  = substr $old_atom_line, 7, 5;
	    $type_old  =~ s/\s//g;

	    $mass_old  = substr $old_atom_line, 51, 10;
	    $mass_old  =~ s/\s//g;

	    $charge  = substr $old_atom_line, 35, 13;
	    $charge  =~ s/\s//g;

	    if ( $newname eq "-" ) {
		# disappears

		printf VDW_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f  \n",
		$intern+1, $type_old, 1, $res_name, $intern_name[$intern], $intern+1, 
		$charge, $mass_old, "dum", 0.0, $mass_old;
		next;
	    }
	    $maps_to = $newname;

    
	    $atom_type_new = $maps_to;
	    $atom_type_new =~ s/\d//g;
 
	    $type_new  = substr $new_atom_line{$maps_to}, 7, 5;
	    $type_new  =~ s/\s//g;

	    $mass_new  = substr $new_atom_line{$maps_to}, 51, 10;
	    $mass_new  =~ s/\s//g;

	    $charge_new  = substr $new_atom_line{$maps_to}, 35, 13;
	    $charge_new  =~ s/\s//g;


	    if ( $atom_type_new  eq $atom_type_old) {

		printf VDW_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f  \n",
		$intern+1, $type_old, 1, $res_name, $intern_name[$intern], $intern+1, $charge, $mass_old,
		$type_old, $charge, $mass_old;

	    } else {

		printf VDW_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f  \n",
		$intern+1, $type_old, 1, $res_name, $intern_name[$intern], $intern+1, $charge, $mass_old,
		$type_new, 0.0, $mass_new;
		
	    }
	}
    }

    (defined $section{$itp_to}{$section_name}) || return;

    # do we have any new atoms here?
    foreach $new_line (split "\n", $section{$itp_to}{$section_name} ) {

	$type_new  = substr $new_line, 7, 5;
	$type_new  =~ s/\s//g;

	$mass_new  = substr $new_line, 51, 10;
	$mass_new  =~ s/\s//g;


	$charge_new = substr $new_line, 35, 13;
	$charge_new =~ s/\s//g;

	    
	$new_atom_name = substr $new_line, 25, 5;
	$new_atom_name =~ s/\s//g;
	$intern = $newname2internal{$new_atom_name};

	($oldname, $oldnum, $newname, $newnum) = split " ", $internal_rep[$intern];
	

	printf CHG_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f  \n",
	$newnum, $type_new, 1, $res_name, $intern_name[$intern], $newnum, 0.0, $mass_new,
	$type_new, $charge_new, $mass_new;


	($oldname eq "-" ) || next;

 
	printf VDW_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f  \n",
	$intern+1, "dum", 1, $res_name, $intern_name[$intern], $intern+1, 0.0, $mass_new,
	$type_new, 0.0, $mass_new;

    
    }

}


######################################################################


