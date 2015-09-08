#! /usr/bin/perl -w

# take that geometries are the same, apart from the new piece 
# or mutated atom type (otherwise the whole thing
# makes no sense, does it

# problems" box vectors are crap - this
# will work if this goes through gro_concat, which
# will discard them and take box vectors from the protein

@ARGV ||
    die "Usage:  <gro from> <itp from>  <gro to>  <itp to>  <output (root)> \n";

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
	$box_line = $_;
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
# "to" coords
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
 	if ( $d < 0.001) {
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
####################################
foreach $intern ( 0 .. $#internal_rep) {

    ($oldname, $oldnum, $newname, $newnum_i) = split " ", $internal_rep[$intern];
    if ( $oldname eq "-" ) {
	$appears[$intern] = 1;
    } else {
	$appears[$intern] = 0;
    }
    if ( $newname eq "-" ) {
	$disappears[$intern] = 1;
    } else {
	$disappears[$intern] = 0;
    }
}
####################################
=pod
foreach $intern ( 0 .. $#internal_rep) {

    ($oldname, $oldnum, $newname, $newnum_i) = split " ", $internal_rep[$intern];
     
    printf " %4s   %4s  %4s   %4s  %4s  ", 
    $intern_name[$intern], $oldname, $oldnum, $newname, $newnum_i;
    if ( 	$disappears[$intern] ) {
	printf " disappears ";
    }
    if ( 	$appears[$intern] ) {
	printf " appears ";
    }
    print "\n";
}
exit;
=cut
####################################



open ( CHG_OFF_ITP, ">$output.chg_off.itp") ||
    die "Cno $output.chg_off.itp: $!\n";

open ( VDW_ITP, ">$output.vdw.itp") ||
    die "Cno $output.itp: $!\n";

open ( CHG_ON_ITP, ">$output.chg_on.itp") ||
    die "Cno $output.chg_on.itp: $!\n";


################################
# open and get to the beginning of the read of both files
$filename = $itp_from;
open (ITP_FROM, "<$filename" ) 
    || die "Cno $filename: $!.\n";

while ( <ITP_FROM> ) {
    print CHG_OFF_ITP ;
    print VDW_ITP ;
    print CHG_ON_ITP ;
    last if /atomtypes/;
}
$filename = $itp_to;
open (ITP_TO, "<$filename" ) 
    || die "Cno $filename: $!.\n";

while ( <ITP_TO> ) {
    last if /atomtypes/;
}



################################
# read in the atom types
@atom_types = ();
while ( <ITP_FROM> ) {
    if ( /^\;/) {
	print CHG_OFF_ITP ;
	print VDW_ITP ;
	print CHG_ON_ITP ;
	next;
    }
    next if ( !/\S/);
    
    if ( /moleculetype/) {
	last;
    }
    chomp;
    @aux = split;
    $name_bond = "$aux[0] $aux[1]";
    if ( ! defined $atomtypes{$name_bond} ) {
	$atomtypes{$name_bond} = $_;
	push @atom_types, $name_bond;
    }
}
################################
# any other atom types in 
while ( <ITP_TO> ) {
    next if ( /^\;/);
    next if ( !/\S/);
    last if /moleculetype/;
    chomp;
    @aux = split;
    $name_bond = "$aux[0] $aux[1]";
    if ( ! defined $atomtypes{$name_bond} ) {
	$atomtypes{$name_bond} = $_;
	push @atom_types, $name_bond;
    }
}

foreach $name_bond ( @atom_types ) {
    print CHG_OFF_ITP  "$atomtypes{$name_bond}\n";
    print VDW_ITP  "$atomtypes{$name_bond}\n";
    print CHG_ON_ITP  "$atomtypes{$name_bond}\n";
}
# dummy atom type
print VDW_ITP  "dum      dum          0.00000  0.00000   A     0.00000e+00   0.00000e+00 ; 0.00  0.0000\n";

######################################################################
######################################################################
######################################################################
# moleculetype
$reading = 0;
$nrexcl  = 3;

while ( <ITP_FROM> ) {
    last if /atoms/;
    if ( /nrexcl/) {
	$reading = 1;
	next;
    }
    next if ( !/\S/);
    if ( $reading ) {
	chomp;
	@aux = split;
	$nrexcl = pop @aux;
	$reading = 0;
    }
}

print CHG_OFF_ITP  "\n[ moleculetype ]\n";
print CHG_OFF_ITP  ";name            nrexcl\n";
print CHG_OFF_ITP  " $output  $nrexcl \n";

print VDW_ITP  "\n[ moleculetype ]\n";
print VDW_ITP  ";name            nrexcl\n";
print VDW_ITP  " $output  $nrexcl \n";

print CHG_ON_ITP  "\n[ moleculetype ]\n";
print CHG_ON_ITP  ";name            nrexcl\n";
print CHG_ON_ITP  " $output  $nrexcl \n";

while ( <ITP_TO> ) {
    last if /atoms/;
}

######################################################################
######################################################################
# atoms section
# is this fixed format, or acpype format?
#1234567890123456789012345678901234567890123456789012345678901234567890
#    10   c3     1   MOL    C8   10     0.091000     12.01000 ; qtot -1.321

@old_atom_lines = ();
print CHG_OFF_ITP  "\n[ atoms ]\n";
print VDW_ITP  "\n[ atoms ]\n";
print CHG_ON_ITP  "\n[ atoms ]\n";
while ( <ITP_FROM> ) {
    if ( /^\;/) {
	print CHG_OFF_ITP ;
	print VDW_ITP ;
	print CHG_ON_ITP ;
	next;
    }
    next if ( !/\S/);
    last if ( /bonds/ );

    push @old_atom_lines, $_;
    
}

while ( <ITP_TO> ) {
    next if ( /^\;/);
    
    next if ( !/\S/);
    last if ( /bonds/ );

    $atom_name = substr $_, 25, 5;
    $atom_name =~ s/\s//g;

    $new_atom_line{$atom_name} = $_;
    
}


($oldname, $oldnum, $newname, $newnum)  = ();

foreach $old_atom_line (@old_atom_lines) {

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

    $charge_old  = substr $old_atom_line, 35, 13;
    $charge_old  =~ s/\s//g;

   

    if ( $disappears[$intern] ) {

	printf VDW_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f  ; disappears \n",
	$intern+1, $type_old, 1, $res_name, $intern_name[$intern], $intern+1, 
	0.0, $mass_old, "dum", 0.0, $mass_old;
	
	printf CHG_OFF_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f  ; disappears\n",
	$oldnum, $type_old, 1, $res_name, $intern_name[$intern], $oldnum, $charge_old, $mass_old,
	$type_old,  0.0, $mass_old;

	
    } else {

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
	    $intern+1, $type_old, 1, $res_name, $intern_name[$intern], $intern+1, $charge_old, $mass_old,
	    $type_old, $charge_old, $mass_old;

	    printf CHG_OFF_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f  \n",
	    $oldnum, $type_old, 1, $res_name, $intern_name[$intern], $oldnum, $charge_old, $mass_old,
	    $type_old,  $charge_old, $mass_old;

 
	} else {

	    printf VDW_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f ; changes type \n",
	    $intern+1, $type_old, 1, $res_name, $intern_name[$intern], $intern+1, 0.0, $mass_old,
	    $type_new, 0.0, $mass_new;
	    
	    printf CHG_OFF_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f ; changes type \n",
	    $oldnum, $type_old, 1, $res_name, $intern_name[$intern], $oldnum, $charge_old, $mass_old,
	    $type_old,  0.0, $mass_old;
	}
    }
}


# do we have any new atoms here?
$atom_number = scalar @old_atom_lines;
foreach $new_atom_name (@new_atom_names) {



    $new_line = $new_atom_line{$new_atom_name};
   

    $type_new  = substr $new_line, 7, 5;
    $type_new  =~ s/\s//g;

    $mass_new  = substr $new_line, 51, 10;
    $mass_new  =~ s/\s//g;


    $charge_new = substr $new_line, 35, 13;
    $charge_new =~ s/\s//g;

    $intern = $newname2internal{$new_atom_name};

    ($oldname, $oldnum, $newname, $newnum) = split " ", $internal_rep[$intern];


    if ( $appears[$intern]  ) {

	printf CHG_ON_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f ; appears  \n",
	$newnum, $type_new, 1, $res_name, $intern_name[$intern], $newnum, 
	0.0, $mass_new, $type_new, $charge_new, $mass_new;

	printf VDW_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f ; appears \n",
	$intern+1, "dum", 1, $res_name, $intern_name[$intern], $intern+1, 0.0, $mass_new,
	$type_new, 0.0, $mass_new;

    } else {

	printf CHG_ON_ITP  "%6d%5s%6d   %3s%6s%5d  %11.6f %11.5f  %5s %11.6f %11.5f  \n",
	$newnum, $type_new, 1, $res_name, $intern_name[$intern], $newnum, 
	$charge_new, $mass_new, $type_new, $charge_new, $mass_new;

	# we wrote it already to VDW file

	
    }


    
}


######################################################################
######################################################################
######################################################################
######################################################################
##############################
#   bonds
@old_bond_lines = ();
print CHG_OFF_ITP  "\n[ bonds ]\n";
print VDW_ITP      "\n[ bonds ]\n";
print CHG_ON_ITP   "\n[ bonds ]\n";
while ( <ITP_FROM> ) {
    if ( /^\;/) {
	print CHG_OFF_ITP ;
	print VDW_ITP ;
	print CHG_ON_ITP ;
	next;
    }
    next if ( !/\S/);
    last if ( /pairs/ );

    push @old_bond_lines, $_;    
}

@new_bond_lines = ();
while ( <ITP_TO> ) {
    next if ( /^\;/);
    
    next if ( !/\S/);
    last if ( /pairs/ );

    $atom_name = substr $_, 25, 5;
    $atom_name =~ s/\s//g;

    ($ai,     $aj, $funct,   $r,  $k) = split;

    # make sure that the numbering between the old 
    # and the new works
    $r_new{$ai}{$aj} = $r_new{$aj}{$ai}  = $r;
    $k_new{$ai}{$aj} = $k_new{$aj}{$ai}  = $k;

    push @new_bond_lines, $_;
}

($ai, $aj, $funct,   $r,  $k)  = ();
($oldname_i, $oldnum_i, $newname_i, $newnum_i) =  ();
($oldname_j, $oldnum_j, $newname_j, $newnum_j) = ();


foreach  (@old_bond_lines) {

    ($ai,     $aj, $funct,   $r,  $k) = split;
    $intern_i = $number2internal{$ai};
    $intern_j = $number2internal{$aj};

    ($oldname_i, $oldnum_i, $newname_i, $newnum_i) = split " ", $internal_rep[$intern_i];
    ($oldname_j, $oldnum_j, $newname_j, $newnum_j) = split " ",$internal_rep[$intern_j];
    
    printf CHG_OFF_ITP  "%6d%5d %4d   %10.4e  %10.4e  %10.4e  %10.4e  ; %s - %s  \n",
    $oldnum_i, $oldnum_j, $funct, $r, $k, $r, $k,  
    $intern_name[$intern_i],  $intern_name[$intern_j];


    printf VDW_ITP  "%6d%5d %4d  ", $intern_i+1, $intern_j+1,  $funct;
    

    if ( $oldname_i eq "-" ||   $oldname_j eq "-" ) {
	printf VDW_ITP  " %10.4e  %10.4e ", $r, 0.0; # how could this happen ?!

    } else {
	printf VDW_ITP  " %10.4e  %10.4e ", $r, $k;
    }

    if ( $newname_i eq "-" ||   $newname_j eq "-" ) {
	printf VDW_ITP  " %10.4e  %10.4e ", $r, 0.0; 

    } else {

	if ( defined $r_new{$newnum_i}{$newnum_j} ) {
	    $r_new = $r_new{$newnum_i}{$newnum_j};
	} elsif  (defined $r_new{$newnum_j}{$newnum_i} ) {
	    $r_new = $r_new{$newnum_j}{$newnum_i};
	} else {
	    die;
	}


	if ( defined $k_new{$newnum_i}{$newnum_j} ) {
	    $k_new = $k_new{$newnum_i}{$newnum_j};
	} elsif  (defined $k_new{$newnum_j}{$newnum_i} ) {
	    $k_new = $k_new{$newnum_j}{$newnum_i};
	} else {
	    die;
	}
	printf VDW_ITP  " %10.4e  %10.4e ", $r_new, $k_new;;

    }
    printf VDW_ITP  " ; %s - %s \n", $intern_name[$intern_i],  $intern_name[$intern_j];

}


foreach  (@new_bond_lines) {

    ($ai,     $aj, $funct,   $r,  $k) = split;

    $is_new = 0;

    $intern_i = $newnumber2internal{$ai};
    ($oldname, $oldnum, $newname, $newnum_i) = split " ", $internal_rep[$intern_i];
    ( $oldname eq "-" )  && ($is_new = 1);


    $intern_j = $newnumber2internal{$aj};
    ($oldname, $oldnum, $newname, $newnum_j) = split " ", $internal_rep[$intern_j];
    ( $oldname eq "-" )  && ($is_new = 1);

    printf CHG_ON_ITP  "%6d%5d %4d   %10.4e  %10.4e  %10.4e  %10.4e  ; %s - %s  \n",
    $newnum_i, $newnum_j, $funct, $r, $k, $r, $k,  
    $intern_name[$intern_i],  $intern_name[$intern_j];

    $is_new || next;


    printf VDW_ITP  "%6d%5d %4d   %10.4e  %10.4e  %10.4e  %10.4e  ; %s - %s  \n",
    $intern_i+1, $intern_j+1, $funct, $r, 0.0, $r, $k,  
    $intern_name[$intern_i],  $intern_name[$intern_j];
	
    
}
##############################
##############################
##############################
##############################
#   pairs - beats me what this crap is
@old_pair_lines = ();
print CHG_OFF_ITP  "\n[ pairs ]\n";
print VDW_ITP      "\n[ pairs ]\n";
print CHG_ON_ITP   "\n[ pairs ]\n";
while ( <ITP_FROM> ) {
    if ( /^\;/) {
	print CHG_OFF_ITP ;
	print VDW_ITP ;
	print CHG_ON_ITP ;
	next;
    }
    next if ( !/\S/);
    last if ( /angles/ );

    push @old_pair_lines, $_;
    
}

@new_pair_lines = ();
while ( <ITP_TO> ) {
    next if ( /^\;/);
    
    next if ( !/\S/);
    last if ( /angles/ );


    push @new_pair_lines, $_;
}


foreach  (@old_pair_lines) {
    ($ai,     $aj, $funct) = split;
    $intern_i = $number2internal{$ai};
    $intern_j = $number2internal{$aj};

    ($oldname_i, $oldnum_i, $newname_i, $newnum_i) = split " ", $internal_rep[$intern_i];
    ($oldname_j, $oldnum_j, $newname_j, $newnum_j) = split " ", $internal_rep[$intern_j];


    printf CHG_OFF_ITP  "%6d%5d %4d   ; %s - %s  \n",
    $oldnum_i, $oldnum_j, $funct,
    $intern_name[$intern_i],  $intern_name[$intern_j];

    printf VDW_ITP  "%6d%5d %4d  ", $intern_i+1, $intern_j+1,  $funct;
    printf VDW_ITP  " ; %s - %s \n", $intern_name[$intern_i],  $intern_name[$intern_j];
}


foreach  (@new_pair_lines) {

    ($ai,     $aj, $funct) = split;

 
    $is_new = 0;

    $intern_i = $newnumber2internal{$ai};
    ($oldname, $oldnum, $newname, $newnum_i) = split " ", $internal_rep[$intern_i];
    ( $oldname eq "-" )  && ($is_new = 1);


    $intern_j = $newnumber2internal{$aj};
    ($oldname, $oldnum, $newname, $newnum_j) = split " ", $internal_rep[$intern_j];
    ( $oldname eq "-" )  && ($is_new = 1);

    printf CHG_ON_ITP  "%6d%5d %4d   ; %s - %s  \n",
    $newnum_i, $newnum_j, $funct,
    $intern_name[$intern_i],  $intern_name[$intern_j];


    $is_new || next;


    printf VDW_ITP  "%6d%5d %4d   ; %s - %s  \n",
    $intern_i+1, $intern_j+1, $funct,
    $intern_name[$intern_i],  $intern_name[$intern_j];
	
    
}
##############################
##############################
##############################
##############################
#   angles
@old_angle_lines = ();
print CHG_OFF_ITP  "\n[ angles ]\n";
print VDW_ITP     "\n[ angles ]\n";
print CHG_ON_ITP  "\n[ angles ]\n";
while ( <ITP_FROM> ) {
    if ( /^\;/) {
	print CHG_OFF_ITP ;
	print VDW_ITP ;
	print CHG_ON_ITP ;
	next;
    }
    next if ( !/\S/);
    last if ( /dihedrals/ );

    push @old_angle_lines, $_;
    
}

@new_angle_lines = ();
while ( <ITP_TO> ) {
    next if ( /^\;/);
    
    next if ( !/\S/);
    last if ( /dihedrals/ );
    ($ai, $aj, $ak, $funct, $theta, $cth) = split;
    $theta_new{$ai}{$aj}{$ak} = $theta_new{$aj}{$ai}{$ak}  = $theta;
    $cth_new{$ai}{$aj}{$ak} = $cth_new{$aj}{$ai}{$ak}  = $cth;

    push @new_angle_lines, $_;
}

foreach  (@old_angle_lines) {
    ($ai, $aj, $ak, $funct, $theta, $cth) = split;

    $intern_i = $number2internal{$ai};
    $intern_j = $number2internal{$aj};
    $intern_k = $number2internal{$ak};

    ($oldname_i, $oldnum_i, $newname_i, $newnum_i) = split " ", $internal_rep[$intern_i];
    ($oldname_j, $oldnum_j, $newname_j, $newnum_j) = split " " ,$internal_rep[$intern_j];
    ($oldname_k, $oldnum_k, $newname_k, $newnum_k) = split " " ,$internal_rep[$intern_k];


    printf CHG_OFF_ITP  "%6d%6d%6d %4d   %10.4e  %10.4e   %10.4e  %10.4e   ; %s - %s  - %s  \n",
    $oldnum_i, $oldnum_j, $oldnum_k, $funct, $theta, $cth, $theta, $cth, 
    $intern_name[$intern_i],  $intern_name[$intern_j],  $intern_name[$intern_k];

    printf VDW_ITP  "%6d%6d%6d %4d  ", $intern_i+1, $intern_j+1, $intern_k+1,  $funct;


    if ( $oldname_i eq "-" ||   $oldname_j eq "-" ||   $oldname_k eq "-" ) {
	printf VDW_ITP  " %10.4e  %10.4e ", $theta, 0.0; # how could this happen ?!

    } else {
	printf VDW_ITP  " %10.4e  %10.4e ", $theta, $cth;
    }

    if ( $newname_i eq "-" ||   $newname_j eq "-" ||   $newname_k eq "-" ) {
	printf VDW_ITP  " %10.4e  %10.4e ", $theta, 0.0; 

    } else {

	if ( defined $theta_new{$newnum_i}{$newnum_j}{$newnum_k} ) {
	    $theta_new = $theta_new{$newnum_i}{$newnum_j}{$newnum_k};
	} else {
	    die;
	}

	if ( defined $cth_new{$newnum_i}{$newnum_j}{$newnum_k} ) {
	    $cth_new = $cth_new{$newnum_i}{$newnum_j}{$newnum_k};
	} else {
	    die;
	}

	printf VDW_ITP  " %10.4e  %10.4e ", $theta_new, $cth_new;;

    }

    printf VDW_ITP  " ; %s - %s  - %s \n", $intern_name[$intern_i],  
    $intern_name[$intern_j],  $intern_name[$intern_k];
}



foreach  (@new_angle_lines) {
    ($ai, $aj, $ak, $funct, $theta, $cth) = split;

 

    $is_new = 0;

    $intern_i = $newnumber2internal{$ai};
    ($oldname, $oldnum, $newname, $newnum_i) = split " ", $internal_rep[$intern_i];
    ( $oldname eq "-" )  && ($is_new = 1);


    $intern_j = $newnumber2internal{$aj};
    ($oldname, $oldnum, $newname, $newnum_j) = split " ", $internal_rep[$intern_j];
    ( $oldname eq "-" )  && ($is_new = 1);

    $intern_k = $newnumber2internal{$ak};
    ($oldname, $oldnum, $newname, $newnum_k) = split " ", $internal_rep[$intern_k];
    ( $oldname eq "-" )  && ($is_new = 1);

    printf CHG_ON_ITP  "%6d%6d%6d %4d   %10.4e  %10.4e   %10.4e  %10.4e   ; %s - %s  - %s  \n",
    $newnum_i, $newnum_j, $newnum_k, $funct, $theta, $cth, $theta, $cth, 
    $intern_name[$intern_i],  $intern_name[$intern_j],  $intern_name[$intern_k];

    $is_new || next;


    printf VDW_ITP  "%6d%6d%6d %4d   %10.4e  %10.4e   %10.4e  %10.4e   ; %s - %s  - %s  \n",
    $intern_i+1, $intern_j+1, $intern_k+1, $funct, $theta, 0.0, $theta, $cth, 
    $intern_name[$intern_i],  $intern_name[$intern_j],  $intern_name[$intern_k];

}



######################################################################
######################################################################
#   dihedrals
@old_dihedral_lines = ();
print CHG_OFF_ITP  "\n[ dihedrals ]\n";
print VDW_ITP  "\n[ dihedrals ]\n";
print CHG_ON_ITP  "\n[ dihedrals ]\n";
while ( <ITP_FROM> ) {
    if ( /^\;/) {
	print CHG_OFF_ITP ;
	print VDW_ITP ;
	print CHG_ON_ITP ;
	next;
    }
    next if ( !/\S/);
    last if ( /dihedrals/ );

    push @old_dihedral_lines, $_;
    
}

@new_dihedral_lines = ();
while ( <ITP_TO> ) {
    next if ( /^\;/);
    
    next if ( !/\S/);
    last if ( /dihedrals/ );
    ($ai, $aj, $ak, $al, $funct, $C0, $C1, $C2, $C3, $C4, $C5) = split;
    @C = ($C0, $C1, $C2, $C3, $C4, $C5);

    @{$C_new{$ai}{$aj}{$ak}{$al}} = @C;

    push @new_dihedral_lines, $_;
}

foreach  (@old_dihedral_lines) {
    ($ai, $aj, $ak, $al, $funct, $C0, $C1, $C2, $C3, $C4, $C5) = split;
     @C = ($C0, $C1, $C2, $C3, $C4, $C5);
   
    $intern_i = $number2internal{$ai};
    $intern_j = $number2internal{$aj};
    $intern_k = $number2internal{$ak};
    $intern_l = $number2internal{$al};

    ($oldname_i, $oldnum_i, $newname_i, $newnum_i) = split " ", $internal_rep[$intern_i];
    ($oldname_j, $oldnum_j, $newname_j, $newnum_j) = split " ",$internal_rep[$intern_j];
    ($oldname_k, $oldnum_k, $newname_k, $newnum_k) = split " ",$internal_rep[$intern_k];
    ($oldname_l, $oldnum_l, $newname_l, $newnum_l) = split " ",$internal_rep[$intern_l];

    printf CHG_OFF_ITP  "%6d%6d%6d%6d %4d   ",
    $oldnum_i, $oldnum_j, $oldnum_k, $oldnum_l, $funct;

    printf CHG_OFF_ITP  "%11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  ", @C;
    printf CHG_OFF_ITP  " %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  ", @C;
    printf CHG_OFF_ITP  " ; %s - %s - %s - %s  \n",
    $intern_name[$intern_i],  $intern_name[$intern_j], $intern_name[$intern_k],  $intern_name[$intern_l];

    printf VDW_ITP  "%6d%6d%6d%6d %4d  ", $intern_i+1, $intern_j+1, $intern_k+1,, $intern_l+1,  $funct;


    if ( $oldname_i eq "-" ||   $oldname_j eq "-" ||   $oldname_k eq "-" ||   $oldname_l eq "-" ) {
	printf VDW_ITP  " %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  ", 
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0; # how could this happen ?!

    } else {
	printf VDW_ITP   " %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  ",  @C;
    }

    if ( $newname_i eq "-" ||   $newname_j eq "-" ||   $newname_k eq "-" ||   $newname_l eq "-" ) {
	printf VDW_ITP   " %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  ", 
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    } else {

	if (  defined $C_new{$newnum_i}{$newnum_j}{$newnum_k}{$newnum_l} ) {
	    printf VDW_ITP   " %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  ", 
	    @{$C_new{$newnum_i}{$newnum_j}{$newnum_k}{$newnum_l}};

	} elsif (  defined $C_new{$newnum_l}{$newnum_k}{$newnum_j}{$newnum_i} ) {
	    printf VDW_ITP   " %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  ", 
	    @{$C_new{$newnum_l}{$newnum_k}{$newnum_j}{$newnum_i}};

	} else {
	    die "C new not defined for $newnum_i  $newnum_j  $newnum_k  $newnum_l ";
	}
    }



    printf VDW_ITP  " ; %s - %s - %s - %s \n", $intern_name[$intern_i],  $intern_name[$intern_j],
    $intern_name[$intern_k], $intern_name[$intern_l] ;


 
}

foreach  (@new_dihedral_lines) {

    ($ai, $aj, $ak, $al, $funct, $C0, $C1, $C2, $C3, $C4, $C5) = split;
     @C = ($C0, $C1, $C2, $C3, $C4, $C5);
 

    $is_new = 0;

    $intern_i = $newnumber2internal{$ai};
    ($oldname, $oldnum, $newname, $newnum_i) = split " ", $internal_rep[$intern_i];
    ( $oldname eq "-" )  && ($is_new = 1);


    $intern_j = $newnumber2internal{$aj};
    ($oldname, $oldnum, $newname, $newnum_j) = split " ", $internal_rep[$intern_j];
    ( $oldname eq "-" )  && ($is_new = 1);

    $intern_k = $newnumber2internal{$ak};
    ($oldname, $oldnum, $newname, $newnum_k) = split " ", $internal_rep[$intern_k];
    ( $oldname eq "-" )  && ($is_new = 1);

    $intern_l = $newnumber2internal{$al};
    ($oldname, $oldnum, $newname, $newnum_l) = split " ", $internal_rep[$intern_l];
    ( $oldname eq "-" )  && ($is_new = 1);

    printf CHG_ON_ITP  "%6d%6d%6d%6d %4d   ",
    $newnum_i, $newnum_j, $newnum_k, $newnum_l, $funct;

    printf CHG_ON_ITP  "%11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  ", @C;
    printf CHG_ON_ITP  " %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  ", @C;
    printf CHG_ON_ITP  " ; %s - %s - %s - %s  \n",
    $intern_name[$intern_i],  $intern_name[$intern_j], $intern_name[$intern_k],  $intern_name[$intern_l];

    $is_new || next;


    printf VDW_ITP  "%6d%6d%6d%6d %4d   ",
    $intern_i+1, $intern_j+1, $intern_k+1, $intern_l+1, $funct;

    printf VDW_ITP   "%11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  ", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    printf VDW_ITP   " %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  ", @C;

    printf VDW_ITP  " ; %s - %s - %s - %s  \n",
    $intern_name[$intern_i],  $intern_name[$intern_j], $intern_name[$intern_k],  $intern_name[$intern_l];
	
    
}

######################################################################
######################################################################
# one last round of "dihedrals"
#   dihedrals
@old_dihedral_lines = ();
print CHG_OFF_ITP  "\n[ dihedrals ]\n";
print VDW_ITP  "\n[ dihedrals ]\n";
print CHG_ON_ITP  "\n[ dihedrals ]\n";
while ( <ITP_FROM> ) {
    if ( /^\;/) {
	print CHG_OFF_ITP ;
	print VDW_ITP ;
	print CHG_ON_ITP ;
	next;
    }
    next if ( !/\S/);
    last if ( /dihedrals/ );

    push @old_dihedral_lines, $_;
    
}

@new_dihedral_lines = ();
while ( <ITP_TO> ) {
    next if ( /^\;/);
    
    next if ( !/\S/);
    last if ( /dihedrals/ );
    ($ai, $aj, $ak, $al, $funct, $phase, $kd, $pn) = split;
 
    $phase_new{$ai}{$aj}{$ak}{$al} = $phase;
    $kd_new{$ai}{$aj}{$ak}{$al}    = $kd;
    $pn_new{$ai}{$aj}{$ak}{$al}    = $pn;

    push @new_dihedral_lines, $_;
}

foreach  (@old_dihedral_lines) {
    ($ai, $aj, $ak, $al, $funct, $phase, $kd, $pn) = split;
   
    $intern_i = $number2internal{$ai};
    $intern_j = $number2internal{$aj};
    $intern_k = $number2internal{$ak};
    $intern_l = $number2internal{$al};

    ($oldname_i, $oldnum_i, $newname_i, $newnum_i) = split " ", $internal_rep[$intern_i];
    ($oldname_j, $oldnum_j, $newname_j, $newnum_j) = split " ",$internal_rep[$intern_j];
    ($oldname_k, $oldnum_k, $newname_k, $newnum_k) = split " ",$internal_rep[$intern_k];
    ($oldname_l, $oldnum_l, $newname_l, $newnum_l) = split " ",$internal_rep[$intern_l];

    printf CHG_OFF_ITP  "%6d%6d%6d%6d %4d   ",
    $oldnum_i, $oldnum_j, $oldnum_k, $oldnum_l, $funct;
    printf CHG_OFF_ITP   "%11.5f  %11.5f  %3d ",  $phase, $kd, $pn;
    printf CHG_OFF_ITP   "%11.5f  %11.5f  %3d ",  $phase, $kd, $pn;
    printf CHG_OFF_ITP  " ; %s - %s - %s - %s  \n",
    $intern_name[$intern_i],  $intern_name[$intern_j], $intern_name[$intern_k],  $intern_name[$intern_l];


    printf VDW_ITP  "%6d%6d%6d%6d %4d  ", $intern_i+1, $intern_j+1, $intern_k+1,, $intern_l+1,  $funct;


    if ( $oldname_i eq "-" ||   $oldname_j eq "-" ||   $oldname_k eq "-" ||   $oldname_l eq "-" ) {
	printf VDW_ITP  " %11.5f  %11.5f %3d ", 
	$phase, 0.0, $pn; # how could this happen ?!

    } else {
	printf VDW_ITP    " %11.5f  %11.5f %3d ", $phase, $kd, $pn;
    }

    if ( $newname_i eq "-" ||   $newname_j eq "-" ||   $newname_k eq "-" ||   $newname_l eq "-" ) {
	printf VDW_ITP  " %11.5f  %11.5f %3d ", 
	$phase, 0.0,$pn; 

    } else {

	if ( defined $phase_new{$newnum_i}{$newnum_j}{$newnum_k}{$newnum_l} ) {
	    printf VDW_ITP   " %11.5f  ", $phase_new{$newnum_i}{$newnum_j}{$newnum_k}{$newnum_l};
	} elsif ( defined $phase_new{$newnum_l}{$newnum_k}{$newnum_j}{$newnum_i} ) {
	    printf VDW_ITP   " %11.5f  ", $phase_new{$newnum_l}{$newnum_k}{$newnum_j}{$newnum_i};
	} else { 
	    die;
	}


	if ( defined $kd_new{$newnum_i}{$newnum_j}{$newnum_k}{$newnum_l} ) {
	    printf VDW_ITP   " %11.5f  ", $kd_new{$newnum_i}{$newnum_j}{$newnum_k}{$newnum_l};
	} elsif ( defined $kd_new{$newnum_l}{$newnum_k}{$newnum_j}{$newnum_i} ) {
	    printf VDW_ITP   " %11.5f  ", $kd_new{$newnum_l}{$newnum_k}{$newnum_j}{$newnum_i};
	} else { 
	    die;
	}

	if ( defined $pn_new{$newnum_i}{$newnum_j}{$newnum_k}{$newnum_l} ) {
	    printf VDW_ITP   " %3d  ", $pn_new{$newnum_i}{$newnum_j}{$newnum_k}{$newnum_l};
	} elsif ( defined $pn_new{$newnum_l}{$newnum_k}{$newnum_j}{$newnum_i} ) {
	    printf VDW_ITP   " %3d  ", $pn_new{$newnum_l}{$newnum_k}{$newnum_j}{$newnum_i};
	} else { 
	    die;
	}



    }



    printf VDW_ITP  " ; %s - %s - %s - %s \n", $intern_name[$intern_i],  $intern_name[$intern_j],
    $intern_name[$intern_k], $intern_name[$intern_l] ;


 
}


foreach  (@new_dihedral_lines) {

    ($ai, $aj, $ak, $al, $funct, $phase, $kd, $pn) = split;

 
    $is_new = 0;

    $intern_i = $newnumber2internal{$ai};
    ($oldname, $oldnum, $newname, $newnum_i) = split " ", $internal_rep[$intern_i];
    ( $oldname eq "-" )  && ($is_new = 1);


    $intern_j = $newnumber2internal{$aj};
    ($oldname, $oldnum, $newname, $newnum_j) = split " ", $internal_rep[$intern_j];
    ( $oldname eq "-" )  && ($is_new = 1);

    $intern_k = $newnumber2internal{$ak};
    ($oldname, $oldnum, $newname, $newnum_k) = split " ", $internal_rep[$intern_k];
    ( $oldname eq "-" )  && ($is_new = 1);

    $intern_l = $newnumber2internal{$al};
    ($oldname, $oldnum, $newname, $newnum_l) = split " ", $internal_rep[$intern_l];
    ( $oldname eq "-" )  && ($is_new = 1);


    printf CHG_ON_ITP  "%6d%6d%6d%6d %4d   ",
    $newnum_i, $newnum_j, $newnum_k, $newnum_l, $funct;
    printf CHG_ON_ITP   "%11.5f  %11.5f  %3d ",  $phase, $kd, $pn;
    printf CHG_ON_ITP   "%11.5f  %11.5f  %3d ",  $phase, $kd, $pn;
    printf CHG_ON_ITP  " ; %s - %s - %s - %s  \n",
    $intern_name[$intern_i],  $intern_name[$intern_j], $intern_name[$intern_k],  $intern_name[$intern_l];


    $is_new || next;


    printf VDW_ITP  "%6d%6d%6d%6d %4d   ",
    $intern_i+1, $intern_j+1, $intern_k+1, $intern_l+1, $funct;

    printf VDW_ITP   "%11.5f  %11.5f  %3d ",  $phase, 0.0, $pn;

    printf VDW_ITP   "%11.5f  %11.5f  %3d ",  $phase, $kd, $pn;

    printf VDW_ITP  " ; %s - %s - %s - %s  \n",
    $intern_name[$intern_i],  $intern_name[$intern_j], $intern_name[$intern_k],  $intern_name[$intern_l];
	
   
}

close CHG_OFF_ITP;
close VDW_ITP;
close CHG_ON_ITP;

######################################################################
# make new gro to go with this


######################################################################
open ( CHG_OFF_GRO, ">$output.chg_off.gro") ||
    die "Cno $output.itp: $!\n";

printf CHG_OFF_GRO "; blah blah\n";
printf CHG_OFF_GRO "%5d \n", scalar @atom_names;

foreach $old_atom_name ( @atom_names) {

    $intern = $name2internal{$old_atom_name};

    ($oldname, $oldnum, $newname, $newnum_i) = split " ", $internal_rep[$intern];


    ($oldname eq "-") && next;

    $x = $x{$oldname};
    $y = $y{$oldname};
    $z = $z{$oldname};
    
    printf CHG_OFF_GRO "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
    1, $res_name, $intern_name[$intern], $oldnum_i, $x, $y, $z; 

}

printf CHG_OFF_GRO $box_line;
close  CHG_OFF_GRO;
######################################################################

open ( VDW_GRO, ">$output.vdw.gro") ||
    die "Cno $output.itp: $!\n";

printf VDW_GRO "; blah blah\n";
printf VDW_GRO "%5d \n", scalar @internal_rep;

foreach $intern ( 0 .. $#internal_rep) {
    ($oldname, $oldnum, $newname, $newnum) = split " ", $internal_rep[$intern];
    if ( defined $x{$oldname} ) {
	$x = $x{$oldname};
	$y = $y{$oldname};
	$z = $z{$oldname};

    } elsif( defined $new_x{$newname}) {
	$x = $new_x{$newname};
	$y = $new_y{$newname};
	$z = $new_z{$newname};

    } else {

	die "$intern  $oldname  $oldnum,  $newname $newnum\n" ;
    }

    printf VDW_GRO "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
    1, $res_name, $intern_name[$intern], $intern+1, $x, $y, $z; 
}
printf VDW_GRO $box_line;
close VDW_GRO;


######################################################################
open ( CHG_ON_GRO, ">$output.chg_on.gro") ||
    die "Cno $output.itp: $!\n";

printf CHG_ON_GRO "; blah blah\n";
printf CHG_ON_GRO "%5d \n", scalar @new_atom_names;

foreach $new_atom_name ( @new_atom_names) {

    $intern = $newname2internal{$new_atom_name};

    ($oldname, $oldnum, $newname, $newnum_i) = split " ", $internal_rep[$intern];


    ($newname eq "-") && next;

    $x = $new_x{$newname};
    $y = $new_y{$newname};
    $z = $new_z{$newname};
    
    printf CHG_ON_GRO "%5d%5s%5s%5d%8.3f%8.3f%8.3f\n",
    1, $res_name, $intern_name[$intern], $newnum_i, $x, $y, $z; 

}

printf CHG_ON_GRO $box_line;
close  CHG_ON_GRO;



######################################################################
######################################################################
######################################################################

close ITP_FROM;
close ITP_TO;

