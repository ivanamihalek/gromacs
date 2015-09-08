#! /usr/bin/perl -w

use Math::Trig;

sub unit   (@);
sub vect_product (@);
sub affine (@);
sub quat_to_rotation (@);
sub line2coords (@);
sub dot_product (@);

$dot_tol = 1.e-4;


@ARGV ||
   die "Usage:  $0  <compound pdb>  <group pdb> \n";

($cmpd, $group) = @ARGV;

foreach ( $cmpd, $group) {
   ( -e $_) || die  "$_ not found\n";
}

$filename = $cmpd;
open (IF, "<$filename" )
   || die "Cno $filename: $!.\n";
@cmpd_carbs     = ();
@cmpd_hydrogens = ();
@cmpd           = ();

$carb_id = 0;
$h_id = 0;
$line_ctr = 0;
while ( <IF> ) {

   next if ( ! /^ATOM/ && ! /^HETATM/ );
   push @cmpd, $_;

   $name = substr $_,  12, 4 ;  $name =~ s/\s//g;

   if ( $name =~ /^C/ ) {
       ($x, $y, $z) = line2coords ($_);
       @{$cmpd_carbs[$carb_id]} = ($x, $y, $z);
       $carb_id ++;
       #$cmpd_carb_line[$carb_id] = $line_ctr;

   } elsif ( $name =~ /^H/ ) {
       ($x, $y, $z) = line2coords ($_);
       @{$cmpd_hydrogens[$h_id]} =  ($x, $y, $z);
       $cmpd_line[$h_id] = $line_ctr;
       $h_id++;

   }
   $line_ctr++;

}
close IF;

foreach $carb_id ( 0 ..  $#cmpd_carbs ) {

   @c = @{$cmpd_carbs[$carb_id]};

   foreach $h_id ( 0 .. $#cmpd_hydrogens ) {
       @h = @{$cmpd_hydrogens[$h_id]};
       $d =  0;
       for $i  ( 0 ..2 ) {
           $aux = $c[$i] - $h[$i];
           $d  += $aux*$aux;
       }
       $d = sqrt ($d);
       if ( $d < 1.2 ) {
           push @{$neighbor[$carb_id]}, $h_id;
       }
   }
}

##############################################
# now the same for the group
$filename = $group;
open (IF, "<$filename" )
   || die "Cno $filename: $!.\n";
@group_carbs = ();
@group_hydrogens = ();

$carb_id =  0;
$h_id    =  0;
@group   = ();
$line_ctr = 0;
while ( <IF> ) {
   next if ( ! /^ATOM/ && ! /^HETATM/ );
   push @group, $_;

   $name = substr $_,  12, 4 ;  $name =~ s/\s//g;

 
   if ( $name =~ /^C/ ) {
       ($x, $y, $z) = line2coords ($_);
       @{$group_carbs[$carb_id]}  = ($x, $y, $z);
       #$group_carb_line[$carb_id] = $line_ctr;
       $carb_id ++;

   } elsif ( $name =~ /^H/ ) {
       ($x, $y, $z) = line2coords ($_);
       @{$group_hydrogens[$h_id]} = ($x, $y, $z);
       $group_line[$h_id] = $line_ctr;
       $h_id++;

   }

   $line_ctr++;

}
close IF;

foreach $carb_id ( 0 ..  $#group_carbs ) {

   @c = @{$group_carbs[$carb_id]};

   foreach $h_id ( 0 .. $#group_hydrogens ) {
       @h = @{$group_hydrogens[$h_id]};
       $d =  0;
       for $i  ( 0 ..2 ) {
           $aux = $c[$i] - $h[$i];
           $d  += $aux*$aux;
       }
       $d = sqrt ($d);
       if ( $d < 1.2 ) {
           push @{$group_neighbor[$carb_id]}, $h_id;
       }
   }
}



##############################################
#
@sockets = ();
foreach $carb_id ( 0 ..  $#cmpd_carbs ) {
   defined $neighbor[$carb_id] || next;

   if ( scalar  @{$neighbor[$carb_id]} == 1 ||
        scalar  @{$neighbor[$carb_id]} >= 3 ) {
       $h_id = $neighbor[$carb_id][0];
       push @sockets, "$carb_id $h_id";
   } else {
       foreach $h_id  ( @{$neighbor[$carb_id]} ) {
           push @sockets, "$carb_id $h_id";
       }

   }

}

#print join "\n", @sockets;


@plugs = ();
foreach $carb_id ( 0 ..  $#group_carbs ) {
   defined $group_neighbor[$carb_id] || next;

   if ( scalar  @{$group_neighbor[$carb_id]} == 1 ||
        scalar  @{$group_neighbor[$carb_id]} >= 3 ) {
       $h_id = $group_neighbor[$carb_id][0];
       push @plugs, "$carb_id $h_id";
   } else {
       foreach $h_id  ( @{$group_neighbor[$carb_id]} ) {
           push @plugs, "$carb_id $h_id";
       }

   }

}

#print join "\n", @plugs;



##############################################
# plug each plug into each socket and output
# the coordinates

$cmpd_no = 0;
@R = ();

foreach $socket (@sockets) {
   ($carb_socket, $h_socket) = split " ", $socket;
   # new position for the $carb from the plug
   @socket_carb    = @{$cmpd_carbs[$carb_socket]};
   @sock_h_directn     = unit (@socket_carb,  @{$cmpd_hydrogens[$h_socket]});
   @neg_sock_h_directn = map { -$_} @sock_h_directn;

   for $i ( 0 .. 2) {
       $new_plug_carb[$i] = $socket_carb[$i]+1.5*$sock_h_directn[$i];
   }

   foreach $plug (@plugs) {

       $cmpd_no ++;
       
 
       ($carb_plug, $h_plug) = split " ", $plug;
       @plug_carb      =  @{$group_carbs[$carb_plug]};
       @plug_hydrogen  =  @{$group_hydrogens[$h_plug]};
       @plug_h_directn =  unit (@plug_carb,  @plug_hydrogen);

       
       # vector product
       @R = ();

       # first check whther we are alreayd (anti) parallel
       $dot = dot_product(@plug_h_directn, @neg_sock_h_directn);
       if ( abs ($dot - 1) < $dot_tol ) { # cos of 1 deg is approx 1-1.e-4
	   # we are already in pretty much the same direction
	   @R =( (1,0,0), (0,1,0), (0,0,1) );

       } elsif ( abs ($dot + 1) < $dot_tol) {

	   # rotate by 180
 	   # find any vector perp to plug_carb and rotate by pi
	   if ( abs($plug_h_directn[0]-1) < $dot_tol ) {
	       @rot_axis = (0,1,0);
	   } elsif ( abs($plug_h_directn[1]-1) < $dot_tol ) {
	       @rot_axis = (0,0,1);
	   } elsif ( abs($plug_h_directn[2]-1) < $dot_tol ) {
	       @rot_axis = (1,0,0);
	   } else  {
	       # any non-colinear
	       $largest_comp = -1;
	       $largest_val = -1;
	       for $i ( 0 ..2) {
		   if ( abs($plug_h_directn[$i]) > $largest_val) {
		       $largest_val  = abs($plug_h_directn[$i]);
		       $largest_comp = $i;
		   } 
	       }
	       @non_col = @plug_h_directn;
	       $non_col[$largest_comp] = 0;

	       # take out the comp along plug_h_directn
	       $dot =  dot_product(@plug_h_directn, @non_col);
	       for $i ( 0 ..2) {
		   $perp[$i] = $non_col[$i] - $dot*$plug_h_directn[$i];
	       }

	       # normalize
	       $dot  = dot_product (@perp, @perp);
	       $dot  = sqrt ($dot);
	       @perp = map { $_/$dot } @perp;

	       ($angle, @rot_axis) = (3.24, @perp);

	       # turn into rot matrix
	       $half  = $angle/2;
	       $shalf = sin ($half);
	       @quat  = (cos($half), $rot_axis[0]*$shalf, 
			 $rot_axis[1]*$shalf,
			 $rot_axis[2]*$shalf);
	       @R = ();
	       quat_to_rotation (\@R, @quat);
	   }

       } else {

	   ($angle, @rot_axis) = vect_product (@plug_h_directn, @neg_sock_h_directn);
	   # turn into rot matrix
	   $half  = $angle/2;
	   $shalf = sin ($half);
	   @quat = ( cos($half),$rot_axis[0]*$shalf, $rot_axis[1]*$shalf,
		     $rot_axis[2]*$shalf);


	   @R = ();
	   quat_to_rotation (\@R, @quat);

       }
       #exit;
       # affine tfm
       @new_group = affine (\@R, @plug_carb, @new_plug_carb, @group);
       #output coordinates, skipping both plug and socket hydrogens
 
       $no_atoms_in_cmpd = 0;
       $filename = "cmpd$cmpd_no.pdb";
       open (OF, ">$filename" )
	   || die "Cno $filename: $!.\n";

       $line_ctr = 0;
       foreach $line ( @cmpd ) {
	   if ($line_ctr != $cmpd_line[$h_socket]) {
	       (substr $line,  0, 6) = "HETATM";
	       (substr $line, 17, 3) = "LIG";
	       print OF $line;
	   }
 	   $line_ctr ++;
      }
       $no_atoms_in_cmpd = $line_ctr;

       $atom_id  = $no_atoms_in_cmpd + 1;
       $line_ctr = 0;
       foreach  $line ( @new_group) {
	   if ($line_ctr != $group_line[$h_plug]) {
	       (substr $line,  0, 6) = "HETATM";
	       (substr $line,  6, 5) = sprintf "%5d", $atom_id;
	       (substr $line, 17, 3) = "LIG";
	       (substr $line, 20, 1) = " ";
	       print OF $line;
	       $atom_id++;
	   }
 
  	   $line_ctr ++;
       }

       close OF;
 
       next;
  }
}




##############################################
##############################################
##############################################
sub unit (@) {
   my @from = @_[0..2];
   my @to   = @_[3..5];
   my @unit;
   my $norm = 0;
   my $i;

   for $i ( 0 .. 2) {
       $unit[$i] = $to[$i]- $from[$i];
       $norm += $unit[$i]*$unit[$i];
   }
   $norm = sqrt ($norm);

   for $i ( 0 .. 2) {
       $unit[$i] /= $norm;
   }

   return  @unit;
}


##############################################
sub affine (@) {
    my $R_ref = shift @_;
    my @tr1 = @_[0 ..2];
    shift @_;
    shift @_;
    shift @_;
    my @tr2 = @_[0 ..2];
    shift @_;
    shift @_;
    shift @_;
    my @in_line = @_;
    my ($crap, $crap2, $x, $y, $z, $new_line);
    my @out_line = ();

    foreach  (  @in_line ) {
	if ( ! ( /^ATOM/ || /^HETATM/) ){
	    push @out_line, $_;
	    next;
	}
	chomp;
	$crap = substr ($_, 0, 30);
	$crap2 = substr ($_, 54);
	$x = substr $_, 30, 8;  $x=~ s/\s//g;
	$y = substr $_, 38, 8;  $y=~ s/\s//g;
	$z = substr $_, 46, 8;  $z=~ s/\s//g;

	$xtr  = $x;
	$ytr  = $y;
	$ztr  = $z;

	$xtr  -= $tr1[0];
	$ytr  -= $tr1[1];
	$ztr  -= $tr1[2];
	

	# rotate
        $xnew = $$R_ref[0][0]*$xtr +   $$R_ref[0][1]*$ytr  +  $$R_ref[0][2]*$ztr;
	$ynew = $$R_ref[1][0]*$xtr +   $$R_ref[1][1]*$ytr  +  $$R_ref[1][2]*$ztr;
	$znew = $$R_ref[2][0]*$xtr +   $$R_ref[2][1]*$ytr  +  $$R_ref[2][2]*$ztr;


	$xnew += $tr2[0];
	$ynew += $tr2[1];
	$znew += $tr2[2];


	$new_line = sprintf  "%30s%8.3f%8.3f%8.3f%s \n",
	$crap,  $xnew, $ynew, $znew, $crap2;

	push @out_line, $new_line;

    }


    return @out_line;

}

############################################33
sub dot_product (@) {
   my @x = @_[0..2];
   my @y = @_[3..5];
   my $dot = 0;

   for ($i=0; $i<3; $i++ ) {
       $dot += $x[$i]*$y[$i];
   }
   return $dot;
}

############################################33
sub vect_product (@) {

   my @x = @_[0..2];
   my @y = @_[3..5];
   my @v;
   my $norm = 0;
   my $theta;

   $v[0]  = $x[1]*$y[2] - $x[2]*$y[1];
   $norm += $v[0]*$v[0];

   $v[1]  = $x[2]*$y[0] - $x[0]*$y[2];
   $norm += $v[1]*$v[1];

   $v[2]  = $x[0]*$y[1] - $x[1]*$y[0];
   $norm += $v[2]*$v[2];

   ($norm) || return (0,0,0,0);

   $norm  = sqrt ($norm);

   for ($i=0; $i<3; $i++ ) {
       $v[$i] /= $norm;
   }
   
   if ( dot_product ( @x,@y ) > 0 ) {
       $theta = asin($norm);
   } else {
       $theta = 3.14 - asin($norm);
   }
   #norm = nor_x*norm_y*sin theta - we are assuming x and y normalized
   return ( $theta, @v);
}


############################################33
sub  quat_to_rotation (@) {

    my $R_ref = shift @_;
    my @q =  @_;

   $$R_ref[0][0] = 1 -2*$q[2]*$q[2] -2*$q[3]*$q[3];
   $$R_ref[1][1] = 1 -2*$q[3]*$q[3] -2*$q[1]*$q[1] ;
   $$R_ref[2][2] = 1 -2*$q[1]*$q[1] -2*$q[2]*$q[2] ;

   $$R_ref[0][1] = 2*$q[1]*$q[2] - 2*$q[0]*$q[3];
   $$R_ref[1][0] = 2*$q[1]*$q[2] + 2*$q[0]*$q[3];

   $$R_ref[0][2] = 2*$q[1]*$q[3] + 2*$q[0]*$q[2];
   $$R_ref[2][0] = 2*$q[1]*$q[3] - 2*$q[0]*$q[2];

   $$R_ref[1][2] = 2*$q[2]*$q[3] - 2*$q[0]*$q[1];
   $$R_ref[2][1] = 2*$q[2]*$q[3] + 2*$q[0]*$q[1];


}
############################################33
sub line2coords (@) {

    my $line = shift @_;
    my ($x, $y, $z);

    $x = substr $line,30, 8;  $x=~ s/\s//g;
    $y = substr $line,38, 8;  $y=~ s/\s//g;
    $z = substr $line,46, 8;  $z=~ s/\s//g;

    return ($x, $y, $z);
}
###################################
# junkyard
=pod
	   #printf " %8.3f  %8.3f  %8.3f  %8.3f  %8.3f \n", $angle, @quat; 
	   printf "sock h direction      %8.3f  %8.3f  %8.3f \n", @neg_sock_h_directn; 
	   printf " old plug h direction %8.3f  %8.3f  %8.3f \n", @plug_h_directn; 
	   printf "dot:  %8.3f \n", dot_product(@plug_h_directn, @neg_sock_h_directn);
	   for $i ( 0 .. 2) {
	       $test[$i] = 0;
	       for $j ( 0 .. 2) {
		   $test[$i] += $R[$i][$j]*$plug_h_directn[$j];
	       }
	   }
	   printf "               test   %8.3f  %8.3f  %8.3f \n", @test; 

	   #printf " %8.3f  %8.3f  %8.3f \n", @{$R[0]};
	   #printf " %8.3f  %8.3f  %8.3f \n", @{$R[1]};
	   #printf " %8.3f  %8.3f  %8.3f \n", @{$R[2]};


       @plug_carb      =  line2coords ($group[$group_carb_line[$carb_plug]]);
       @plug_hydrogen  =  line2coords ($group[$group_line[$h_plug]]);

       @plug_h_directn =  unit (@plug_carb,  @plug_hydrogen);
       printf " old plug h direction %8.3f  %8.3f  %8.3f \n", 
       @plug_h_directn; 
	   
 
       @plug_carb      =  line2coords ($new_group[$group_carb_line[$carb_plug]]);
       @plug_hydrogen  =  line2coords ($new_group[$group_line[$h_plug]]);

       @new_plug_h_directn =  unit (@plug_carb,  @plug_hydrogen);
       printf " new plug h direction %8.3f  %8.3f  %8.3f \n", 
       @new_plug_h_directn; 
	   
=cut
