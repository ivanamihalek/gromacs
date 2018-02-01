#! /usr/bin/perl -w

sub output_frame ();
sub output_last_frame ();

(@ARGV >= 2) ||
    die "Usage:  $0  <games stddout>  <output name (root)>\n";

($filename, $pdb) = @ARGV;
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";


# first find the atom types from
# the top of the file:
%atom_name = ();
$reading = 0;
while ( <IF> ) {
    last if ( /INTERNUCLEAR DISTANCES/);
    if ( /ATOM      ATOMIC                      COORDINATES (BOHR)/) {
	$reading = 1;
    } elsif ($reading) {
	next if ( !/\S/);
	next if ( !/CHARGE/);
	($atom_name, $charge) = split " ", $_;
	$atom_name{$charge} = $atom_name;
    }
}

@lines =();
$reading = 0;
$frame_no = 0;

open (OF, ">$pdb\_traj.pdb") ||
    die "Cno  $pdb\_traj.pdb: $!\n";

while ( <IF> ) {
    if  ( /ATOM   CHARGE       X              Y              Z/) {
	$reading = 1;
    } elsif ($reading) {
	
	next if (/\-\-\-\-/);

	if ( !/\S/ ) {
	    output_frame ();
	    @last_frame = @lines;
	    @lines =();
	    $reading = 0;
	} else {
	    push @lines,  $_;
	}
    }

}

close OF;
# find partial charges
seek IF, 0,0;

$ctr = 0;
$reading = 0;
while ( <IF> ) {
    if ( /ATOM         MULL\.POP\.    CHARGE          LOW\.POP\.     CHARGE/ ) {
	$ctr++;
	($ctr == 2) && ($reading = 1);
    } elsif ( $reading ) {
	last if ( !/\S/ );
	@aux = split;
	$mull_charge{$aux[0]} = $aux[3];
    }
}
close IF;


#######################################
# output the last frame withe the B-factor hacked
@lines = @last_frame;
open (OF, ">$pdb\_charges.pdb") ||
    die "Cno  $pdb\_traj.pdb: $!\n";

output_last_frame (); 

close OF;

#######################################
# output the last frame withe the B-factor hacked


#######################################

sub output_frame () {

    my $line;
    $frame_no ++;

    $record = "HETATM";
    $serial = 0;
    $atom_name = "";
    $alt_loc = " ";
    $res_name = "LIG";
    $chain_id = "A";  
    $res_ctr  = "1";
    $i_code   = " ";   
    printf  OF  "MODEL %3d\n", $frame_no;
    foreach $line (@lines ) {
       $serial++;
       chomp;
       ($atom_name, $charge, $x, $y, $z) = split  " ", $line;
       $atom_name .= $serial;
       printf  OF "%-6s%5d  %-3s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f\n", 
       $record,   $serial,  $atom_name,   $alt_loc,   $res_name,
       $chain_id,  $res_ctr ,   $i_code ,   $x,  $y,   $z;
    }
    print OF "ENDMDL\n";
}
#######################################

sub output_last_frame () {

    my $line;
    $frame_no ++;

    $record = "HETATM";
    $serial = 0;
    $atom_name = "";
    $alt_loc = " ";
    $res_name = "LIG";
    $chain_id = "A";  
    $res_ctr  = "1";
    $i_code   = " ";   
    foreach $line (@lines ) {
       $serial++;
       chomp;
       ($atom_name, $charge, $x, $y, $z) = split  " ", $line;
       $atom_name .= $serial;
       printf  OF "%-6s%5d  %-3s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f\n", 
       $record,   $serial,  $atom_name,   $alt_loc,   $res_name,
       $chain_id,  $res_ctr ,   $i_code ,   $x,  $y,   $z, 1.0, $mull_charge{$serial};
    }
    print OF "END\n";
}
