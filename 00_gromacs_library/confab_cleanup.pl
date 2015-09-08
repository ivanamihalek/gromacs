#! /usr/bin/perl -w

sub process ();

(@ARGV ==2) ||
    die "Usage:  $0  <template sqm.in>  <template pdb> \n";

$confab = "/usr/local/bin/confab";
$pdb2gro = "/home/ivanam/perlscr/translation/pdb2gro.pl";

foreach ( @ARGV, $confab,  $pdb2gro) {
    (-e $_) || die "$_ not found.\n";
}

($template_coords, $template_pdb) = @ARGV;

###################################

$filename = $template_coords; # template header & coords
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";

@header = ();
@template_coords = ();
$reading_hdr = 1;
while ( <IF> ) {
    next if ( !/\S/);
    if ( $reading_hdr ) {
	push @header, $_;
	/\// && ($reading_hdr = 0);
    } else {
	push @template_coords, $_;
    }
}

close IF;

print @header;

###################################
$confab_out   = "confab.pdb";;
$confab_cmd = "$confab $template_pdb $confab_out";

foreach $param_choice (  " -c 100 -r 0.1",  " -c 1000 -r 0.3") {

    (system $confab_cmd.$param_choice) && die "Error running $cmd\n";


    ###################################

    $out_ctr = 0;

    $filename = $confab_out; # pdb file
    open (IF, "<$filename" ) 
	|| die "Cno $filename: $!.\n";

    @conformer = ();

    while ( <IF> ) {
	next if ( !/\S/);

	if ( /^END/ ) {
	    process ();
	    @conformer = ();
	} elsif ( /^ATOM/ || /^HETATM/ ) {
	    push  @conformer, $_;

	}
    }
    close IF;


}


###################################

sub process () {

    my $line_ctr = 0;
    my @new_cords = ();
    my @new_sqm = @header;
    my @conformer_matching_lines = ();

    foreach $line_ctr  ( 0 .. $#template_coords ) {

	$template_line   = $template_coords[$line_ctr];
	$conformer_line = $conformer[$line_ctr];

	$template_type = substr $template_line, 6, 5;
	$template_type =~ s/\s//g;

	$conformer_type = substr $conformer_line, 12, 4;
	$conformer_type =~ s/\s//g;

	if ( substr ($conformer_type, 0, 1)  ne substr ($template_type, 0, 1)) {
	    print "type mismatch:\n";
	    print "template_line $template_line";
	    print "conformer_line $conformer_line";
	    die "\n";
	}

	# insert new coords into template
	$x = substr $conformer_line, 30, 8;
	$y = substr $conformer_line, 38, 8;
	$z = substr $conformer_line, 46, 8;


	$new_x = sprintf "%10.4lf", $x;
	substr ($template_line, 13, 10) = $new_x;

	$new_y = sprintf "%10.4lf", $y;
	substr ($template_line, 28, 10) = $new_y;

	$new_z = sprintf "%10.4lf", $z;
	substr ($template_line, 42, 10) = $new_z;

	push  @new_sqm, $template_line;
	
	push @conformer_matching_lines, $conformer_line;
    }

    $out_ctr++;

    $newpdb = "confab.$out_ctr.pdb";
    open (NEWPDB, ">$newpdb") || die "Cannot open $newpdb: $!\n";
    print NEWPDB @conformer_matching_lines;
    close NEWPDB;

    

    ($out_ctr >= 10 ) && exit;

}
