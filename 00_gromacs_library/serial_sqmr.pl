#! /usr/bin/perl -w

sub process ();

(@ARGV ==2) ||
    die "Usage:  $0  <template sqm.in>  <template pdb> \n";

$confab = "/usr/local/bin/confab";
$sqm    = "/home/ivanam/docking/02_ligand_prep/amber11/bin/sqm";
$acpype = "/home/ivanam/docking/02_ligand_prep/acpype/acpype.py";
$pdb2gro = "/home/ivanam/perlscr/gromacs/pdb2gro.pl";

foreach ( @ARGV, $confab, $sqm, $acpype, $pdb2gro) {
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

   # if at any point the sqm calculation converges,
   # we exit from process()

   # if we did not find the conformation for which the
   # calculation would converge, try the larger number of
   # conformers

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
    print "\n\nrunning sqm for conformer $out_ctr\n";

    ( -e "sqm.$out_ctr.out") &&  `rm sqm.$out_ctr.out`;

    open (SQM, ">sqm.in") || die "Cno sqm.in: $!\n";
    print SQM @new_sqm;
    close SQM;

    $cmd = "$sqm -i sqm.in -o sqm.$out_ctr.out";
    system $cmd;

    print `tail -n7 sqm.$out_ctr.out`;

    $ret = "" || `grep \'Calculation Completed\' sqm.$out_ctr.out`;


    $ret || return;

    $newpdb = "confab.$out_ctr.pdb";
    open (NEWPDB, ">$newpdb") || die "Cannot open $newpdb: $!\n";
    print NEWPDB @conformer_matching_lines;
    close NEWPDB;

    # read the charge from the header
    $cmd = "$acpype -i $newpdb ";
    system $cmd;

    # hack the gro file 

    $template_gro = $template_pdb;
    $template_gro =~ s/\.pdb/.gro/;
    $cmd = "$pdb2gro < $template_pdb > $template_gro";
    system $cmd;

    my @my_gro = split "\n", `cat  $template_gro`;
    $line_ctr = 0;

    $acpype_dir = "confab.$out_ctr.acpype";
    chdir $acpype_dir;
    
    $filename = "$template_gro"; # new gro (note that we are in different dir now
    open (NEW_GRO, ">$filename" ) 
	|| die "Cno $filename: $!.\n";

    $filename = "confab.$out_ctr\_GMX.gro"; # gro produced my 
    open (GRO, "<$filename" ) 
	|| die "Cno $filename: $!.\n";

    foreach $acpype_gro_line ( <GRO> ) {
	chomp  $acpype_gro_line;
	if ( $line_ctr >= 2 && defined $my_gro[$line_ctr-2]) {
	    $my_line = $acpype_gro_line;
	    substr ($my_line, 20 ) = substr ($my_gro[$line_ctr-2], 20 );
	} else {
	    $my_line = $acpype_gro_line;
	}
	print NEW_GRO "$my_line\n";
	$line_ctr ++;
    }


    close GRO;
    close NEW_GRO;
    $template_itp = $template_gro;
    $template_itp =~ s/\.gro/.itp/;

    # hack the itp file - the molecule name there now is confab.$out_ctr
    `sed \'s/confab.$out_ctr/atp/g\' confab.$out_ctr\_GMX.itp > $template_itp`;
    exit;

}
