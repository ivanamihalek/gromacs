#! /usr/bin/perl -w


# from spc216.gro
#   1SOL     OW    1    .230    .628    .113
#    1SOL    HW1    2    .137    .626    .150
#    1SOL    HW2    3    .231    .589    .021

@template_O  = ( .230,    .628,    .113);
@template_H1 = ( .137,    .626,    .150);
@template_H2 = ( .231,    .589,    .021);

$new_serial = 0;

while ( <> ) {

    next if ( $_ !~ /^ATOM/ &&  $_ !~/^HETATM/);
    $serial = substr $_, 6, 5;  $serial =~ s/\s//g;
    $name = substr $_,  12, 4 ;  $name =~ s/\s//g;
    $res_name = substr $_,  17, 3; $res_name=~ s/\s//g;
    $res_seq  = substr $_, 22, 4;  $res_seq=~ s/\s//g;
    $x = substr $_,30, 8;  $x=~ s/\s//g;
    $y = substr $_,38, 8;  $y=~ s/\s//g;
    $z = substr $_, 46, 8; $z=~ s/\s//g;

    if ( $res_name eq "HOH" || $res_name eq "SOL") {

	# add my own hydrogens
	($name =~ /^H/) && next;

	$res_name = "SOL";
	$new_serial++;
	$serial = $new_serial;

	@my_O = ( $x/10, $y/10, $z/10);
	$name = "OW";
	printf "%5d%5s  %-3s%5d%8.3f%8.3f%8.3f\n",
	$res_seq, $res_name,  $name, $serial, @my_O;

	for $i ( 0 .. 2 ) {
	    $transl[$i] = $my_O[$i] - $template_O[$i];
	    $H1[$i] =  $template_H1[$i] + $transl[$i];
	    $H2[$i] =  $template_H2[$i] + $transl[$i];
	}
	$name = "HW1";
	printf "%5d%5s  %-3s%5d%8.3f%8.3f%8.3f\n",
	$res_seq, $res_name,  $name, $serial, @H1;
	$name = "HW2";
	printf "%5d%5s  %-3s%5d%8.3f%8.3f%8.3f\n",
	$res_seq, $res_name,  $name, $serial, @H2;

	
    } else {
	printf "%5d%5s  %-3s%5d%8.3f%8.3f%8.3f\n",
	$res_seq, $res_name,  $name, $serial, $x/10, $y/10, $z/10;
    }


}
