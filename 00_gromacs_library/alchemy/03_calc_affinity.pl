#! /usr/bin/perl -w


@tables = ("with_protein/charge.table", "with_protein/vdw.table",
    "no_protein/charge.table", "no_protein/vdw.table");

foreach ( @tables) {
    (-e $_) || die "$_ not found\n";
}

foreach $table ( @tables) {

    $filename = $table;
    open (IF, "<$filename" ) 
	|| die "Cno $filename: $!.\n";
    
    $ctr = 0;
    while ( <IF> ) {
	next if (/^#/ ) ;
	next if (!/\S/);
	chomp;
	($from{$table}[$ctr], $to{$table}[$ctr], $dg{$table}[$ctr]) = split;
	 $ctr++;
    }
    $no_lines{$table} = $ctr;
}

close IF;

foreach $t ( 1 .. $#tables) {
    ($no_lines{$tables[$t]} == $no_lines{$tables[$t-1]})
	|| die "no lines in $tables[$t] not the same as in $tables[$t-1]\n";
}

$sum = 0;
foreach $t ( 0 .. $#tables){
    $sum[$t] = 0;
    foreach $ctr ( 0 .. $no_lines{$tables[0]}-1) {
    
	$table = $tables[$t];
	if ( $t) {
	    ($from{$table}[$ctr] == $from{$tables[$t-1]}[$ctr]) || die;
	    ($to{$table}[$ctr] == $to{$tables[$t-1]}[$ctr]) || die;
	}
	#$sum[$t] += ($to{$table}[$ctr]-$from{$table}[$ctr])*$dg{$table}[$ctr];
	if ( $tables[$t] =~ /with/ ) {
	    $sum[$t] += $dg{$table}[$ctr];
	} else {
	    $sum[$t] -= $dg{$table}[$ctr];
	}
    }
    printf "%25s  %8.3f\n", $tables[$t], $sum[$t]; 
}

$total = 0;
foreach $t ( 0 .. $#tables){
    $total += $sum[$t];
}

printf "  total %8.2f    exp(-sum) %10.4f \n", 
    $total,  exp (-$total);
