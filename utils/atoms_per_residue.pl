#!/usr/bin/perl -w
use strict;
use warnings;


my $ff_file = "/usr/local/gromacs/share/gromacs/top/amber99.ff/aminoacids.rtp";

open(INFILE, "<$ff_file") || die "Cno $ff_file: $!\n";

my ($reading_res, $reading_atoms) = (0, 0);
my @residues = ('ALA', 'VAL', 'ILE', 'LEU', 'MET', 'GLY',
                'GLU', 'GLN', 'ASP', 'ASN',
                'LYS', 'ARG', 'HIE', # histidine with H on epsilon carbon
                'SER', 'THR', 'CYS',
                'TRP', 'PHE', 'TYR', 'PRO');

my $residues = "_".(join("_",@residues))."_";

my %atoms = ();

my $res;
while (<INFILE>) {
    if (/\[\s*(\S*)\s*\]/) {
        my $kwd = $1;
        if ($residues =~ "_$kwd\_") {
            $reading_res = 1;
            $res = $kwd;
        } elsif ($kwd eq "atoms") {
            $reading_atoms = 1;
        } else {
            ($reading_res, $reading_atoms) = (0, 0);
        }
    } elsif ($reading_res && $reading_atoms) {
        defined $atoms{$res} ||  ($atoms{$res} = ());
        my @field = split;
        my $atomname = $field[0];
        (substr $atomname, 0, 1)eq'H' || push @{$atoms{$res}}, $atomname;
    }
}


for my $res(keys %atoms) {
    #print "$res  @{$atoms{$res}}\n";
    my $atomct = scalar @{$atoms{$res}};
    print "\"$res\":$atomct, ";
}
print "\n";