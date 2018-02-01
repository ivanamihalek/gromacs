#! /usr/bin/perl -w
sub fix_impl_params (@);
 
(@ARGV > 1 )  ||
    die "Usage: $0 <ligand name>  <new name>.\n";

($ligand, $new_name) = @ARGV;
$gromacs_path = "/usr/local/bin";
$grompp      = "$gromacs_path/grompp";
foreach  ( "../00_input", "../05_mpt_eq",  "../01_topology", $grompp) {
    ( -e $_ ) || die "\n$_ not found.\n\n";

}

# find top file
$ret = "" || `ls ../01_topology/*top`;
$ret || die "no top file found in  ../01_topology\n";
@topfile = split "\n", $ret;
(@topfile > 1 ) && die "more than one top file in  ../01_topology\n";
print "top file: $topfile[0]\n";

# find gro file
$ret = "" || `ls ../05_mpt_eq/*.gro`;
$ret || die "no gro file found in  ../05_mpt_eq/\n";
@grofile = split "\n", $ret;
(@grofile > 1 ) && die "more than one gro file in ../05_mpt_eq/ \n";
print "gro file: $grofile[0]\n";

# find mdp file
(-e "../00_input/md.mdp") || die "../00_input/md.mdp not found\n";
print "mdp file: ../00_input/md.mdp\n";

# find itp file
$itpfile = "../00_input/$ligand.itp";
(-e "$itpfile") || die "$itpfile not found\n";
print "itp file: $itpfile\n";


# add the genborn params to the *itp file
$ret = "" || `grep implicit_genborn_params $itpfile`;
$ret || fix_impl_params ("$itpfile", "$new_name.itp");
$pwd = `pwd`; chomp $pwd;

# make new top  file to include the genborn params
$new_top = "$new_name.top";
open (OF, ">$new_top") || die "Cno $new_name.top.\n";
@lines = split "\n", `cat  $topfile[0]`;
foreach $line (@lines) {
    if ( $line =~ /include/ && $line =~ /$ligand/) {
	print OF "#include \"$pwd/$new_name.itp\"\n";
    } else {
	print OF $line."\n",
    }
}
close OF;

# make the corresponding new tpr  file
$program = "$grompp";
$log     = "grompp.log";
$command = "$program  -f ../00_input/md.mdp -c  $grofile[0] ".
    " -p  $new_top  -o $new_name.tpr > $log";
(system $command) && die "Error:\n$command\n"; 



############################################################################## 
##############################################################################
sub fix_impl_params (@){


    #($forcefield  eq "amber99sb") 
    #	|| die "GB params assumes  amber99sb forcefield ...\n";
 
    my $itpfile     = $_[0];
    my $new_itpfile = $_[1];
    my %params = ();

    # sp2, all-atom, aromatic
    $params{"C"}  = "         0.172    1      1.554    0.1875    0.72 ; C";
    $params{"C*"} = "         0.172    0.012  1.554    0.1875    0.72 ; C";
    $params{"CA"} = "         0.18     1      1.037    0.1875    0.72 ; C";
    $params{"CB"} = "         0.172    0.012  1.554    0.1875    0.72 ; C";
    $params{"CC"} = "         0.172    1      1.554    0.1875    0.72 ; C";

    $params{"CN"} = "           0.172    0.012  1.554    0.1875    0.72 ; C";
    $params{"CR"} = "           0.18     1      1.073    0.1875    0.72 ; C";
    $params{"CV"} = "           0.18     1      1.073    0.1875    0.72 ; C";
    $params{"CW"} = "           0.18     1      1.073    0.1875    0.72 ; C";

    # sp3 all-atom 
    $params{"CT"} = "           0.18     1      1.276    0.190     0.72 ; C";

 
    $params{"H"} =  "           0.1      1      1        0.115     0.85 ; H";
    $params{"HC"} = "           0.1      1      1        0.125     0.85 ; H";
    $params{"H1"} = "           0.1      1      1        0.125     0.85 ; H";
    $params{"HA"} = "           0.1      1      1        0.125     0.85 ; H";
    $params{"H4"} = "           0.1      1      1        0.115     0.85 ; H";
    $params{"H5"} = "           0.1      1      1        0.125     0.85 ; H";
    $params{"HO"} = "           0.1      1      1        0.105     0.85 ; H";
    $params{"HS"} = "           0.1      1      1        0.125     0.85 ; H";
    $params{"HP"} = "           0.1      1      1        0.125     0.85 ; H";

    $params{"N"} =  "           0.155    1      1.028    0.17063   0.79 ; N";
    $params{"NA"} = "           0.155    1      1.028    0.17063   0.79 ; N";
    $params{"NB"} = "           0.155    1      1.215    0.17063   0.79 ; N";
    $params{"N2"} = "           0.16     1      1.215    0.17063   0.79 ; N";
    $params{"N3"} = "           0.16     1      1.215    0.1625    0.79 ; N";

    $params{"O"}  = "           0.15     1      0.926    0.148     0.85 ; O";
    $params{"OH"} = "           0.152    1      1.080    0.1535    0.85 ; O";
    $params{"O2"} = "           0.17     1      0.922    0.148     0.85 ; O";
    $params{"S"}  = "           0.18     1      1.121    0.1775    0.96 ; S";
    $params{"SH"} = "           0.18     1      1.121    0.1775    0.96 ; S";
    $params{"BR"} = "           0.1      1      1        0.125     0.85 ; BR";
    $params{"F"}  = "           0.1      1      1        0.156     0.85 ; F";
    $params{"CL"} = "           0.1      1      1        0.70      0.85 ; CL";
    $params{"I"}  = "           0.1      1      1        0.206     0.85 ; I";
    $params{"P5"} = "           0.1      1      1        0.190     0.85 ; P5";
         
    # take as identical (Ivana);
    # these are "united" how did they end up in the same list as all-atom?
    # sp3
    $params{"C2"} =  $params{"CT"}; 
    $params{"C3"} =  $params{"CT"}; 
    $params{"CH"} =  $params{"CT"}; 
    $params{"CS"} =  $params{"CT"}; 
    $params{"CD"} =  $params{"CA"}; # should united atome appear here at all ...?
    $params{"CP"} =  $params{"CA"};
    $params{"CX"} =  $params{"CT"}; # CX - tip of that funny three-membered ring
    $params{"CE"} =  $params{"C*"};
    $params{"CG"} =  $params{"C*"};
    $params{"H2"} =  $params{"H"};

    $params{"NT"} =  $params{"N3"}; # sp3 nitrogen with 4 substituents
    $params{"NH"} =  $params{"N2"}; # sp2 nitrogen in base NH2 group or arginine NH2
    $params{"N*"} =  $params{"N2"}; # sp2 nitrogen in base NH2 group or arginine NH2
    $params{"ND"} =  $params{"N"};  # sp2 nitrogen in amide
    $params{"HN"} =  $params{"H"};  # amide or imino hydrogen


    $params{"OS"} =  $params{"OH"}; # sther or esther O params  replaced by alcohol
    $params{"SS"} =  $params{"S"};    
    $params{"N1"} =  $params{"N2"}; #triple bond in CN?
    $params{"NC"} =  $params{"NB"};
    $params{"NO"} =  $params{"NB"}; # nitrobenzyl?
    $params{"SY"} =  $params{"S"}; #sulfur in 5 ring member
    $params{"SS"} =  $params{"S"}; # sulfur dioxide
    $params{"C1"} =  $params{"C"}; #sp carbon

    my @lines = split "\n", `cat $itpfile`;
    my $new_itp = "";
    my $new_field  =   "[ implicit_genborn_params ]\n".
	";name    sar      st     pi       gbr       hct\n";
    my $reading = 0;
    my $name;
   
    foreach my $line (@lines) {

	if ( $line =~ /atomtype/ ) {
	    $reading = 1;

	} elsif ( $reading && $line =~ /\[/ ) {
	    
	    # add the hacked gb field;
	    $new_itp .= $new_field."\n";
	    $reading = 0;

	} elsif  ($reading  && $line =~ /\S/ && $line !~ /name/) {  # not empty or a header a header line
		($name) = split " ", $line;
		$new_field .= " $name ";
		if ( defined $params{$name} ) {
		    $new_field .= $params{$name};

		} elsif ( defined $params{uc $name} ) {
		    $new_field .= $params{uc $name};

		#} elsif (  defined $params{uc substr $name, 0, 1}  ) {
		#    $new_field .= defined $params{uc substr $name, 0, 1} ;

		} else {
		    die "GB params for atom type $name not found.\n";
		}
		$new_field .= "\n";

		
	}

	$new_itp .= $line."\n";

	

    }

    open (NEW_ITP, ">$new_itpfile") || die "Cno $new_itpfile: $!\n";
    print NEW_ITP $new_itp;
    close NEW_ITP;

 
    return;
}

