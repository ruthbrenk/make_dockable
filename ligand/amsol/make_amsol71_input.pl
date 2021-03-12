#!/usr/bin/perl

# Niu Huang, Lab of Comp. Chem, UMAB, Dec, 2002
#
# Modified by Qiuyu 2020
# $[ = 1;                 # set array base to 1
$/ = "\n";    # set input record separator
$, = ' ';     # set output field separator
$\ = "\n";    # set output record separator

if ( $#ARGV < 1 ) {
    die "missing arguments\nusage: amsol.pl input.mol2 MolName\n";
}

$inmol2file = "$ARGV[0]";
$name = "$ARGV[1]";

open( MOL2,     "<$inmol2file" ) || die "can't open $inmol2file\n";  #goto label
open( AMSOLIN1, ">temp.in-wat" ) || die "can't create *.in-wat\n";
open( AMSOLIN2, ">temp.in-hex" ) || die "can't create *.in-hex\n";

# read in a mol2 record, store all data

$_ = <MOL2>;

if (/$^\@\<TRIPOS\>MOLECULE/) {
    $_  = <MOL2>;
    chomp;
    $id = $_ ;
    $_ = <MOL2>;
    @Fld = split( ' ', $_, 999 );
    $natom = $Fld[0];
    print "Natom:",$natom;
    $nbond = $Fld[1];
    print "Nbond",$nbond;
    $heavy = $natom;    # initialize heavy atom counter
    $h     = 0;         # initialize hydrogen counter
    $chrg  = 0;         # initialize charge counter

    while (1) {
        $_ = <MOL2>;
        last if (/$^\@\<TRIPOS\>ATOM/);
    }

    unless ( /$^\@\<TRIPOS\>ATOM/ ) {
        warn "Invalid mol2 format. @<TRIPOS>ATOM not found.\n";
    }

    # determine atom type, also count Hydrogens
    for ( $i = 0 ; $i < $natom ; $i++ ) {
        $_ = <MOL2> ;
        $_ =~ s/^\s+//;
        $_ =~ s/\s+$//;
        @Fld = split /\s+/ ;
        print @Fld;

        if ( $Fld[5] =~ m/Br/ ) {
            $atom = '35';
        }
        elsif ( $Fld[5] =~ m/O/ ) {
            $atom = '8';
        }
        elsif ( $Fld[5] =~ m/Cl/ ) {
            $atom = '17';
        }
        elsif ( $Fld[5] =~ m/N/ ) {
            $atom = '7';
        }
        elsif ( $Fld[5] =~ m/H/ ) {
            $atom = '1';
            $h++;
            $heavy--;
        }
        elsif ( $Fld[5] =~ m/S/ ) {
            $atom = '16';
        }
        elsif ( $Fld[5] =~ m/F/ ) {
            $atom = '9';
        }
        elsif ( $Fld[5] =~ m/C/ ) {
            $atom = '6';
        }
        elsif ( $Fld[5] =~ m/I/ ) {
            $atom = '53';
        }
        elsif ( $Fld[5] =~ m/P/ ) {
            $atom = '15';
        }
        else {
            $atom = '99';    # delete troublesome molecules;
            warn "unkown element type: $Fld[5] in $id\n";
            close(MOL2);
            close(AMSOLIN1);
            close(AMSOLIN2);
        }
        $chrg += $Fld[5];
        $Type[$i] = $atom;
        $X[$i]    = $Fld[2];
        $Y[$i]    = $Fld[3];
        $Z[$i]    = $Fld[4];
    }


    $_ = <MOL2>;
    unless ( /^\@\<TRIPOS\>BOND/ ) {
        warn "Invalid bond line in $id";
    }
    for ( $i = 0 ; $i < $nbond ; $i++ ) {
        $_ = <MOL2>;
    }
    if ( $chrg < 0 ) {
        $netchg = int $chrg - 0.5;    # adjust formal charge
    }
    else {
        $netchg = int $chrg + 0.5;    # adjust formal charge
    }

    # write amsol input with Cartesian coordinates
    $keyword1 =
      "SOLVNT=WATER CHARGE=$netchg AM1 SM5.42R 1SCF GEO-OK CART TLIMIT=100"
      ;                               # increase time limit to 100 s.
    $keyword2 =
"SOLVNT=GENORG IOFR=1.4345 ALPHA=0.00 BETA=0.00 GAMMA=38.93\n& DIELEC=2.06 FACARB=0.00 FEHALO=0.00 CHARGE=$netchg AM1 SM5.42R 1SCF\n& GEO-OK CART TLIMIT=100";

    if ( 1 ) {

#      print AMSOLIN1 " ";      # one blank line precede each record, not working for amsol6.8
#      print AMSOLIN2 " ";
        print AMSOLIN1 $keyword1;
        print AMSOLIN2 $keyword2;
        print AMSOLIN1 $name, $natom,
          "\n";    #blank line to separate keyword and Z matrix, required
        print AMSOLIN2 $name, $natom, "\n";

        for ( $i = 0 ; $i < $natom ; $i++ ) {
            printf AMSOLIN1 "%1s %10.4f %1s %10.4f %1s %10.4f %1s\n",
              $Type[$i], $X[$i], "1", $Y[$i], "1", $Z[$i], "1";
            printf AMSOLIN2 "%1s %10.4f %1s %10.4f %1s %10.4f %1s\n",
              $Type[$i], $X[$i], "1", $Y[$i], "1", $Z[$i], "1";
        }

        print AMSOLIN1 " ";    # one blank line trailing each record
        print AMSOLIN2 " ";

    }
}

close(MOL2)     || die "can't close $inmol2file.mol2\n";
close(AMSOLIN1) || die "can't close temp.in-wat\n";
close(AMSOLIN2) || die "can't close temp.in-hex\n";

