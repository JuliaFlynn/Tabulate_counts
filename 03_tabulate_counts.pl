#!/usr/bin/perl
#
# read in uniq files of barcodes for different timepoints
# read in barcode-ORF assembly files
# write counts file organized by position, amino acid, codon

if ($#ARGV != 2) {
  print "usage: script.pl cleaned_bc_ORF_assembly_file 03_description_file output_file\n";
  exit;
}

$bc_ORF_in = $ARGV[0];
$description_file = $ARGV[1];
$out_file = $ARGV[2];

print "bc-ORF assembly file: $bc_ORF_in\n" ;
print "tabulate description file: $description_file\n" ;

open(DF, $description_file) ;
$line = <DF> ;
chomp($line) ;
@spline = split (/,/, $line) ;
if ($spline[0] eq "timepoints") {
  $ntimepoints = $spline[1] ;
} else {
  print "Error: expected timepoints statement at line 1 of description file\n";
  exit ;
}
print "Number of timepoints to tabulate: $ntimepoints\n" ;

$line = <DF> ;

for ($i=0; $i<$ntimepoints; $i++) {
  $line = <DF> ;
  chomp($line) ;
  @spline = split (/,/, $line) ;
  $uniq_filename[$i]=$spline[0] ;
  $hour[$i] = $spline[1] ;
}

close(DF) ;

open(BCF, $bc_ORF_in) ;

open(OUTF, ">$out_file") ;

print OUTF "Position, aa mutation, codon mutation, barcode for hours:" ;
for ($i=0; $i<$ntimepoints; $i++) {
  print OUTF ",$hour[$i]" ;
}
print OUTF "\n" ;

$line = <BCF> ;

while ($line = <BCF>) {
  chomp($line) ;
  @spline = split (/,/, $line) ;
  $nspline = @spline ;
  $pos = $spline[0];
  $codon_mut = $spline[1];
  $aa_mut = $spline[2];
  if ($nspline > 3) {
    for ($i=3; $i<$nspline; $i++) {
      $bc_cur = $spline[$i] ;

      print OUTF "$pos,$aa_mut,$codon_mut,$bc_cur" ;

      for ($j=0; $j<$ntimepoints; $j++) {
        $inseq = $uniq_filename[$j] ;
        $searchresult = `grep $bc_cur $inseq` ;
        chomp($searchresult) ;
        @splresult = split(/,/, $searchresult) ;
        $searchcount = $splresult[1] ;
        print OUTF ",$searchcount" ;
      }
      print OUTF "\n" ;
    }
  }
}

close(BCF) ;
close(OUTF) ;
