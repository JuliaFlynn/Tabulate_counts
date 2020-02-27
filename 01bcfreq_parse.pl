#!/usr/bin/perl
#
# read in single end fastq file
# read in file containing information on how to parse samples
#   including sequence for constant region
# for sequences that match constant region and have PHRED>cutoff
#   output bc and index sequences to parsed files (one for each 
#   constant region
#   NOTE: match updated to tolerate up to 2 snp - output as separate files
#         *.snp0 *.snp1 and *.snp2
# output unparsed sequences as leftover.fastq

if ($#ARGV != 2) {
  print "usage: script.pl fastq parse_description_file PHRED_cutoff\n";
  exit;
}

$fastq_in = $ARGV[0];
$description_file = $ARGV[1];
$phred_cut = $ARGV[2];
print "phred cut: $phred_cut\n" ;

print "fastq file: $fastq_in\n" ;
print "parse description file: $description_file\n" ;

open(DF, $description_file) ;
$line = <DF> ;
chomp($line) ;
@spline = split (/,/, $line) ;
$nparse = $spline[1] ;
print "Number of files to output: $nparse * 3\n" ;
for ($i=0; $i<$nparse; $i++) {
  $line = <DF> ;
  chomp($line) ;
  @spline = split (/,/, $line) ;
  $prefix[$i]=$spline[0] ;

  $temp = $spline[2] ;
  @stemp = split (/_/, $temp) ;
  $qcheck1_start[$i] = $stemp[1] - 1 ;
  $qcheck1_len[$i] = $stemp[2] - $stemp[1] + 1 ;

  $temp = $spline[4] ;
  @stemp = split (/_/, $temp) ;
  $qcheck2_start[$i] = $stemp[1] - 1 ;
  $qcheck2_len[$i] = $stemp[2] - $stemp[1] + 1 ;
 
#  print "qcheck2 $i: $qcheck2_start[$i] - $qcheck2_len[$i]\n" ;

  $temp = $spline[3] ;
  @stemp = split (/_/, $temp) ;
  $constant_start[$i] = $stemp[1] - 1 ;
  $constant_len[$i] = $stemp[2] - $stemp[1] + 1 ;

#  print "$i $constant_start[$i] $constant_len[$i]\n" ;

  $line = <DF> ;
  chomp($line) ;
  @spline = split (/,/, $line) ;
  $constant_seq[$i] = $spline[1] ;

#  print "constant: $constant_seq[$i]\n" ;
}

$con_start = $constant_start[0] ;
$con_len = $constant_len[0] ;
$q1_start = $qcheck1_start[0] ;
$q1_len = $qcheck1_len[0] ;
$q2_start = $qcheck2_start[0] ;
$q2_len = $qcheck2_len[0] ;
#print "q2: $q2_start -- $q2_len\n" ;
$sanity_check = 1 ;
for ($i=0; $i<$nparse; $i++) {
  if ($con_start != $constant_start[$i]) {$sanity_check=0} ;
  if ($con_len != $constant_len[$i]) {$sanity_check=0} ;
  if ($q1_start != $qcheck1_start[$i]) {$sanity_check=0} ;
  if ($q1_len != $qcheck1_len[$i]) {$sanity_check=0} ;
  if ($q2_start != $qcheck2_start[$i]) {$sanity_check=0} ;
  if ($q2_len != $qcheck2_len[$i]) {$sanity_check=0} ;
}

if ($sanity_check !=1) {
  print "failed sanity check! constant locations differ!\n" ;
  exit ;
}

close(DF) ;

$ntrac = $q1_len + $q2_len ;
for ($i=0; $i<$ntrac; $i++) {
  for ($j=0; $j<41; $j++) {
    $qtrac[$i][$j] = 0 ;
  }
}

open(INF, $fastq_in) ;
open(DUMPF, ">>", "leftover.fastq") ;

my @file_h ;
for my $file (0..($nparse-1)) {
  open($file_h[$file*3], ">>", "$prefix[$file]_bc_ind_snp0.seq") ;
  open($file_h[$file*3+1], ">>", "$prefix[$file]_bc_ind_snp1.seq") ;
  open($file_h[$file*3+2], ">>", "$prefix[$file]_bc_ind_snp2.seq") ;
}

$countlines = 0;
$check = 0 ;

while ($line = <INF>) {
  chomp($line) ;
  $check++ ;
  if ($check == 1) {$line1=$line} ;
  if ($check == 2) {
    $seq = $line ;
    $line2 = $line ;
  }
  if ($check == 3) {$line3=$line} ;
  if ($check == 4) {
    $line4 = $line ;
    $qual = $line ;
    $check = 0 ;
# determine if phred scores of bc and index pass cutoff
    $q1 = substr($qual,$q1_start,$q1_len) ;
    $q2 = substr($qual,$q2_start,$q2_len) ;
    $qtest = $q1 . $q2 ;
    @qc = split (//, $qtest) ;
    $testposition=0 ;
    $quality_ok = 1 ;
    foreach (@qc) {
      $qval = ord($_);
      $qval = $qval - 33 ;
      if ($qval < $phred_cut) {$quality_ok=0} ;
      $qtrac[$testposition][$qval]++ ;
      $testposition++ ;
    }
    if ($quality_ok == 1) {
#      print "quality ok, $seq\n" ;
# check to see if constant region is a match

      $parsed_ok = 0 ;

      $testseq = substr($seq,$con_start,$con_len) ;

#      print "testseq: $testseq\n" ;

      @list_int = (0..($nparse-1)) ;
      for $i (@list_int) {

#        print "constseq: $constant_seq[$i]\n" ;

        $n_snp = 0 ;
        @list_j = (0..($con_len-1)) ;
        for $j (@list_j) {
          $let1 = substr($testseq,$j,1) ;
          $let2 = substr($constant_seq[$i],$j,1) ;
          if ($let1 ne $let2) {$n_snp++} ;
          if ($n_snp > 2) {

#            print "blasted\n" ;

            last ; 
          }
        }

        if ($n_snp < 3) {
# output bc and index sequences to appropriate file
          $outstring = substr($seq,$q1_start,$q1_len) ;
#          print "seq $q2_start $q2_len $seq\n" ;
          $iseq = substr($seq,$q2_start,$q2_len) ;
          $outstring = $outstring . "," . $iseq ;
          $temp = $i*3+$n_snp ;
#          print "i $i temp $temp\n" ;
          print { $file_h[$temp] } $outstring, "\n";
          $parsed_ok = 1 ;
        }
      } 
      if ($parsed_ok != 1) {
# constant region didn't match, dump to leftover.fastq 
        print DUMPF "$line1\n" ;
        print DUMPF "$line2\n" ;
        print DUMPF "$line3\n" ;
        print DUMPF "$line4\n" ;
      }
    } else {
# quality failed, dumb seq to leftover.fastq
      print DUMPF "$line1\n" ;
      print DUMPF "$line2\n" ;
      print DUMPF "$line3\n" ;
      print DUMPF "$line4\n" ;
    }
  }
}


close(INF) ;
close(DUMPF) ;
for ($i=0; $i<$nparse; $i++) {
  {$file_h[$i*3]} close() ;
  {$file_h[$i*3+1]} close() ;
  {$file_h[$i*3+2]} close() ;
}
open(SUMF, ">", "phred.sum") ;
for ($i=1; $i<41; $i++) {
  print SUMF "qscore_$i\," ;
  for ($j=0; $j<$ntrac; $j++) {
    print SUMF "$qtrac[$j][$i]\," ;
  }
  print SUMF "\n" ;
}
close(SUMF) ;
