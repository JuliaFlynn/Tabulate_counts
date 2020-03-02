#!/usr/bin/perl
#
# read in sorted sequence file (sequences only)
# count number of each unique sequence
# output each unique sequence and its number of occurences
#
if ($#ARGV !=1) {
  print "usage: sort2unique.pl infile outfile\n" ;
  exit ;
}

$infile = $ARGV[0] ;
$outfile = $ARGV[1] ;

open(INF, $infile) ;
open(OUTF, ">$outfile") ;

$line2 = "ZZZ" ;
$count = 1 ;
while ($line = <INF>) {
  chomp($line) ;
  if ($line eq $line2) {
    $count++;
  } else {
    print OUTF "$line2\,$count\n" ;
    $line2 = $line ;
    $count = 1 ;
  }
}

print OUTF "$line\,$count\n" ;

close(INF);
close(OUTF);
