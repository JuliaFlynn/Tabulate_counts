#!/usr/bin/perl
#
# read in single 01seq file containing bc,index for samples matching constant
#   region (e.g. parsed based on parent plasmid/construct)
# read in file containing information on how to parse indexes
#   and output a separate file for each index with bc only

sub basecomp {
  my($inchar, $outchar) ;
  $inchar = uc(substr($_[0],0,1)) ;
  $outchar = "N" ;
  if($inchar eq "A") {
    $outchar = "T" ;
  }
  if($inchar eq "C") {
    $outchar = "G" ;
  }
  if($inchar eq "G") {
    $outchar = "C" ;
  }
  if($inchar eq "T") {
    $outchar = "A" ;
  }
  $outchar ;
}


sub revcomp {
  my($inseq, $outseq, $seqlen, $p, $q, $inchar, $outchar);
  $inseq = $_[0] ;
  $seqlen = length($inseq); 
  $outseq = "" ;
  for ($p=$seqlen-1;$p>=0;$p=$p-1) {
    $inchar = substr($inseq,$p,1) ;
    $outchar = "N" ;
    $outchar = &basecomp($inchar) ;
    $outseq = $outseq . $outchar ;
  }
  $outseq ;
}


if ($#ARGV != 1) {
  print "usage: script.pl 01_file_input parse2_description_file\n";
  exit;
}

$seq_in = $ARGV[0];
$description_file = $ARGV[1];

print "seq file: $seq_in\n" ;
print "parse2 description file: $description_file\n" ;

open(DF, $description_file) ;
$line = <DF> ;
chomp($line) ;
@spline = split (/,/, $line) ;
$nparse = $spline[1] ;
$index_len = $spline[2] ;
print "Number of files to output: $nparse\n" ;
$line = <DF> ;
for ($i=0; $i<$nparse; $i++) {
  $line = <DF> ;
  chomp($line) ;
  @spline = split (/,/, $line) ;
  $prefix[$i]=$spline[0] ;
  $index_seq[$i] = &revcomp($spline[1]) ;
# reverse complement because read from opposite strand as end2 primer
}

close(DF) ;

open(INF, $seq_in) ;
open(DUMPF, ">", "leftover.seq") ;

my @file_h ;
for my $file (0..($nparse-1)) {
  open($file_h[$file], ">", "$prefix[$file].bc") ;
}

while ($line = <INF>) {
  chomp($line) ;
  @spline = split (/,/, $line) ;
  $index_cur = $spline[1] ;
  $bc_cur = $spline[0] ;

# update below still ongoing

  @list_int = (0..($nparse-1)) ;
  $parse_check = 0 ;
  for $i (@list_int) {
    @list_pos = (0..($index_len-1)) ;
    $differ = 0 ;
    for $j (@list_pos) {
      $let1 = substr($index_cur,$j,1) ;
      $let2 = substr($index_seq[$i],$j,1) ;
      if ($let1 ne $let2) {
        $differ = 1 ;
        last ;
      }
    }
    if ($differ == 0) {
# send bc to output
      $parse_check = 1 ;
      print { $file_h[$i] } $bc_cur, "\n";
    }
  }
  if ($parse_check == 0) {
# unparsed bc - send line to DUMPF
    print DUMPF "$line\n" ;
  }
}

close(INF) ;
close(DUMPF) ;
for ($i=0; $i<$nparse; $i++) {
  {$file_h[$i]} close() ;
}
