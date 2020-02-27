#!/bin/bash
counter=0
for F in `ls -1 *.gz` ; do
  counter=$(($counter + 1))
  echo $F
  in_file="A$counter.fastq"
  echo $in_file
  gunzip -c $F > $in_file
  perl ./01bcfreq_parse.pl $in_file 01_Hsp90_full.input 10
  rm $in_file
done
