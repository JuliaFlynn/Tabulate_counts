#!/bin/bash
prefix="Hsp90"
for i in {1..12}; do
  echo $i
  in_file="$prefix""_ID_$i.bc"
  out_file="$prefix""_ID_$i.sortbc"
  echo $in_file
  echo $out_file
  sort $in_file > $out_file
  in_file="$prefix""_ID_$i.sortbc"
  out_file="$prefix""_ID_$i.uniq"
  echo $in_file
  echo $out_file
  perl sort2uniq.pl $in_file $out_file
done
