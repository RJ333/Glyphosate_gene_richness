#!/bin/sh

# this script combines the output of the script "xy" 
# (is there are script that contains the "samtools depth" step?) 
# to one list for all genes in a folder
# and sorts it by coverage

PATH_TO_DATA=/data/jwerner/glyphosate/IMP/ 

cd $PATH_TO_DATA
# find only directories
for sample in $(find . -mindepth 1 -maxdepth 1 -type d) 
do 
  name=$(echo $sample | cut -c 3-)
  cd $sample/output_IMP/Analysis
  echo $sample
  echo $name
  # add all 10 samples to one list and add a column with the sample name  
  awk -v name="$name" 'BEGIN{FS = OFS = "\t"}{print name, $0}' mg.assembly.length.txt >> /data/Rene/contig_length.tsv 
  cd $PATH_TO_DATA 
done