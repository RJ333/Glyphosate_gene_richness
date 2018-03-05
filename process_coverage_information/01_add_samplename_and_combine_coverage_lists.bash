#!/bin/sh

# this script combines the output of the script "xy" 
# (is there are script that contains the "samtools depth" step?) 
# to one list for all genes in a folder
# and sorts it by coverage

PATH_TO_DATA=/drives/d/1_Fachliches/gene_richness_selected_genes 

cd $PATH_TO_DATA
# find only directories
for gene in $(find . -mindepth 1 -maxdepth 1 -type d) 
do 
  cd $gene
  # add all 10 samples to one list and add a column with the sample name
  for sample in *.sorted.removeduplicates.bam.cov.tmp
  do
    echo "$gene is being processed"  
    awk 'BEGIN{OFS="\t"}{print $0, FILENAME}' $sample >> ${gene}_tmp.txt 
  done
  # split the long sample name in an array "a" at ".", take only the first part of that array, sort ascendingly by column3 "coverage"
  awk	'{split ($4,a,".");print $1,$2,$3,a[1]|"sort -nk3 "}' ${gene}_tmp.txt > ${gene}_combined.txt 
  # clean up intermediate files
  rm ${gene}_tmp.txt 
  cd .. 
done