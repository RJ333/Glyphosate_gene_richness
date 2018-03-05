#!/bin/sh

#this script removes bases with coverage below threshold from previously created coverage files

PATH_TO_DATA=/drives/d/1_Fachliches/gene_richness_selected_genes
threshold=5  # coverage must be equal or higher to be kept

cd $PATH_TO_DATA
echo "threshold is set to $threshold";
for coverage_file in */*_combined.txt
do
  echo $coverage_file;
  awk -v a="$threshold" '$3 >= a {print $0}' $coverage_file > ${coverage_file}_filtered_${threshold}.txt
done
