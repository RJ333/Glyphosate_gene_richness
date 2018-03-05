#!/bin/sh

###filter out genes and contigs which do not match the quality requirements
#not tested yet
PATH_TO_DATA=/drives/d/1_Fachliches/gene_richness_selected_genes
min_gene_length="400"
min_contig_length="1000"

cd $PATH_TO_DATA
echo "gene length threshold is set to $min_gene_length"
echo -e "contig length threshold is set to $min_contig_length\n"
for length_file in */*_contig_and_gene_length.txt
do 
  echo "processing $length_file"
  awk -v a="$min_gene_length" -v b="$min_contig_length" '$2 >= a && $4 >= b {print $0}' $length_file > ${length_file}_trimmed.cov
  echo -e "file ${length_file}_trimmed.cov has been generated\n"
done
