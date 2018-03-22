#!/bin/sh

<<<<<<< HEAD:single_assembly_gene_richness_and_reads_mapping/03_clean_up_appended_gene_lists.sh
# BEGIN adds table header and sets output field sep to tab
# split() adresses specific column and splits on "_", storing the pieces in array a
# from this leftover a[3], another split is performed to remove the ".tsv", stored in array b
=======
# puts all single files into one file

for nfiles in named_*.tsv
do	
  cat $nfiles >> appended_genes.tsv
  echo "appending $nfiles done"
done

echo "cleaning up table"

# BEGIN adds table header and sets Output field sep to tab
# split adresses specific column and splits on "_", storing the pieces in array a
# from this leftover, another split is performed to remove the ".tsv", stored in array b
>>>>>>> a6bfaabca3ef6b582bbbf357bdd366aaa9e3ef70:single_assembly_gene_richness_and_reads_mapping/02_append_and_clean_samtools_view_output.sh
# the respective fields are printed

awk 'BEGIN { OFS = "\t" ; print "contig\tgene\tsample\treads_per_gene"}
  {
  split ($1, a, "_") 
  split (a[3], b, "\\.") 
  print $2, a[2], b[1], $3
  }' appended_genes*tsv > contig_gene_sample_reads.tsv
  
rm appended_genes*tsv

echo "cleaning is done, output written to contig_gene_sample_reads.tsv"
