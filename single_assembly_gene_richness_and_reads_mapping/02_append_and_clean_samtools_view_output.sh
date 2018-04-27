#!/bin/sh

# BEGIN adds table header and sets output field sep to tab
# split() adresses specific column and splits on "_", storing the pieces in array a
# from this leftover a[3], another split is performed to remove the ".tsv", stored in array b
# puts all single files into one file
cd /data/Rene/glyph/named_reads_per_gene

for nfiles in named_*.tsv
do	
  cat $nfiles >> appended_genes.tsv
  echo "appending $nfiles done"
done

echo "cleaning up table"

# BEGIN adds table header and sets Output field sep to tab
# split adresses specific column and splits on "_", storing the pieces in array a
# from this leftover, another split is performed to remove the ".tsv", stored in array b
# the respective fields are printed

awk 'BEGIN { OFS = "\t" ; print "contig\tgene\tsample\treads_per_gene"}
  {
  split ($1, a, "_|__") 
  split (a[3], b, "\\.") 
  print $2, a[2], b[1], $3
  }' appended_genes*tsv > contig_gene_sample_reads.tsv
  
rm appended_genes*tsv

echo "cleaning is done, output written to contig_gene_sample_reads.tsv"
