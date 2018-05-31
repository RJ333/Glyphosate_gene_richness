#!/bin/sh

OUTPUT_DIR=/data/Rene/glyph/prokka 
mkdir -p $OUTPUT_DIR/named

# adds the name of the file in folder "unnamed" as first column and copies to "named" folder 
cd $OUTPUT_DIR/renamed
for files in *.tsv
do	
  awk '{a=FILENAME;}{print a"\t"$0}' $files > $OUTPUT_DIR/named/named_${files}
done

rm -r $OUTPUT_DIR/renamed

cd $OUTPUT_DIR/named

for nfiles in named_*.tsv
do	
  cat $nfiles >> appended_genes.tsv
done

echo "all information appended into appended_genes.tsv"
echo "cleaning up appended_genes.tsv"

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