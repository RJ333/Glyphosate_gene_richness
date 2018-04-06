#!/bin/bash
#get all product names from annotation:
cd /data/jwerner/glyphosate/IMP/
for sample in $(find . -mindepth 1 -maxdepth 1 -type d | sed "s|^\./||")
do
  #echo $sample
  for prokka_annot in /data/jwerner/glyphosate/IMP/$sample/output_IMP/Analysis/annotation/annotation.filt.gff
  do
    #echo $prokka_annot
    awk -v SAMPLE=$sample\
	  '\
	  BEGIN { FS = "\t|;"; OFS = "\t"}\
	  {match($0,/product=[^;]*/); product_value=substr($0,RSTART,RLENGTH); gsub(" ","@",product_value); {print product_value, SAMPLE}}' $prokka_annot |\
	  awk '{gsub("/","@")} 1' >> /data/Rene/unique_products_per_sample.tsv
  done
done
cd /data/Rene/
sort unique_products_per_sample.tsv| tr "/" "@" | uniq | cut -f 1 | sort| uniq -c| awk '$1 > 6 {print $0}' | awk '{split($2,a,"="); print a[2]}' > unique_products_greater_6_all_samples.tsv
