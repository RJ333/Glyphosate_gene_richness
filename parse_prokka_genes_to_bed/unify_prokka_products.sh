# get all product names from annotation in lower case and combine them with the original prokka file:
cd /data/jwerner/glyphosate/IMP/
for sample in $(find . -mindepth 1 -maxdepth 1 -type d | sed "s|^\./||")
do
  echo $sample
  for prokka_annot in /data/jwerner/glyphosate/IMP/$sample/output_IMP/Analysis/annotation/annotation.filt.gff
  do
    #echo $prokka_annot
    awk -v SAMPLE=$sample\
	  '\
	  BEGIN { FS = "\t|;"; OFS = "\t"}\
	  {match($0,/product=[^;]*/); product_value=substr($0,RSTART,RLENGTH); {print tolower(product_value), SAMPLE}}' $prokka_annot | tr '/()[:blank:]-' '@@@@@@@'  > /data/Rene/unified_${sample}.tsv
	  # pasting the modified products column to the original prokka file
	  paste -d '\t' $prokka_annot /data/Rene/unified_${sample}.tsv > /data/Rene/${sample}_prokka_modified.tsv
	  cat /data/Rene/unified_${sample}.tsv >> /data/Rene/unified_all_samples.tsv
	  rm /data/Rene/unified_${sample}.tsv
  done
done

# use the list of products from all samples as list of queries for samtools. remove duplicates and set minimum threshold
cd /data/Rene/
sort unified_all_samples.tsv| rev | sed 's/@/\t/' | rev | cut -f 1 |sort| uniq -c| awk '$1 > 1 {print $0}' | awk '{split($2,a,"="); print a[2]}' > unique_products_greater_1_all_samples_new.tsv

