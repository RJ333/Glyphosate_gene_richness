#!/bin/sh

# this script adds the annotated gene region the coverage was calculated for. 
# It became complicated as it was supposed to distuingish between two genes placed on the same contig, only differing by their annotated region.
# it now continously counts all unique occurrences of the gene based on the data of the bed file
#
#with a lot of help from
#https://stackoverflow.com/questions/49050245/how-to-use-awk-to-add-specific-values-to-a-column-based-on-numeric-ranges/49054710?noredirect=1#comment85120455_49054710


PATH_TO_DATA=/drives/d/1_Fachliches/gene_richness_selected_genes

cd $PATH_TO_DATA
for gene in $(find . -mindepth 1 -maxdepth 1 -type d) 
do 
  cd $gene
  echo "$gene is being processed"
  awk 'NR==FNR{start[NR]=$2; end[NR]=$3; key[$1,$2]=$4 sprintf("_%03d",NR); next}
           {for(i in start)
              {s=start[i];
               if(s<=$2 && $2<=end[i] && ($1,s) in key) print $0,key[$1,s] | "sort -nk1 "}}' intersect.bed ${gene}_combined.txt_filtered_5.txt > ${gene}_counts.cov
  echo "file ${gene}_counts.cov generated, leaving folder $gene"
  cd ..
done
