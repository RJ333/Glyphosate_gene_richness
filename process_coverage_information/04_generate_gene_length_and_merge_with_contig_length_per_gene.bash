#!/bin/sh

#extract gene length from prokka annotation 
#this works only for genes which have an annotated ec number, otherwise it picks the wrong column
PATH_TO_DATA=/drives/d/1_Fachliches/gene_richness_selected_genes

cd $PATH_TO_DATA
for gene in $(find . -mindepth 1 -maxdepth 1 -type d) 
do
  echo "$gene is stored in the for loop"
  var=$(echo $gene | cut -c 3- |awk -vFS= -vOFS= '{$NF=toupper($NF)}1')
  var2=$(echo $gene | cut -c 3-)
  echo "$var is used for the awk pattern matching"
  echo "$var2 is used to match and create filenames"
  cd $gene
  echo "grepping $var"
  grep $var prokka_annotation.txt > ${var}_temp
  echo "awking $var"
  awk -F '\t|;' '{for(i=5;i<=NF;i++){if($i~/^gene=/){column=$i}} print $1,$3-$2,column}' ${var}_temp > gene_length_${var2}.txt

  #remove first column from contig length files "##sequence ID blabla"
  awk '{print $2,$3,$4}' contig_lengths.txt > contig_lengths_cut.txt 
  #both files can be used for all different genes, should be stored on above level

  #generating a file which combines gene length and contig length
  awk 'NR==FNR {contig[$1]=$3; next}                    
       $1 in contig {print $0, contig[$1]}' contig_lengths_cut.txt gene_length_${var2}.txt > ${var2}_contig_and_gene_length.txt
  cd ..
done
#what it does: for first file, use first column as index and third column as value	 
#then check if the first column of the file, first column