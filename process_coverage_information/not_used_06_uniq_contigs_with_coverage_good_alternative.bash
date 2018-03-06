#!/bin/sh

#counting occurrences of $gene and $contigs in columns 5 and 1 respectively in all cov files found
#printing information including previously applied parameters
#all output except for the last input gene will be printed by  FNR==1 { prt() }, last input file uses END { prt() }
#adding step to calculate duplicates? not tested yet
PATH_TO_DATA=/drives/d/1_Fachliches/gene_richness_selected_genes

cd $PATH_TO_DATA
for gene in $(find . -mindepth 1 -maxdepth 1 -type d) 
do 
  cd $gene
  echo "start counting for $gene"
  awk ' function prt() { if (fname != "") {
      printf "gene was found %d times on %d contigs with minimum coverage of %d in %s",       
	  length(genes), length(contigs), threshold, fname
      delete genes
      delete contigs
      }
    fname = FILENAME
  }
  
  BEGIN {
    threshold = "5"
  }
  FNR==1 { prt() }
  {
    genes[$5]
    contigs[$1]
  }
  END { prt() }' *counts.cov > 06_${gene}_on_contig_counts.txt
  cd ..
done


###thanks to Ed Morton on Stack Overflow
#https://stackoverflow.com/questions/49067678/print-output-of-user-defined-function-in-awk-gives-unexpected-token-error