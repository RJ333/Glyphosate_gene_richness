#!/bin/bash

#extract gene length from prokka annotation 

genename="phnM"

awk -F '\t|;' -v a = "$genename" ' $0 ~ a {print $1,$3-$2,$6}'  prokka_annotation.txt > gene_length_${genename}.txt

#remove first column from contig length files "##sequence ID blabla"
awk '{print $2,$3,$4}' contig_lengths.txt > contig_lengths_cut.txt #both files can be used for all different genes, should be stored on above level

#generating a file which combines gene length and contig length
awk 'NR==FNR {contig[$1]=$3; next}                    
     $1 in contig {print $0, contig[$1]}' contig_lengths_cut.txt gene_length_${genename}.txt > ${genename}_contig_and_gene_length.txt

#what it does: for first file, use first column as index and third column as value	 
#then check if the first column of the file, first column
