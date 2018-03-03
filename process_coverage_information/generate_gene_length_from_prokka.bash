#!/bin/bash

#extract gene length from prokka annotation 

genename="phnM"

awk -F '\t|;' -v gene="$genename" ' $0 ~ gene {print $3-$2,$6}' 

awk -F '\t|;' ' $0 ~ "phnM" {print $1,$3-$2,$6}'  prokka_annotation.txt > gene_length_phnM


#extract contig length from contig_length and match with gene_length
awk '{print $2,$3,$4}' contig_lengths.txt > contig_lengths_cut.txt


awk 'NR==FNR { contig[NR]=$1; contig_length[$1]=$3; next}; (for k in contig) if (0){print $0 "," contig[$1]}' contig_lengths_cut.txt gene_length_phnM

awk 'NR==FNR { contig[$1]=$3; next}; {for (k in contig) if ($3 ~ contig[k]) print $0, contig[$1] | "uniq"}' contig_lengths_cut.txt gene_length_phnM

awk 'NR==FNR { contig[$1]=$3; next}; {for (k in contig) if ($3~ contig[k]) print $0, contig[$1] }' contig_length_test.txt gene_length_test.txt
###working:
awk 'NR==FNR {contig[$1]=$3; next}                    
     $1 in contig {print $0, contig[$1]}' contig_lengths_cut.txt gene_length.txt	

#what it means: for first file, use first column as index and third column as value	 
#then check if the first column of the file, first column
	
	#this script adds the annotated gene region the coverage was calculated for. 
#It became complicated as it was supposed to distuingish between two genes placed on the same contig, only differing by their annotated region.
#it now continously counts all unique occurrences of the gene based on the data of the bed file
#####but I still don't exactly understand how it does it!!!
#
#with a lot of help from
#https://stackoverflow.com/questions/49050245/how-to-use-awk-to-add-specific-values-to-a-column-based-on-numeric-ranges/49054710?noredirect=1#comment85120455_49054710

awk 'NR==FNR{start[NR]=$2; end[NR]=$3; key[$1,$2]=$4 sprintf("_%03d",NR); next}
           {for(i in start)
              {s=start[i];
               if(s<=$2 && $2<=end[i] && ($1,s) in key) print $0,key[$1,s] | "sort -nk1 "}}' intersect.bed phnm_combined.txt_filtered_5.txt > test.cov
			   
	

awk -F',' 'NR==FNR{label[$1]=$1;date[$1]=$2;next}; ($2==label[$2]){print $0 "," date[$2]}' <(sort -k1 file2.csv) <(sort -k2 file1.csv) &> file3.csv
Explanation

    -F',': 						#Sets the field separator to ,
    NR==FNR : 					#NR is the current input line number and FNR the current file's line number. The two will be equal only while the 1st file is being read.
    label[$1]=$1: 				#Save column 1 (label) from the first file (argument) file2.csv in hash-array using column 1 as the key
    date[$1]=$2: 				#Save column 2 (date) from the first file file2.csv in hash-array using column 1 as the key
    next: 						#Then, skip to the next line so that this is only applied on the 1st file.
    ($2==label[$2]): 			#the else block will only be executed if this is the second file (file1.csv), so we check whether field 2 of this file is in array label with the field $2 as the key ($2==label[$2]).
    {print $0 "," date[$2]}: 	#If that's true print the entire file1.csv and append the date column from file2.csv.
    sort -k1 file2.csv: 		#Sorts file2.csv column 1 by setting k1. k2 sorts column 2

	
#File 1
dbID    labnumber    myID    Status
CMV_1235    LAB06    56-1    Fail
CMV_1236    LAB14    57-1    Fail
CMV_2137    LAB84    54-4    Pass
CMV_2238    LAB85    50-3    
CMV_C131    LAB21    51-2    Pass

#File 2
labnumber    date
LAB06    18/01/2016
LAB14    27/04/2016
LAB18    10/01/2016
LAB21    9/02/2016
LAB69    4/03/2016
LAB84    18/02/2016
LAB22    18/03/2016
LAB85    27/03/2016

(Not totally overlapping: there may be samples in file 1 but not file 2 and vice versa)

I want to print to file 3:
dbID    labnumber    myID    Status    date
CMV_1235    LAB06    56-1    Fail      18/01/2016
CMV_1236    LAB14    57-1    Fail      27/04/2016
CMV_2137    LAB84    54-4    Pass      18/02/2016
CMV_2238 LAB85 50-3 27/03/2016
