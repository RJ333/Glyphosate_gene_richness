#!/bin/sh

# extract gene length from prokka annotation 

# outer loop to get sample names
cd /data/jwerner/glyphosate/IMP/
for sample in $(find . -mindepth 1 -maxdepth 1 -type d | sed "s|^\./||")
do
  echo $sample
  for prokka_annot in /data/jwerner/glyphosate/IMP/$sample/output_IMP/Analysis/annotation/annotation.filt.gff
  do
    echo "$prokka_annot is stored in the for loop"
    for genes in phnM
	do
	  grep "$genes" $prokka_annot |\
	  #awk -F '\t|;' '{$1=$1;print $0}'
	  #awk -F '\t|;' -v SAMPLE=$sample -v GENES=$genes '{for(i=5;i<=NF;i++){if($i~/^gene=/){column=$i}} print $1,$5-$4,column,GENES,SAMPLE}'
	  #awk -v SAMPLE=$sample -v GENES=$genes 'BEGIN { FS = "\t|;" ; OFS = "\t" }' '{match($0,/gene=[^;]*/); gene_value=substr($0,RSTART,RLENGTH); match($0,/eC_number=[^;]*/); eC_number_value=substr($0,RSTART,RLENGTH)} END {print $1,$5-$4,gene_value,eC_number_value,GENES,SAMPLE}'
      awk -v SAMPLE=$sample -v GENES=$genes'{match($0,/gene=[^;]*/); gene_value=substr($0,RSTART,RLENGTH); match($0,/eC_number=[^;]*/); eC_number_value=substr($0,RSTART,RLENGTH); print(gene_value,eC_number_value)}' 
	done
  done
done

## there are some problems with the last awk step at the moment, it is not printing what is it supposed to print
	
	
  #var=$(echo $genes | cut -c 3- |awk -vFS= -vOFS= '{$NF=toupper($NF)}1')
 # var2=$(echo $genes | cut -c 3-)
  #echo "$var is used for the awk pattern matching"
 # echo "$var2 is used to match and create filenames"
 # cd $genes
 # echo "grepping $var"
 # grep $var genesation.txt > ${var}_temp
 # echo "awking $var"
 # awk -F '\t|;' '{for(i=5;i<=NF;i++){if($i~/^genes=/){column=$i}} print $1,$3-$2,column}' ${var}_temp > prokka_annot_length_${var2}.txt

  #generating a file which combines gene length and contig length
 # awk 'NR==FNR {contig[$1]=$3; next}                    
      # $1 in contig {print $0, contig[$1]}' contig_lengths_cut.txt gene_length_${var2}.txt > ${var2}_contig_and_gene_length.txt
  #cd ..
#done
#what it does: for first file, use first column as index and third column as value	 
#then check if the first column of the file, first column