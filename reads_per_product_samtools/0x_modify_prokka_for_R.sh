#!/bin/sh


OUTPUT_DIR=/data/Rene/glyph/prokka  
rm $OUTPUT_DIR/unified_all_samples.tsv

# path to prokka files:
#/data/jwerner/glyphosate/IMP/*/output_IMP/Analysis/annotation

# create output folder
#mkdir -p /data/Rene/glyph/prokka

# go to folder with the different samples
#cd /data/jwerner/glyphosate/IMP

# put for loop in a function to simply copy into the terminal
# function adds a column with the sample name and saves the modified file somewhere else
#function hihi {
#for sample in *
#do
#  echo $sample
#  cd /data/jwerner/glyphosate/IMP/$sample/output_IMP/Analysis/annotation
#  awk -v sample=$sample '{a=sample;}{print a"\t"$0}' annotation.filt.gff > /data/Rene/glyph/prokka/${sample}_prokka.gff
#done
#}
dieser schritt muss angepasst werden, nachdem ich prokka daten verändere 
--> parse_prokka_genes_to_bed/* inkl der GNU parallel version des Skriptes checken

# run function
hihi

# append the 10 modified prokka files to one list

cd /data/Rene/glyph

function hoho {
for prokkas in *_prokka.gff
do	
  cat $prokkas >> all_prokkas.gff
  echo "appending $prokkas done"
done
}

# run function
hoho

die folgenden Schritte sollten nicht mehr in Excel durchgeführt werden, 
viele Ersetzungen werden nun vorher schon durchgeführt

# resulting prokka file with more than 1.200.000 lines is to big for excel. for convenience, some cleaning tasks should be tested in excel
# before writing code. So, the file is split roughly in half, creating two new files (xaa and xab)

split -l 700000 all_prokkas.gff

mv xaa all_prokkas_part1.gff
mv xab all_prokkas_part2.gff

# in excel
add headers to the columns: sample contig_id start_pos end_pos gene_length eC_number gene product (also including note=...) products_adjusted
add column sample_contig_id to uniquely identify contigs from different samples
add sample_contig_id_product for same reason
split single column with annotations based on ";" as sep
order annotation columns, use "gene" etc as header
delete columns with prodigal, CDS, "0", "+/-", inference, locus tag, PROKKA id
add column gene length
remove counts from genE_01, genE_02 etc --> genE
replace in product column > adjusted_products: 
  "/", "-", "_","'", """ " and whitespace with @, 
  everything lower case (=klein())
  (or later in R, if forgotten in excel)