#!/bin/sh

OUTPUT_DIR=/data/Rene/glyph/prokka
ORIGINAL_BASE_DIR=/data/jwerner/glyphosate/IMP
PROKKA_PATH=output_IMP/Analysis/annotation


# this script removes any previous versions of the combined prokka file. then is collects prokka
# files from all samples, extracts the important columns and adds them to prokka_all_modified
# including the sample name as first column
# the output can then be used as input for R with splitstackshape 
mkdir -p $OUTPUT_DIR
rm -f $OUTPUT_DIR/prokka_all_modified.tsv 

for sample in A1 A2 A3 A4 A5 A6 A7 B8 B9 B10
do
  echo $sample
  for prokka_annot in ${ORIGINAL_BASE_DIR}/${sample}/${PROKKA_PATH}/annotation.filt.gff
  do
	echo $prokka_annot
	# be careful with linebreks within awk
	awk -v sample=$sample \
    'BEGIN {OFS = FS = "\t"}
    {a=sample;} 
    {print a, $1, $4, $5, $9}' $prokka_annot >> $OUTPUT_DIR/prokka_all_modified.tsv
    echo "prokka annotation of sample $sample modified and appended"
  done
done
echo "output written to $OUTPUT_DIR/prokka_all_modified.tsv"