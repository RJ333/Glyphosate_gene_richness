#!/bin/sh

PATH_TO_DATA=/data/jwerner/glyphosate/IMP/ 
OUTPUT_DIR=/data/Rene/glyph/

# delete older files
rm $OUTPUT_DIR/contig_length.tsv

for name in A1 A2 A3 A4 A5 A6 A7 B8 B9 B10 
do 
  cd $PATH_TO_DATA/$name/output_IMP/Analysis
  echo $name
  # add all 10 samples to one list and add a column with the sample name  
  awk -v name="$name" 'BEGIN{FS = OFS = "\t"}{print name, $0}' mg.assembly.length.txt >> $OUTPUT_DIR/contig_length.tsv 
  cd $PATH_TO_DATA
done