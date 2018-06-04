#!/bin/sh

# add contig length and average coverage of all glyphosate samples to one list and add a column with the sample name  
BASE_DIR=/data/jwerner/glyphosate/IMP/
OUTPUT_DIR=/data/Rene/glyph/

cd $BASE_DIR
# removing old versions if existing
rm ${OUTPUT_DIR}/contig_length_coverage.tsv 
for contig_coverage in */metagenome_coverage/covstats.txt
do
  sample=$(echo $contig_coverage |  sed "s/\/.*$//")
  echo "processing $sample"
  # selecting columns with contig_id, average coverage and contig_length 
  awk -v SAMPLE=$sample 'BEGIN {FS = OFS= "\t"} NR > 1 {print SAMPLE, $1, $2, $3}' $contig_coverage >> ${OUTPUT_DIR}/contig_length_coverage.tsv
done
echo "output contig_length_coverage.tsv generated in $OUTPUT_DIR"  