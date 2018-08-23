#!/usr/bin/bash

# set up the general folder structure and get your sample data

# directories of original data
BASE_DIR=/data/Rene/glyph/tree/
PROKKA_DIR=/data/jwerner/glyphosate/IMP/${SAMPLE}/output_IMP/Analysis/annotation/

SAMPLE_ARRAY=(A1 A2 A3 A4 A5 A6 A7 B8 B9 B10)

mkdir -p ${BASE_DIR}/prokka_original/
mkdir -p ${BASE_DIR}/prokka_modified/

# copy own data into tree base folder
for SAMPLE in "${SAMPLE_ARRAY[@]}"
do 
    echo $SAMPLE
	cp /data/jwerner/glyphosate/IMP/${SAMPLE}/output_IMP/Analysis/annotation/prokka.faa \
	  ${BASE_DIR}/prokka_original/${SAMPLE}_prokka.faa
done

echo "prokka files have been copied"

# original Prokka data was put into a single line, 
# sample names were included in header and whitespaces replaced with @

for SAMPLE in "${SAMPLE_ARRAY[@]}"
do 
    echo $SAMPLE
	awk '{if (!/>/) {printf "%s",$0;next} else {printf "%s%s%s","\n",$0,"\n"}}' \
	  prokka_original/${SAMPLE}_prokka.faa | \
	  tail -n +2 | \
	  sed "s/ /_${SAMPLE} /" | \
	  sed "s/ /@/g" > prokka_modified/${SAMPLE}_prokka.oneline.sampleheader.nospaces.faa
done

echo "prokka files have been modified"