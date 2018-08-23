#!/usr/bin/bash

# set up the general folder structure and get your sample data

# directories of original data
BASE_DIR=/data/Rene/glyph/tree/
REF_DIR=${BASE_DIR}/uniprot_sequences
GIT_REPO=/data/Rene/Glyphosate_gene_richness/tree_analysis
SUB_DIR=${BASE_DIR}/trees_degradation
PROKKA_DIR=/data/jwerner/glyphosate/IMP/${sample}/output_IMP/Analysis/annotation/

mkdir -p ${BASE_DIR}/prokka_original/
mkdir -p ${BASE_DIR}/prokka_modified/

# copy own data into tree working folder
for SAMPLE in A1 A2 A3 A4 A5 A6 A7 B8 B9 B10
do 
    echo $SAMPLE
	cp ${PROKKA_DIR}/prokka.faa ${BASE_DIR}/prokka_original/${SAMPLE}_prokka.faa
done

# original Prokka data was put into a single line, 
# sample names were included in header and whitespaces replaced with @

for SAMPLE in A1 A2 A3 A4 A5 A6 A7 B8 B9 B10
do 
    echo $SAMPLE
	awk '{if (!/>/) {printf "%s",$0;next} else {printf "%s%s%s","\n",$0,"\n"}}' \
	  prokka_original/${SAMPLE}_prokka.faa | \
	  tail -n +2 | \
	  sed "s/ /_${SAMPLE} /" | \
	  sed "s/ /@/g" > prokka_modified/${SAMPLE}_prokka.oneline.sampleheader.nospaces.faa
done