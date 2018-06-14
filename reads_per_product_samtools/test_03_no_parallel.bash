#!/usr/bin/bash

# directories of original data
ORIGINAL_BASE_DIR=/data/jwerner/glyphosate/IMP/
PROKKA_DIR=/data/Rene/glyph/prokka 
SAMTOOLS_BIN=/data/jwerner/tools/samtools-1.7/samtools

# working directory
OUTPUT_DIR=/data/Rene/glyph/prokka/samtools_test

# results directories
BED_DIR=${OUTPUT_DIR}/bed
RESULTS_DIR=${OUTPUT_DIR}/results
TMP_DIR=${OUTPUT_DIR}/tmp

mkdir -p $OUTPUT_DIR
mkdir -p $BED_DIR
mkdir -p $RESULTS_DIR
mkdir -p $TMP_DIR

SAMPLE=A1
	PRODUCT_NAME="?@d@glucose@1@phosphatase"
	BED_FILENAME=${BED_DIR}/intersect_${PRODUCT_NAME}_${SAMPLE}.bed

	grep "$PRODUCT_NAME" ${PROKKA_DIR}/prokka_for_bed_${SAMPLE}.tsv \
	  > ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff

	cat ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff |\
	awk -F 'BEGIN {FS = OFS = "\t"} NR > 1 {print $1,$2,$3,$4,$7}' > ${BED_FILENAME}

	$SAMTOOLS_BIN view -L ${BED_FILENAME} ${ORIGINAL_BASE_DIR}/$SAMPLE/output_IMP/Assembly/mg.reads.sorted.bam |\
      grep -v -P "^\@" | cut -f 1,3 | sort | uniq | cut -f 2  | sort | uniq -c |\
	  perl -ane '$_=~/^\s+(\d+) (.+)$/;chomp($2); print "$2\t$1\n"; '\
	  > ${RESULTS_DIR}/mg.reads.per.gene_${PRODUCT_NAME}_${SAMPLE}.tsv

	rm ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff