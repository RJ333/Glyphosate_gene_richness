#!/usr/bin/bash

# directories of original data
ORIGINAL_BASE_DIR=/data/jwerner/glyphosate/IMP
PROKKA_DIR=/data/Rene/glyph/prokka 
SAMTOOLS_BIN=/data/jwerner/tools/samtools-1.7/samtools
echo "input directories set"

# working directory
OUTPUT_DIR=/data/Rene/glyph/prokka/samtools_test
echo "output directory $OUTPUT_DIR set"

# results directories
BED_DIR=${OUTPUT_DIR}/bed
RESULTS_DIR=${OUTPUT_DIR}/results
TMP_DIR=${OUTPUT_DIR}/tmp

mkdir -p $OUTPUT_DIR
mkdir -p $BED_DIR
mkdir -p $RESULTS_DIR
mkdir -p $TMP_DIR

echo "result dirs created"

SAMPLE="A1"
PRODUCT_NAME="?@d@glucose@1@phosphatase"
BED_FILENAME=${BED_DIR}/intersect_${PRODUCT_NAME}_${SAMPLE}.bed
echo "$SAMPLE"
echo "$PRODUCT_NAME"
echo "$BED_FILENAME"

	grep "$PRODUCT_NAME" ${PROKKA_DIR}/prokka_for_bed_${SAMPLE}.tsv 
	
	echo "grep output"?
	
	grep "$PRODUCT_NAME" ${PROKKA_DIR}/prokka_for_bed_${SAMPLE}.tsv \
	  > ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff
	
	echo "grepped product in ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff"
	
	cat ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff |\
	awk 'BEGIN {FS = OFS = "\t"} {print $1,$2,$3,$4,$7}'
	
	echo "awk output?"
	
	cat ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff |\
	awk 'BEGIN {FS = OFS = "\t"} {print $1,$2,$3,$4,$7}' > ${BED_FILENAME}
	
	echo "output written to ${BED_FILENAME}"
	echo "starting samtools step with ${ORIGINAL_BASE_DIR}/$SAMPLE/output_IMP/Assembly/mg.reads.sorted.bam"
	
	$SAMTOOLS_BIN view -L ${BED_FILENAME} ${ORIGINAL_BASE_DIR}/$SAMPLE/output_IMP/Assembly/mg.reads.sorted.bam |\
      grep -v -P "^\@" | cut -f 1,3 | sort | uniq | cut -f 2  | sort | uniq -c |\
	  perl -ane '$_=~/^\s+(\d+) (.+)$/;chomp($2); print "$2\t$1\n"; '\
	  > ${RESULTS_DIR}/mg.reads.per.gene_${PRODUCT_NAME}_${SAMPLE}.tsv
	
	echo "samtools output written to ${RESULTS_DIR}/mg.reads.per.gene_${PRODUCT_NAME}_${SAMPLE}.tsv"
	

	rm ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff
	echo "${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff removed"