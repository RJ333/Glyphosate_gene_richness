#!/usr/bin/bash

## NOT TESTED YET!!  use a small part of the gene list for testing

# this script may take days, so use a tmux-session

# tmux 		# session will stay logged in
# tmux a 	# regain tmux session, must be logged in to the login node (sshbio-49)
# ctrl+b, c # new tab in tmux session for jobinfo etc
# ctrl+b, n # next
# ctrl+b, p # previous

# how many cores should GNU Parallel use for samtools?
GNU_CORES=45

# directories of original data
ORIGINAL_BASE_DIR=/data/jwerner/glyphosate/IMP/
PROKKA_DIR=/data/Rene/glyph/prokka 
SAMTOOLS_BIN=/data/jwerner/tools/samtools-1.7/samtools

# working directory
OUTPUT_DIR=/data/Rene/glyph/prokka/samtools  

# results directories
BED_DIR=${OUTPUT_DIR}/bed
RESULTS_DIR=${OUTPUT_DIR}/results
TMP_DIR=${OUTPUT_DIR}/tmp

mkdir -p $OUTPUT_DIR
mkdir -p $BED_DIR
mkdir -p $RESULTS_DIR
mkdir -p $TMP_DIR

samtools_view_parallel() {
	SAMPLE="$1"
	cd $OUTPUT_DIR
	PRODUCT_NAME="$2"
	BED_FILENAME=${BED_DIR}/intersect_${PRODUCT_NAME}_${SAMPLE}.bed

	grep "$PRODUCT_NAME" $PROKKA_DIR/prokka_for_bed_${SAMPLE}.tsv \
	  > ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff

	cat ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff |\
	awk -F 'BEGIN {FS = OFS = "\t"} NR > 1 {print $1,$2,$3,$4,$7}' > ${BED_FILENAME}

	$SAMTOOLS_BIN view -L ${BED_FILENAME} ${ORIGINAL_BASE_DIR}/$SAMPLE/output_IMP/Assembly/mg.reads.sorted.bam |\
      grep -v -P "^\@" | cut -f 1,3 | sort | uniq | cut -f 2  | sort | uniq -c |\
	  perl -ane '$_=~/^\s+(\d+) (.+)$/;chomp($2); print "$2\t$1\n"; '\
	  > ${RESULTS_DIR}/mg.reads.per.gene_${PRODUCT_NAME}_${SAMPLE}.tsv

	rm ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff
}


export -f samtools_view_parallel
# These variables are defined outside the function and must be exported to be visible
export BED_DIR
export TMP_DIR
export OUTPUT_DIR
export RESULTS_DIR
export SAMTOOLS_BIN
export ORIGINAL_BASE_DIR

parallel -j $GNU_CORES samtools_view_parallel ::: A1 A2 A3 A4 A5 A6 A7 B8 B9 B10 :::: $PROKKA_DIR/test_products.tsv
