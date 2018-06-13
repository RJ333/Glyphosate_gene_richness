
## NOT TESTED YET!!

#!/usr/bin/bash

# this script may take days, so use a tmux-session

# tmux 		# session will stay logged in
# tmux a 	# regain tmux session, must be logged in to the login node (sshbio-49)
# ctrl+b, c # new tab in tmux session for jobinfo etc
# ctrl+b, n # next
# ctrl+b, p # previous

# how many cores should GNU Parallel use for samtools?
GNU_CORES=44

# directories of original data
ORIGINAL_BASE_DIR=/data/jwerner/glyphosate/IMP/
ORIGINAL_PROKKA_DIR=  #/data/jwerner/glyphosate/IMP/${sample}/output_IMP/Analysis/annotation/ 
SAMTOOLS_BIN=/data/jwerner/tools/samtools-1.7/samtools

# working directory
OUTPUT_DIR=  # /data/Rene/glyph/prokka  

# results directories
BED_DIR=${OUTPUT_DIR}/bed
RESULTS_DIR=${OUTPUT_DIR}/results
TMP_DIR=${OUTPUT_DIR}/tmp

mkdir -p $BED_DIR
mkdir -p $RESULTS_DIR


## TODO:

## include sample vector?
## include unique list argument?
## sam_parallel shouldn't use original prokka tables anymore --> R modified ones
# check if all variables are needed
# which columns are the needed ones?

sam_parallel() {
	SAMPLE="$1"
	cd $OUTPUT_DIR
	PRODUCT_NAME_NO_SPACES="$2"
	BED_FILENAME=${BED_DIR}/intersect_${PRODUCT_NAME_NO_SPACES}_${SAMPLE}.bed

	grep "$PRODUCT_NAME_NO_SPACES" $SAMPLE/{SAMPLE}_prokka_modified.tsv > ${TMP_DIR}/${PRODUCT_NAME_NO_SPACES}_${SAMPLE}_tmp.gff

	cat ${TMP_DIR}/${PRODUCT_NAME_NO_SPACES}_${SAMPLE}_tmp.gff | \
	awk -F '\t|;' '{for(i=9;i<=NF;i++){if($i~/^product=/){column=$i}} print $1,$4,$5,$5-$4,column}' > ${BED_FILENAME}

	$SAMTOOLS_BIN view -L ${BED_FILENAME} ${ORIGINAL_BASE_DIR}/$SAMPLE/output_IMP/Assembly/mg.reads.sorted.bam | grep -v -P "^\@" | cut -f 1,3 | sort |\
	  uniq | cut -f 2  | sort | uniq -c | perl -ane '$_=~/^\s+(\d+) (.+)$/;chomp($2); print "$2\t$1\n"; '\
	  > ${RESULTS_DIR}/mg.reads.per.gene_${PRODUCT_NAME_NO_SPACES}_${SAMPLE}.tsv

	rm ${TMP_DIR}/${PRODUCT_NAME_NO_SPACES}_${SAMPLE}_tmp.gff
}


export -f sam_parallel
# These variables are defined outside the function and must be exported to be visible
export BED_DIR
export TMP_DIR
export OUTPUT_DIR
export RESULTS_DIR
export SAMTOOLS_BIN
export ORIGINAL_BASE_DIR

parallel -j $GNU_CORES sam_parallel ::: A1 A2 A3 A4 A5 A6 A7 B8 B9 B10 :::: $OUTPUT_DIR/unified_unique_prokka_products_greater_${THRESHOLD}.tsv
