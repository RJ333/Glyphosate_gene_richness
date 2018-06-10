
## NOT TESTED YET!!

#!/usr/bin/bash

# this script may take days, so use a tmux-session

# tmux 		# session will stay logged in
# tmux a 	# regain tmux session, must be logged in to the login node (sshbio-49)
# ctrl+b, c # new tab in tmux session for jobinfo etc
# ctrl+b, n # next
# ctrl+b, p # previous

# how often needs a product be present per sample
THRESHOLD=1
# how many cores should GNU Parallel use for samtools?
GNU_CORES=44

# directories of original data
ORIGINAL_BASE_DIR=/data/jwerner/glyphosate/IMP/
ORIGINAL_PROKKA_DIR=/data/jwerner/glyphosate/IMP/${sample}/output_IMP/Analysis/annotation/ 
SAMTOOLS_BIN=/data/jwerner/tools/samtools-1.7/samtools

# working directory
OUTPUT_DIR=/data/Rene/glyph/prokka  

# results directories
BED_DIR=${OUTPUT_DIR}/bed
RESULTS_DIR=${OUTPUT_DIR}/results
TMP_DIR=${OUTPUT_DIR}/tmp

mkdir -p $BED_DIR
mkdir -p $RESULTS_DIR

## include sample vector?
## include unique list argument?


# get all product names from annotation in lower case and combine them with the original prokka file:
cd $OUTPUT_DIR
rm $OUTPUT_DIR/unified_all_samples.tsv
for sample in $(find . -mindepth 1 -maxdepth 1 -type d | sed "s|^\./||")
do
  for prokka_annot in ${OUTPUT_DIR}/${sample}/annotation.filt.gff
  do
    awk -v SAMPLE=$sample \
	  '\
	  BEGIN { FS = "\t|;"; OFS = "\t"}\
	  {match($0,/product=[^;]*/); product_value=substr($0,RSTART,RLENGTH); {print tolower(product_value), SAMPLE}}' $prokka_annot | tr '/()[:blank:]-_' '@@@@@@'  > $OUTPUT_DIR/unified_${sample}.tsv
	  # pasting the modified products column to the original prokka file
	  paste -d '\t' $prokka_annot $OUTPUT_DIR/unified_${sample}.tsv > $OUTPUT_DIR/${sample}/${sample}_prokka_modified.tsv
	  cat $OUTPUT_DIR/unified_${sample}.tsv >> $OUTPUT_DIR/unified_all_samples.tsv
	  # remove the for-pasting tables
	  rm $OUTPUT_DIR/unified_${sample}.tsv
	  rm $prokka_annot
  done
done

# generate a combined prokka modified output for all samples
cd $OUTPUT_DIR
for prokka_table in *_prokka_modified.tsv
do
  sample=$(echo $prokka_table | sed "s/_.*$//")
  awk -v sample=$sample '{a = sample;} {print a"\t"$0}' $prokka_table >> $OUTPUT_DIR/all_prokka_modified.gff
  echo "$sample appended"
done
cat $OUTPUT_DIR/all_prokka_modified.gff | rev | sed 's/^@//' | rev > $OUTPUT_DIR/all_prokka_modified2.gff
rm $OUTPUT_DIR/all_prokka_modified.gff

# use the list of products from all samples as list of queries for samtools. remove duplicates and set minimum threshold
cd $OUTPUT_DIR
sort unified_all_samples.tsv| rev | sed 's/@/\t/' | rev | cut -f 1 |\
sort| uniq -c| awk -v THRESHOLD=$THRESHOLD '$1 > THRESHOLD {print $0}' |\
awk '{split($2,a,"="); print a[2]}' > unified_unique_prokka_products_greater_${THRESHOLD}.tsv

## sam_parallel shouldn't use original prokka tables anymore --> R modified ones
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
