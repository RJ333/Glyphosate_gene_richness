#!/bin/sh

# the samtools part of the script may take days, so use a tmux-session

## TODO: script that starts this script as TMUX session?
# https://gist.github.com/swaroopch/728896

# tmux 		# session will stay logged in
# tmux a 	# regain tmux session, must be logged in to the login node (sshbio-49)
# ctrl+b, c # new tab in tmux session for jobinfo etc
# ctrl+b, n # next
# ctrl+b, p # previous


PROKKA_DIR=/data/Rene/glyph/prokka
ORIGINAL_BASE_DIR=/data/jwerner/glyphosate/IMP
PROKKA_PATH=output_IMP/Analysis/annotation

# for samtools step
SAMTOOLS_BIN=/data/jwerner/tools/samtools-1.7/samtools
OUTPUT_DIR=$PROKKA_DIR/samtools
BED_DIR=${OUTPUT_DIR}/bed
RESULTS_DIR=${OUTPUT_DIR}/results
TMP_DIR=${OUTPUT_DIR}/tmp

mkdir -p $PROKKA_DIR
mkdir -p $OUTPUT_DIR
mkdir -p $BED_DIR
mkdir -p $RESULTS_DIR
mkdir -p $TMP_DIR
mkdir -p $OUTPUT_DIR/named

THRESHOLD=1

# this script removes any previous versions of the combined prokka file. then is collects prokka
# files from all samples, extracts the important columns and adds them to prokka_all_modified
# including the sample name as first column
# the output can then be used as input for R with splitstackshape 

rm -f $PROKKA_DIR/prokka_all_modified.tsv 

for sample in A1 A2 A3 A4 A5 A6 A7 B8 B9 B10
do
  echo $sample
  for prokka_annot in ${ORIGINAL_BASE_DIR}/${sample}/${PROKKA_PATH}/annotation.filt.gff
  do
	echo "processing $prokka_annot"
	# be careful with linebreks within awk
	awk -v sample=$sample \
    'BEGIN {OFS = FS = "\t"}
    {a=sample;} 
    {print a, $1, $4, $5, $9}' $prokka_annot >> $PROKKA_DIR/prokka_all_modified.tsv
    echo "prokka annotation of sample $sample modified and appended"
  done
done
echo "output written to $PROKKA_DIR/prokka_all_modified.tsv"


# the conda environment "Renv" needs to activated. So far it only works on bio-48 and only if 
# activated before running the script. Use:

# conda activate Renv 

# calling the R script:
Rscript /data/Rene/git/Glyphosate_gene_richness/reads_per_product_samtools/02_modify_prokka_products.r \
  -p="$PROKKA_DIR/prokka_all_modified.tsv" -o="$PROKKA_DIR" -t=$THRESHOLD


# how many cores should GNU Parallel use for samtools? 
# bio-48 has 24 cores, 
# bio-49 has 50 cores
GNU_CORES=15

# results directories

echo "starting samtools parallel"
samtools_view_parallel() {
	SAMPLE="$1"
	PRODUCT_NAME="$2"
	BED_FILENAME=${BED_DIR}/intersect_${PRODUCT_NAME}_${SAMPLE}.bed
	grep "$PRODUCT_NAME" ${PROKKA_DIR}/prokka_for_bed_${SAMPLE}.tsv \
	  > ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff
	
	cat ${TMP_DIR}/${PRODUCT_NAME}_${SAMPLE}_tmp.gff |\
	awk 'BEGIN {FS = OFS = "\t"} {print $1,$2,$3,$4,$7}' > ${BED_FILENAME}

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
export PROKKA_DIR

# not sure how to break this line?
parallel -j $GNU_CORES samtools_view_parallel ::: A1 A2 A3 A4 A5 A6 A7 B8 B9 B10 :::: $PROKKA_DIR/unified_unique_prokka_products_greater_${THRESHOLD}.tsv
echo "output written to $RESULTS_DIR"

echo "processing of samtools output"

# adds the name of the file in folder "results" as first column and copies to "named" folder 
cd $RESULTS_DIR
echo "adding file name as column to files"
for files in *.tsv
do	
  awk '{a=FILENAME;}{print a"\t"$0}' $files > $OUTPUT_DIR/named/named_${files}
done

cd $OUTPUT_DIR/named

echo "appending genes"
for nfiles in named_*.tsv
do	
  cat $nfiles >> appended_genes.tsv
done

echo "all information appended into appended_genes.tsv"
echo "cleaning up appended_genes.tsv and ..."

# BEGIN adds table header and sets Output field sep to tab
# split adresses specific column and splits on "_", storing the pieces in array a
# from this leftover, another split is performed to remove the ".tsv", stored in array b
# the respective fields are printed

echo "... formatting appended genes with awk"

awk 'BEGIN { OFS = "\t" ; print "contig\tgene\tsample\treads_per_gene"}
  {
  split ($1, a, "_") 
  split (a[3], b, "\\.") 
  print $2, a[2], b[1], $3
  }' appended_genes*tsv > contig_gene_sample_reads.tsv
  
rm appended_genes*tsv

echo "cleaning is done, output written to $OUTPUT_DIR/named/contig_gene_sample_reads.tsv"