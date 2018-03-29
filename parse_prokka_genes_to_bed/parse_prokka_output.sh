#!/usr/bin/env bash

BASE_DIR=/data/jwerner/glyphosate/IMP
SAMTOOLS_BIN=/data/jwerner/tools/samtools-1.7/samtools # path on bio-49

usage()
{
cat << EOF
usage:
$0 -o <OUTPUT_DIR>

The script creates a bed file based on prokka annotations.

OPTIONS:
   -o          Output directory (mandatory)

   -h/--help   Show this message

EOF
}

while getopts ":o:" o; do
    case "${o}" in
        o)
            OUTPUT_DIR=${OPTARG}
            ;;
        *)
            # usage
            ;;
    esac
done
shift $((OPTIND-1))

if [[ -z $OUTPUT_DIR ]]
then
    echo Parameter -o is not set. This parameter is required. Exit code: 1. Exiting ...
    echo -en "\n"
    usage
    echo -en "\n"
    exit 1
fi

if [ -d $OUTPUT_DIR ]
then
	echo -en "\n"
    echo Output directory already exists. Aborting to avoid overriding. Exit code: 2. Exiting ...
    echo -en "\n"
    usage
    echo -en "\n"
    exit 2
else
    BED_DIR=${OUTPUT_DIR}/bed
    RESULTS_DIR=${OUTPUT_DIR}/results
    TMP_DIR=${OUTPUT_DIR}/tmp
    mkdir -p $BED_DIR
    mkdir -p $RESULTS_DIR
    mkdir -p $TMP_DIR
fi

cd $BASE_DIR

#for SAMPLE_FOLDER in *
#do
#    cd ${SAMPLE_FOLDER}/output_IMP/Analysis/annotation
##    for GENE_NAME in gyrA gyrB gox phnA phnC phnD phnE phnF phnG phnH phnI phnJ phnK phnL phnM phnN phnX phnW phnR \
##      phnS phnT phnU phnP pspE pstC purU purB pphA phoR phoB phoP soxA soxB thiO recA rpoC ktrB ydiF nuoF mntB araD \
##      tfdB artI arfA malF qedA mlhB uvrB fdhA fdm recG
##     for GENE_NAME in norB nirS nosZ soxR soxS soxC katG* katE czcD merA tmoS chiA* chtA xynB bcsZ celB cenC
#    for GENE_NAME in dnaJ dnaK dmlR ftsY ftsZ rplB polA gap* sdhA
#    do
#        GREP_GENE_NAME="${GENE_NAME}_"
#        BED_FILENAME=${BED_DIR}/intersect_${GENE_NAME}_${SAMPLE_FOLDER}.bed
#        grep $GREP_GENE_NAME annotation.filt.gff > ${TMP_DIR}/tmp.gff
#
#        cat ${TMP_DIR}/tmp.gff | \
#          awk -F '\t|;' '{for(i=9;i<=NF;i++){if($i~/^gene=/){column=$i}} print $1,$4,$5,$5-$4,column}' > ${BED_FILENAME}
#
#        $SAMTOOLS_BIN view -L ${BED_FILENAME} ../../Assembly/mg.reads.sorted.bam | grep -v -P "^\@" | cut -f 1,3 | \
#          sort | uniq | cut -f 2  | sort | uniq -c | perl -ane '$_=~/^\s+(\d+) (.+)$/;chomp($2); print "$2\t$1\n"; ' \
#          > ${RESULTS_DIR}/mg.reads.per.gene_${GENE_NAME}_${SAMPLE_FOLDER}.tsv &
#
#        rm ${TMP_DIR}/tmp.gff
#    done
#    cd ../../../..
#done

for SAMPLE_FOLDER in *
do
    cd ${SAMPLE_FOLDER}/output_IMP/Analysis/annotation
    # for loops not recommended to read lines from file, set IFS= for while loop, if you have e.g. spaces per line
	cat /data/Rene/test_products.tsv | while read PRODUCT_NAME_NO_SPACES 
    do
	# tr only works on single characters. use echo '§' | xxd -c 1 to check possible delimiter
	PRODUCT_NAME=`echo ${PRODUCT_NAME_NO_SPACES} | tr "@" " "`
        BED_FILENAME=${BED_DIR}/intersect_${PRODUCT_NAME_NO_SPACES}_${SAMPLE_FOLDER}.bed
        grep "$PRODUCT_NAME" annotation.filt.gff > ${TMP_DIR}/tmp.gff

        cat ${TMP_DIR}/tmp.gff | \
          awk -F '\t|;' '{for(i=9;i<=NF;i++){if($i~/^product=/){column=$i}} print $1,$4,$5,$5-$4,column}' > ${BED_FILENAME}

        parallel -j 25 -k $SAMTOOLS_BIN view -L ${BED_FILENAME} ../../Assembly/mg.reads.sorted.bam | grep -v -P "^\@" | cut -f 1,3 | \
          sort | uniq | cut -f 2  | sort | uniq -c | perl -ane '$_=~/^\s+(\d+) (.+)$/;chomp($2); print "$2\t$1\n"; ' \
          > ${RESULTS_DIR}/mg.reads.per.gene_${PRODUCT_NAME_NO_SPACES}_${SAMPLE_FOLDER}.tsv &
          ### TODO (WARNING): currently all processes are forked at the same time. this needs to be resolved. limit to a certain number of processes.

        rm ${TMP_DIR}/tmp.gff
    done
    cd ../../../..
done

### TODO: to map product name back to genes: all lower case, unify white spaces, dashes and  !!  


rmdir -p ${OUTPUT_DIR}/tmp

