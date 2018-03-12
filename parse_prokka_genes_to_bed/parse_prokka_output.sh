#!/usr/bin/env bash

BASE_DIR=/data/jwerner/glyphosate/IMP
SAMTOOLS_BIN=/data/jwerner/tools/samtools-1.7/samtools # path on bio-49

usage()
{
cat << EOF
usage:
$0 -i <INPUT_DIR> -o <OUTPUT_DIR>

The script creates a bed file based on prokka annotations.

OPTIONS:
   -i          Input directory (prokka annotation directory, mandatory)
   -o          Output directory (mandatory)

   -h/--help   Show this message

EOF
}

while getopts ":i:o:" o; do
    case "${o}" in
        i)
            INPUT_DIR=${OPTARG}
            ;;
        o)
            OUTPUT_DIR=${OPTARG}
            ;;
        *)
            # usage
            ;;
    esac
done
shift $((OPTIND-1))

if [[ -z $INPUT_DIR ]] && [[ -z $OUTPUT_DIR ]]
then
    echo Parameter -i and or -o is not set. These parameters are required. Exit code: 1. Exiting ...
    echo -en "\n"
    usage
    echo -en "\n"
    exit 1
fi

if [ -d $OUTPUT_DIR ]
then
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

for SAMPLE_FOLDER in *
do
    cd ${SAMPLE_FOLDER}/output_IMP/Analysis/annotation
    for GENE_NAME in gyrA gyrB gox phnA phnC phnD phnE phnF phnG phnH phnI phnJ phnK phnL phnM phnN phnX phnW phnR \
      phnS phnT phnU phnP pspE pstC purU purB pphA phoR phoB phoP soxA soxB thiO recA rpoC ktrB ydiF nuoF mntB araD \
      tfdB artI arfA malF qedA mlhB uvrB fdhA fdm recG
    do
        GREP_GENE_NAME="${GENE_NAME}_"
        BED_FILENAME=${BED_DIR}/intersect_${GENE_NAME}_${SAMPLE_FOLDER}.bed
        # get all lines which carry a phn gene and keep only the columns with contig name, start, stop and gene name
        if [ $GENE_NAME == "pphA" ]
        then
            grep $GREP_GENE_NAME annotation.filt.gff | grep "3.11.1.3" > ${TMP_DIR}/tmp.gff
        elif [ $GENE_NAME == "soxA" ] || [ $GENE_NAME == "soxB" ]
        then
            grep $GREP_GENE_NAME annotation.filt.gff | grep "1.5.3.1" > ${TMP_DIR}/tmp.gff
        else
            grep $GREP_GENE_NAME annotation.filt.gff > ${TMP_DIR}/tmp.gff
        fi

        cat ${TMP_DIR}/tmp.gff | \
          awk -F '\t|;' '{for(i=9;i<=NF;i++){if($i~/^gene=/){column=$i}} print $1,$4,$5,$5-$4,column}' > ${BED_FILENAME}

        $SAMTOOLS_BIN view -L ${BED_FILENAME} ../../Assembly/mg.reads.sorted.bam | grep -v -P "^\@" | cut -f 1,3 | \
          sort | uniq | cut -f 2  | sort | uniq -c | perl -ane '$_=~/^\s+(\d+) (.+)$/;chomp($2); print "$2\t$1\n"; ' \
          > ${RESULTS_DIR}/mg.reads.per.gene_${GENE_NAME}_${SAMPLE_FOLDER}.tsv &

        rm ${TMP_DIR}/tmp.gff
    done
    cd ../../../..
done

rmdir -p ${OUTPUT_DIR}/tmp
