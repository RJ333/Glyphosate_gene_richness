#!/usr/bin/env bash

BASE_DIR=/data_bio/jwerner/glyphosate/A1/output_IMP/Analysis/annotation

usage()
{
cat << EOF
usage:
$0 -i <INPUT_DIR> -o <OUTPUT_DIR> [-g <GENE_NAME>] [-b <BED_FILE>]

The script creates a bed file based on prokka annotations.

OPTIONS:
   -i          Input directory (prokka annotation directory, mandatory)
   -o          Output directory (mandatory)
   -g          Gene name to parse for (default: phnM)
   -b          Output bed file (default: intersect.bed)

   -h/--help   Show this message

EOF
}

GENE_NAME="phnM"
BED_FILENAME=intersect.bed

while getopts ":i:o:g:b:" o; do
    case "${o}" in
        i)
            INPUT_DIR=${OPTARG}
            ;;
        o)
            OUTPUT_DIR=${OPTARG}
            ;;
        g)
            GENE_NAME=${OPTARG}
            ;;
        b)
            BED_FILE=${OPTARG}
            ;;
        *)
            # usage
            ;;
    esac
done
shift $((OPTIND-1))

BED_FILE=${OUTPUT_DIR}/${BED_FILENAME}

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
    mkdir $OUTPUT_DIR
fi

# get all lines which carry a phn gene and keep only the columns with contig name, start, stop and gene name
grep $GENE_NAME ${INPUT_DIR}/annotation.filt.gff | \
  awk -F '\t|;' '{for(i=9;i<=NF;i++){if($i~/^gene=/){column=$i}} print $1,$4,$5,$5-$4,column}'  > ${BED_FILE}
