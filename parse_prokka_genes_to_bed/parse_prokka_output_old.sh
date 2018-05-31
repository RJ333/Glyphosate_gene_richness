#!/usr/bin/env bash

BASE_DIR=/data/IOW/BIO/Environmental_Microbiology/Glyphosate/schweden17/prokka

usage()
{
cat << EOF
usage:
$0 -i <INPUT_DIR> -o <OUTPUT_DIR> [-g <GENE_NAME>] [-b <BED_FILE>]

The script creates a bed file based on prokka annotations.

OPTIONS:
   -i          Input directory (raw fastq files, mandatory)
   -o          Output directory (mandatory)
   -g          Gene name to parse for (default: phnM)
   -b          Output bed file

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

SCRIPT_FOLDER=/data/projects/scripts/Glyphosate_gene_richness/parse_prokka_genes_to_bed
PARSE_PHN_TMP_SCRIPT=${SCRIPT_FOLDER}/parse_annotation.py
PHN_TMP_FILE=${OUTPUT_DIR}/prokka_annotation.txt
CONTIG_LENGTH_FILE=${OUTPUT_DIR}/contig_lengths.txt
BED_FILE=${OUTPUT_DIR}/${BED_FILENAME}

# get all lines which carry a phn gene and keep only the columns with contig name, start, stop and gene name
find ${INPUT_DIR} -regextype sed -regex ".*/PROKKA[0-9].gff" -exec cat {} \; | grep "Prodigal" | cut -f 1,4,5,9 \
  > ${PHN_TMP_FILE}

# create contig_length.txt
cat ${INPUT_DIR}/prokka*/PROKKA*.gff | grep "^##" > ${CONTIG_LENGTH_FILE}

python $PARSE_PHN_TMP_SCRIPT --input_file $PHN_TMP_FILE --contig_length $CONTIG_LENGTH_FILE --gene_name $GENE_NAME --output_file $BED_FILE
