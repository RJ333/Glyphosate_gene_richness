#!/bin/sh

: '
This script is a wrapper around cutadapt to remove oligonucleotides 
e.g. primers. It requires that forward and reverse reads are separate. 
The output is prepended by "cut_"
'

: '
Args:
    -i|--input 			: the absolute input path, containing gzipped reads
	-o|--output 		: the absolute output path
	-c|--cores 			: the number of cores to use
	-p|--primersequence : the oligonucleotide sequence to be removed from the beginning of the read!
	-d|--direction 		: select between forward reads "R1" and reverse reads "R2"

Return:
	the selected read libraries without the specified sequence, gzipped.
'


POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -i|--input)
    INPUT="$2"
    shift # past argument
    shift # past value
    ;;
	-o|--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
	-c|--cores)
    CORES="$2"
    shift # past argument
    shift # past value
    ;;
	-p|--primersequence)
    PRIMER_SEQUENCE="$2"
    shift # past argument
    shift # past value
    ;;
	-d|--direction)
    DIRECTION="$2"
    shift # past argument
    shift # past value
    ;;
	-*) 
	echo "Unknown argument: \"$key\""; 
	exit 1
	;;
	*)
	break
	;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

echo INPUT = "${INPUT}"
echo OUTPUT = "${OUTPUT}"
echo CORES = "${CORES}"
echo PRIMER_SEQUENCE = "${PRIMER_SEQUENCE}"
echo DIRECTION = "${DIRECTION}"

# run cutadapt
cd ${INPUT}
for library in *${DIRECTION}*.gz
do
  cutadapt -j $CORES -g $PRIMER_SEQUENCE -o ${OUTPUT}/cut_${library} $library 
done
