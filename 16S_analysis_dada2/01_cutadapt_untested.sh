#!/bin/sh

: 'multiline comment in bash
This script is a wrapper around cutadapt to remove oligos such as primers.
It assumes that forward and reverse reads are separate. The output is 
prepended by "cut_"
'

: '
Args:
    -i|--input 			: the input path, containing gzipped reads
	-o|--output 		: the output path
	-c|--cores 			: the number of cores
	-p|--primersequence : the oligonucleotide sequence to be removed from the beginning of the read!
	-d|--direction 		: choose the reads to work on: "R1" for forward, "R2" is reverse, "R" for both directions
       --default 		: I dont know what it does, but I will see whether it is required

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
	--default)
    DEFAULT=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
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
  echo $library
  echo $CORES
  echo ${OUTPUT}/cut_${library} 
done

# return unused arguments (I guess?)
# if [[ -n $1 ]]; then
    # echo "Last line of file specified as non-opt/last argument:"
    # tail -1 "$1"
# fi
