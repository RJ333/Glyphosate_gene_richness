#!/bin/bash
# this script assumes that f and r reads are still separated
# and primer set 341f-805r (V3-V4)

# define input and output dir
input=$1
output=$2

# run cutadapt with fixed primer sequences
for xy in $1/*R*.gz
do
  cutadapt -j 20 -g CCTACGGGNGGCWGCAG -g GACTACHVGGGTATCTAATCC -o $2/cut_${xy} ${xy} 
done