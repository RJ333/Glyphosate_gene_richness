#!/bin/bash
# this script assumes that f and r reads are still separated
# and primer set 341f-805r (V3-V4)

input=$1 # argument passing has not been tested yet
output=$2

for xy in $1/*R*.gz
do
  cutadapt -j 20 -g CCTACGGGNGGCWGCAG -g GACTACHVGGGTATCTAATCC -o $2/cut_${xy} ${xy} 
done