#!/bin/bash

# this script removes the -g specified primers from all read libraries in a dir
# it assumes that forward and reverse reads are separate
# -j is the amount of processors
# output is prepended by "cut_"

for library in *R*.gz
do
  cutadapt -j 20 -g CCTACGGGNGGCWGCAG -g GACTACHVGGGTATCTAATCC -o /data/projects/glyphosate/reads/reads_16S_cutadapt/water_dna/cut_${library} $library
done
