#!/bin/bash
for xy in *R*.gz
do
  cutadapt -j 20 -g CCTACGGGNGGCWGCAG -g GACTACHVGGGTATCTAATCC -o /data/projects/glyphosate/reads/reads_16S_cutadapt/water_dna/cut_${xy} $xy
done