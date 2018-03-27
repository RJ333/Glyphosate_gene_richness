#!/bin/sh

# this scripts adds all single read_per_gene files to an big list.

for nfiles in named_*.tsv
do	
  echo $nfiles
  cat $nfiles >> appended_genes.tsv
done
