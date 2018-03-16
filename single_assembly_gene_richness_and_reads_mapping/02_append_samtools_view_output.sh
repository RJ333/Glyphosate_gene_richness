#!/bin/sh

for nfiles in named_*.tsv
do	
  echo $nfiles
  cat $nfiles >> appended_genes.tsv
done
