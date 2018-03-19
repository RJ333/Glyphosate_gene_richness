#!/bin/bash

# this script contains all genes we generated "reads per gene"-values for. It extracts gene and product entry from the prokka annotation
# which can then be used as a dictionary from gene name to function

grep -E 'gene=araD|gene=arfA|gene=artI|gene=bcsZ|gene=celB|gene=cenC|gene=chiA|gene=czcD|gene=fdhA|gene=fdm|gene=gyrA|gene=gyrB|'\
'gene=katE|gene=katG|gene=ktrB|gene=malF|gene=merA|gene=mlhB|gene=mntB|gene=nirS|norB|gene=nosZ|gene=nuoF|gene=phnA|gene=phnC|gene=phnD|'\
'gene=phnE|gene=phnF|gene=phnG|gene=phnH|gene=phnI|gene=phnJ|gene=phnK|gene=phnL|gene=phnM|gene=phnN|gene=phnP|gene=phnR|'\
'gene=phnS|gene=phnT|gene=phnU|gene=phnW|gene=phnX|gene=phoP|gene=phoR|gene=pphA|gene=pspE|gene=pstC|gene=purB|gene=purU|gene=qedA|'\
'gene=recG|gene=rpoC|gene=soxA|gene=soxB|gene=soxC|gene=soxR|gene=soxS|gene=tfdB|gene=thiO|gene=tmoS|gene=uvrB|gene=xynB|gene=ydiF' annotation.filt.gff |\
awk 'BEGIN{OFS = "\t"};
  {match($0,/gene=[^;]*/);
  gene_value=substr($0,RSTART,RLENGTH);
  split(gene_value,a,"_");
  split(a[1],b,"=");
  match($0,/product=.*/);
  print b[2],substr($0,RSTART,RLENGTH)}' |\ 
  tr " " "-" |\ 
  awk 'BEGIN{ OFS = "\t"}{split($2,a,"="); print $1, a[2]}' |\
  awk -F"\t" '!_[$1]++' | sort | uniq
  