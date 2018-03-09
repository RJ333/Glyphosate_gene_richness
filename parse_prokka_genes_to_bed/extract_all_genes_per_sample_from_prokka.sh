#!bin/sh
# extract_all_per_sample_from_prokka

# this script looks for the column that contains "gene" and takes it and the previous one, which is either "eC_number=..." or "ID=PROKKA..."
# the file gets cleaned up so that the numbers as phnM_1 are removed, the "gene=" is removed, all ID=PROKKA is removed. It counts the
# number of the same genes found per assembly with E C number, if present

# for soxA we have e.g.
# 4 soxA 1.5.3.1
# 8 soxA 1.8.2.-

cd /drives/d/data/samples/

for assembly in $(find . -mindepth 1 -maxdepth 1 -type d) 
do 
  cd $assembly
  echo "$assembly is being processed"
  awk -F '\t|;' '{col = 0; for (i=9; i<=NF; i++) {if($i~/^gene=/) {col = i}; } if (col>0) {print $col,$(col-1)}}' annotation.filt.gff | \
    awk '{if($2~/ID=/) {print $1} else {print $0}}' | \
	awk '{sub(/_.*/,"",$1);print}' | \
	awk '{sub(/.*=/,"",$1);print}' | \
	awk '{sub(/.*=/,"",$2);print}' | \
	sort | uniq -c > ${assembly}_gene_list 
  echo "gene list for assembly ${assembly} has been generated, leaving folder $assembly"
  cd ..
done
