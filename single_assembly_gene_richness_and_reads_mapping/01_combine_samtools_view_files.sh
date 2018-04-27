#!/bin/sh

# adds the name of the file in folder "unnamed" as first column and copies to "named" folder 
cd /data/Rene/glyph/renamed
for files in *.tsv
do	
  echo $files
  #awk '{print FILENAME"\t"$0}' $files > $files.bk; mv $files.bk $files
  awk '{a=FILENAME;}{print a"\t"$0}' $files > ../named_reads_per_gene/named_${files}
done


