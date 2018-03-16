#!/bin/sh

# adds the name of the file as first column and deletes the original file

for files in *.tsv
do	
  echo $files
  #awk '{print FILENAME"\t"$0}' $files > $files.bk; mv $files.bk $files
  awk '{a=FILENAME;}{print a"\t"$0}' $files > named_${files}
  rm $files
done


