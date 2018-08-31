#!/bin/bash

cat test1 | while read line
do
  headerline=$(awk '/>/{print $0}' $line)
  fastaline=$(awk '!/>/{print $0}' $line)
  fastaline_length=$(awk -v linelength=$fastaline '{print length(linelength)}')

  echo ${headerline}_${fastaline_length}
  echo $fastaline