#!/bin/bash
cd /data/Rene/glyph/unnamed_reads_per_gene
# Loop over all files in the current directory

# would be better to first select those filenames who actually need renaming
# e.g. by searching for more than 2 "_" with grep?
# find . -mindepth 1 -maxdepth 1 -type f | grep -E "_{3,}"
for i in *; do

  # Extract the part before the first _
  head="${i%%_*}"

  # Get the rest of the string
  tail=${i#*_}

  # Extract the part after the last _
  rhead="${tail##*_}"

  # Extract the "middle" portion
  rtail="${tail%_*}"

  # Substitute _ with @ in the "middle"
  fixedrtail=${rtail//_/@}

  # Rename files
  echo -e "Renaming \"$i\" to \"$head_${fixedrtail}_$rhead\""
  mv $i ../renamed/"${head}_${fixedrtail}_${rhead}"
done