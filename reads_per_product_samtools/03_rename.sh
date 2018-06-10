#!/bin/bash
BASE_DIR=/data/Rene/glyph/prokka/
mkdir -p $BASE_DIR/renamed

# lets first move those files who are NOT having too many "_". 
# The "+" reduces the number of called shells
# positioning of "-not" is important
cd $BASE_DIR/results
find  . -mindepth 1 -maxdepth 1 -type f -not -name '*_*_*_*' -exec mv -t $BASE_DIR/renamed {} +

# then we treat the file names that need treatment
# Loop over all remaining files in the current directory
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
  mv $i $BASE_DIR/renamed/"${head}_${fixedrtail}_${rhead}"
done
