# download silva training set from 
# https://www.dropbox.com/sh/mfcivbudmc21cqt/AAB1l-AUM5uKvjrR33ct-cTXa?dl=0

# using -P, this somehow was downloaded to dada2_tutorial/, but not into db?
wget -O /data/projects/dada2_tutorial/db/silva_nr_v123_train_set.fa.gz\
  'https://www.dropbox.com/sh/mfcivbudmc21cqt/AAAPdfcBQgLSS5JmnjCDMuWfa/silva_nr_v123_train_set.fa.gz?dl=0'
cd /data/projects/dada2_tutorial/
mv silva_nr_v123_train_set.fa.gz ./db

wget -O /data/projects/dada2_tutorial/db/silva_species_assignment_v123.fa.gz\ 
  'https://www.dropbox.com/sh/mfcivbudmc21cqt/AAAvMTycF0dhRWJJt_CJLcbia/silva_species_assignment_v123.fa.gz?dl=0'

# install Bioconductor and dada2 in conda dada2
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n dada2 bioconductor-shortread=1.38.0 r-base=3.5.1 bioconductor-dada2=1.8

# if conda is not available, e.g. on io-49
source ~/.bashrc

conda activate dada2
R