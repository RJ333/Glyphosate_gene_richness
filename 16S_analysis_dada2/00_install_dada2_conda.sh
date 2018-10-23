# install Bioconductor and dada2 in conda dada2
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n dada2 bioconductor-shortread=1.38.0 r-base=3.5.1 bioconductor-dada2=1.8

# if conda is not available, e.g. on io-49
source ~/.bashrc

conda activate dada2
R