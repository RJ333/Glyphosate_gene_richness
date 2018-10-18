# install Bioconductor and dada2 in conda Renv
conda create -n Renv351
conda activate Renv351
conda install -c r r

R

setwd("/data/projects/glyphosate/analysis_16S")
setwd("R_351_files")
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("dada2")


configure: error: zlib not found
ERROR: configuration failed for package ‘ShortRead’
# zlib not found, packages ShortRead and dada2 could not be installed

# need sudo
apt-get install zlib1g-dev 

biocLite("ShortRead")
biocLite("dada2")

# https://benjjneb.github.io/dada2/bigdata.html
# https://benjjneb.github.io/dada2/tutorial.html

# Samples have been demultiplexed, i.e. split into individual per-sample fastq files.
# Non-biological nucleotides have been removed, e.g. primers, adapters, linkers, etc.
# If paired-end sequencing data, the forward and reverse fastq files contain reads in matched order