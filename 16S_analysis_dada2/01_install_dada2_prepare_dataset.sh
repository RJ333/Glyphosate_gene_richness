# install Bioconductor and dada2 in conda dada2
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n dada2 bioconductor-shortread=1.38.0 r-base=3.5.1 bioconductor-dada2=1.8
# https://benjjneb.github.io/dada2/bigdata.html
# https://benjjneb.github.io/dada2/tutorial.html

# Samples have been demultiplexed, i.e. split into individual per-sample fastq files.
# Non-biological nucleotides have been removed, e.g. primers, adapters, linkers, etc.
# If paired-end sequencing data, the forward and reverse fastq files contain reads in matched order

# following tutorial on https://f1000research.com/articles/5-1492/v2


# first installing these packages manuallysource
source("http://bioconductor.org/biocLite.R")
biocLite(c("knitr", "BiocStyle"))

library("knitr")
library("BiocStyle")

.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
   install.packages(.cran_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
   source("http://bioconductor.org/biocLite.R")
   biocLite(.bioc_packages[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

set.seed(100)

miseq_path <- file.path("/data/projects/dada2_tutorial", "MiSeq_SOP")
filt_path <- file.path("/data/projects/dada2_tutorial", "filtered")

if(!file_test("-d", miseq_path)) {
  dir.create(miseq_path)
  download.file("http://www.mothur.org/MiSeqDevelopmentData/StabilityNoMetaG.tar",
                 destfile = file.path(miseq_path, "StabilityNoMetaG.tar"))
  system(paste0("tar -xvf ", file.path(miseq_path, "StabilityNoMetaG.tar"),
                 " -C ", miseq_path, "/"))
}

fns <- sort(list.files(miseq_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }

if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
		      c(filtFs[[i]], filtRs[[i]]),
                      trimLeft=10, truncLen=c(245, 160),
                      maxN=0, maxEE=2, truncQ=2,
                      compress=TRUE)
}

derepFs <- derepFastq(filtsFs)
derepRs <- derepFastq(filtsRs)
sam.names <- sapply(strsplit(basename(filtsFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names