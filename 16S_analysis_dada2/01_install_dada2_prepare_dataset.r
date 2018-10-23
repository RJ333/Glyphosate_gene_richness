# https://benjjneb.github.io/dada2/bigdata.html
# https://benjjneb.github.io/dada2/tutorial.html

# Samples have been demultiplexed, i.e. split into individual per-sample fastq files.
# Non-biological nucleotides have been removed, e.g. primers, adapters, linkers, etc.
# If paired-end sequencing data, the forward and reverse fastq files contain reads in matched order

# following tutorial on https://f1000research.com/articles/5-1492/v2

# define working directory to story RData image
setwd("/data/projects/dada2_tutorial")
load(file = "dada2_tutorial.RData")



# first installing these packages manually
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

# be careful, the original script contains filtFs and filt"s"Fs!!

if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
		      c(filtFs[[i]], filtRs[[i]]),
                      trimLeft = 10, truncLen = c(245, 160),
                      maxN = 0, maxEE = 2, truncQ = 2,
                      compress = TRUE)
}

derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names

ddF <- dada(derepFs[1:40], err = NULL, selfConsist = TRUE)
# Convergence after  8  rounds.  > 6 hours, bio-48
ddR <- dada(derepRs[1:40], err = NULL, selfConsist = TRUE)

dadaFs <- dada(derepFs, err = ddF[[1]]$err_out, pool = TRUE, multithread = TRUE)
dadaRs <- dada(derepRs, err = ddR[[1]]$err_out, pool = TRUE, multithread = TRUE)

save.image(file = "dada2_tutorial.RData")