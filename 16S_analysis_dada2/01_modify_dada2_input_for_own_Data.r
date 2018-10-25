# adjust dada2 for own samples

# check meta data csv

/mnt/d/data/mothur_sample

# define working directory to story RData image
# setwd("/mnt/d/data/mothur_sample")
setwd("D:/data/mothur_sample")

wget -O /mnt/d/data/mothur_sample/metadata/metadata.csv\
  'https://www.dropbox.com/s/wpxnvsui4mjj8z3/MIMARKS_Data_combined.csv?raw=1'

### the outcommented steps are only required the first time

# first installing these packages manually
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("knitr", "BiocStyle"))

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

# here starts the processing
fns <- sort(list.files("D:/data/mothur_sample", full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

# Trim and Filter

# takes 3 random samples and plots their quality value over length
ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }

# generates file names for filtered files
if(!file_test("-d", "D:/data/mothur_sample/filtered")) dir.create("D:/data/mothur_sample/filtered")
filtFs <- file.path("D:/data/mothur_sample/filtered", basename(fnFs))			# be careful, the original script 
filtRs <- file.path("D:/data/mothur_sample/filtered", basename(fnRs))			# contains filtFs and filt"s"Fs!!
# trimming reads
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
		      c(filtFs[[i]], filtRs[[i]]),
                      trimLeft = 10, truncLen = c(245, 160),
                      maxN = 0, maxEE = 2, truncQ = 2,
                      compress = TRUE)
}

# Infer sequence variants
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names


ddF <- dada(derepFs[1:40], err = NULL, 
			selfConsist = TRUE) # Convergence after 8 rounds.> 6 hours, bio-48???? 
								# few minutes on laptop
ddR <- dada(derepRs[1:40], err = NULL, 
			selfConsist = TRUE) # Convergence after 7 rounds.> 1 hour, bio-49
								# few minutes on Laptop

plotErrors(ddF) # this doesn't work yet due to failing x11 forwarding (ubuntu subsystem, bio 49?)
plotErrors(ddR)	# but x11 forwarding works with mobaxterm	

test <- dada(derepFs, err = ddF[[1]]$err_out, pool = TRUE) # takes about 2 - 3 hours without multithread
dadaRs <- dada(derepRs, err = ddR[[1]]$err_out, pool = TRUE, multithread = TRUE) # 20 mins for multiple cores
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)