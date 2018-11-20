# this is the dada2 workflow for the biofilm sequencing run with specific
# trimming parameters

# define working directory to story RData image
setwd("/data/projects/glyphosate/reads/dada2_processed/biofilm")
load(file = "dada2_biofilm.RData")

save.image(file = "dada2_biofilm.RData")

# setup bioconductor software including dada2 and phyloseq
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

miseq_path <- file.path("/data/projects/glyphosate/reads/reads_16S_cutadapt", 
						"biofilm")
filt_path <- file.path("/data/projects/glyphosate/reads/dada2_processed", 
						"biofilm")

# here starts the processing
fns <- sort(list.files(miseq_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

# plot sequence quality of 3 random samples before trimming
# this doesn't work always, probably due to empty fastq lines
# you can still select manually a couple of samples
ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }

# if we trim too much, the pairing is based on few bases and very unspecific
# this can result in many paired reads with very different length
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

# Trim and Filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(280, 235),
                  maxN = 0, maxEE = c(2.5, 5), truncQ = 2, rm.phix = FALSE,
                  trimLeft = c(10, 0), minLen = c(270, 235), 
				  compress = TRUE, multithread = TRUE)
				  
# reads after trimming
head(out)

# as percentage 
out[,2]/out[,1]*100

# now plot the quality after trimming
ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(filtFs[i]) + ggtitle("Fwd filt")) }
for(i in ii) { print(plotQualityProfile(filtRs[i]) + ggtitle("Rev filt")) }

# dereplicate reads, keep abundance and quality information
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
# split character adjusted
sam.names <- sapply(strsplit(basename(filtFs), "_S"), `[`, 1) # split character adjusted
names(derepFs) <- sam.names
names(derepRs) <- sam.names

# Infer sequence variants with training subset should be about 100 M bases
ddF <- dada(derepFs[1:60], err = NULL, selfConsist = TRUE, multithread = TRUE) 
ddR <- dada(derepRs[1:60], err = NULL, selfConsist = TRUE, multithread = TRUE) 

# check the sequencing error distribution
plotErrors(ddF, nominalQ = TRUE)
plotErrors(ddR, nominalQ = TRUE)	

dadaFs <- dada(derepFs, err = ddF[[1]]$err_out, pool = TRUE, multithread = TRUE) 
dadaRs <- dada(derepRs, err = ddR[[1]]$err_out, pool = TRUE, multithread = TRUE)
# Join forward and reverse reads, corrected for seq errors
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# Construct sequence table and remove chimeras
seqtab.all <- makeSequenceTable(mergers)
# number of samples and uniq seqs
dim(seqtab.all)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.all)))
hist(nchar(getSequences(seqtab.all)), 
	 main = "Distribution of sequence lengths", 
	 breaks = (seq(250, 
				   500, 
				   by = 10)))

# pick expected and abundant sequence lengths and check again
seqtab2 <- seqtab.all[,nchar(colnames(seqtab.all)) %in% seq(392, 448)]

table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)), 
	 main = "Distribution of sequence lengths", 
	 breaks = (seq(390, 
				   450, 
				   by = 5)))

# remove chimera based on denoised sequences
# checks if both parts of a contig represents either part of more abundant seqs
seqtab2.nochim <- removeBimeraDenovo(seqtab2, 
									 method = "consensus", 
									 multithread = TRUE, 
									 verbose = TRUE)

# number of samples and uniq seqs
dim(seqtab2.nochim)

# ratio of non-chimeric to total seqs
# chimeras can make up many but low abundant seqs
sum(seqtab2.nochim)/sum(seqtab2)

# include deseq2 here?
# we could probably combine the different sequencing runs here before assigning taxonomy

# Assign taxonomy (I used silva instead of rdp classifier)
ref_fasta <- "/data/db/silva_nr_v132_train_set.fa.gz"
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)