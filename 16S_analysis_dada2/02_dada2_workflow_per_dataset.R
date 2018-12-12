#!/usr/bin/env R

#' This is the general dada2 workflow for primer removed reads. It is intended to be run separately for a characteristic
#' dataset, divided e.g. by sequencing runs or sampling types. The workflow needs to be performed interactively
#' as we supervise the trimming process visually (using ggplot2). Additionally, we check manually the length  
#' distribution of the combined forward and reverse read pairs which gives us information on how well the trimming
#' and the chimera removal performed. We will finally save the output as phyloseq object, which can easily be merged
#' with more phyloseq objects, e.g. from parallel runs of this workflow for other data sets from the same experiment.

# define working directory
setwd("/data/projects/glyphosate/reads/")
#load(file = "dada2_water_dna.RData")

# Load packages into session, and print package version
.cran_packages <- c("ggplot2", 
					"gridExtra")
.bioc_packages <- c("dada2", 
					"phyloseq", 
					"DECIPHER", 
					"phangorn")
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# define paths
cutadapt_reads_path <- file.path("../reads_16S_cutadapt", "water_dna")
dada2_trimmed_path <- file.path("./dada2_processed", "water_dna")

set.seed(100)


# collect reads
amplicon_libraries <- sort(list.files(cutadapt_reads_path, full.names = TRUE))
raw_forward_libraries <- amplicon_libraries[grepl("R1", amplicon_libraries)]
raw_reverse_libraries <- amplicon_libraries[grepl("R2", amplicon_libraries)]

# plot sequence quality of 3 random samples before trimming
# this doesn't work always, probably due to empty fastq lines
# you can still select manually a couple of samples
random_sample_picks <- sample(length(raw_forward_libraries), 3)
for(random_sample in random_sample_picks) {print(plotQualityProfile(raw_forward_libraries[i]) + ggtitle("Fwd"))}
for(random_sample in random_sample_picks) {print(plotQualityProfile(raw_reverse_libraries[i]) + ggtitle("Rev"))}

if(!file_test("-d", dada2_trimmed_path)) dir.create(dada2_trimmed_path)
filtered_forward_libraries <- file.path(dada2_trimmed_path, basename(raw_forward_libraries))
filtered_reverse_libraries <- file.path(dada2_trimmed_path, basename(raw_reverse_libraries))

# Trim and Filter
# if we trim too much, the pairing is based on too few overlapping bases and very unspecific
# this may result in many merged_reads with very different length
# what we will see below when we check the length distribution 

# the values used below are based on testing and the plotQualityProfile function
filter_result <- filterAndTrim(raw_forward_libraries, filtered_forward_libraries, 
							   raw_reverse_libraries, filtered_reverse_libraries, 
							   truncLen = c(280, 225),
							   maxN = 0, 
							   maxEE = c(2.5, 5), 
							   truncQ = 2, 
							   rm.phix = FALSE,
							   trimLeft = c(10, 0), 
							   minLen = c(270, 220), 
							   compress = TRUE, 
							   multithread = TRUE)

# reads after trimming in absolute and relative values
head(filter_result)
filter_result[,2]/filter_result[,1]*100

# now check the quality after trimming
for(random_sample in random_sample_picks) {print(plotQualityProfile(filtered_forward_libraries[i]) + ggtitle("Fwd filt"))} 
for(random_sample in random_sample_picks) {print(plotQualityProfile(filtered_reverse_libraries[i]) + ggtitle("Rev filt"))}

# dereplicate reads, keep abundance and quality information
dereplicated_forward <- derepFastq(filtered_forward_libraries)
dereplicated_reverse <- derepFastq(filtered_reverse_libraries)

# generate sample name from full description
# adjust split character according to your samples
sample_names <- sapply(strsplit(basename(filtered_forward_libraries), "_S"), `[`, 1)
names(dereplicated_forward) <- sample_names
names(dereplicated_reverse) <- sample_names

# Infer sequence variants with training, subset should be about 100 M bases
# the first 60 samples should be enough
error_rates_forward <- dada(dereplicated_forward[1:60], 
							err = NULL, 
							selfConsist = TRUE, 
							multithread = TRUE) 
error_rates_reverse <- dada(dereplicated_reverse[1:60], 
							err = NULL, 
							selfConsist = TRUE, 
							multithread = TRUE) 

# check the sequencing error distribution
plotErrors(error_rates_forward, nominalQ = TRUE)
plotErrors(error_rates_reverse, nominalQ = TRUE)	

# correct seq errors by applying learned error rates on all seqs
corrected_reads_forward <- dada(dereplicated_forward, 
								err = error_rates_forward[[1]]$err_out, 
								pool = TRUE, 
								multithread = TRUE) 
corrected_reads_reverse <- dada(dereplicated_reverse, 
								err = error_rates_reverse[[1]]$err_out, 
								pool = TRUE, 
								multithread = TRUE)

# Join corrected forward and reverse reads, 
merged_reads <- mergePairs(corrected_reads_forward, dereplicated_forward, corrected_reads_reverse, dereplicated_reverse)

# Construct sequence table and remove chimeras
seqtab.all <- makeSequenceTable(merged_reads)

# check number of samples and uniq seqs
dim(seqtab.all)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.all)))
hist(nchar(getSequences(seqtab.all)), 
	 main = "Distribution of sequence lengths", 
	 breaks = (seq(250, 
                   500, 
                   by = 10)))

# keep expected and/or abundant sequence lengths and check distribution again
seqtab2 <- seqtab.all[,nchar(colnames(seqtab.all)) %in% seq(392, 448)]

table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)), 
	 main = "Distribution of picked sequence lengths", 
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

save.image(file = "dada2_water_dna.RData")