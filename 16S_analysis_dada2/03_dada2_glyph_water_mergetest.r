#!/usr/bin/env R

#' In this script we combine the sequence tables from the separately performed dada2 sequence processing runs.
#' dada2 corrected the sequences for sequencing errors and therefore treats each unique reference as own OTU or ASV
#' (Amplicon Sequence Variant). The sequence tables still contain the actual amplicon sequence as taxa. This allows 
#' us to identify same ASVs from different runs (which otherwise could be called OTU3 in run 1 and OTU7 in run 2).
#' The data was saved as RDS object. We will read it, merge it and assign taxonomy to it. Finally, all data will be stored
#' within a phyloseq object.

####### TO DO: so far only tested for the merging step..., go the full way

# define working directory to story RData image
setwd("/data/projects/glyphosate/reads/dada2_processed/")

# define required packages
.cran_packages <- c("ggplot2", 
					"gridExtra")
.bioc_packages <- c("dada2", 
					"phyloseq", 
					"DECIPHER", 
					"phangorn")

# set path to data base for taxonomic annotation
ref_fasta <- "/data/db/silva_nr_v132_train_set.fa.gz"					
					
# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# load the sequence tables
seqtab_dna_water <- readRDS("./water_dna/seqtab_dna_water.RDS")
seqtab_cdna_water <- readRDS("./water_cdna/seqtab_cdna_water.RDS")
seqtab_biofilm <- readRDS("./biofilm/seqtab_biofilm.RDS")

# merge the sequencetable from two different runs into one table
all_seqtables <- mergeSequenceTables(table1 = seqtab_dna_water, 
									  table2 = seqtab_cdna_water,
									  table3 = seqtab_biofilm,
									  repeats = "error", 
									  orderBy = "abundance")
  
# assign taxonomy with specially prepared silva database
taxtab <- assignTaxonomy(all_seqtables, 
						 refFasta = ref_fasta, 
						 multithread = TRUE)
# get unique sequences and name them also with their sequence
seqs <- getSequences(all_seqtables)
names(seqs) <- seqs
seqs <- DNAStringSet(seqs)

# we use a phyloseq object here to conveniently change the headers of all interlinked data
# combine sequences, counts and taxonomy into phyloseq object
ps_dada <- phyloseq(tax_table(taxtab),
					otu_table(all_seqtables,
							  taxa_are_rows = FALSE),
					#phy_tree(fitGTR$tree),
					refseq(seqs))

# this changes the header from the actual sequence to ASV_001, ASV_002 etc
taxa_names(ps_dada)
n_seqs <- seq(ntaxa(ps_dada))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps_dada) <- paste("ASV", 
							formatC(n_seqs, 
									width = len_n_seqs, 
									flag = "0"), 
							sep = "_")
taxa_names(ps_dada)

# write otu table, taxonomy and seqs out
write.table(as.data.frame(otu_table(ps_dada)), file = "dada2_otu.tsv", sep = "\t")
write.table(as.data.frame(tax_table(ps_dada)), file = "dada2_tax.tsv", sep = "\t")
writeXStringSet(refseq(ps_dada), "ASV_seqs.fasta", format = "fasta")