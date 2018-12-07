# in this script I tested how to combine objects from different workspaces
# namely the processed dada2 objects, which were treated separately 


# define working directory to story RData image
setwd("/data/projects/glyphosate/reads/dada2_processed/water_dna")

# load first workspace representing first amplicon sequencing run
load(file = "dada2_water_dna.RData")

# this will be the workspace for the combined dada2 objects
save.image(file = "../dada2_water_comb.RData")

# define required packages
.cran_packages <- c("ggplot2", 
					"gridExtra")
.bioc_packages <- c("dada2", 
					"phyloseq", 
					"DECIPHER", 
					"phangorn")
					
# check if packages are installed
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

# rename respective object from workspace to differentiate from next workspace
dna_water_nochim <- seqtab2.nochim

# add other work spaces to temporary environment
tmp.env <- new.env()
load("../water_cdna/dada2_water_cdna.RData", envir = tmp.env)

# get the objects you need into your globalenv and remove tmp.env
cdna_water_nochim <- get("seqtab2.nochim", pos = tmp.env)
rm(tmp.env)

# merge the sequencetable from two different workspaces into one table
water_seqtable <- mergeSequenceTables(table1 = dna_water_nochim, 
									  table2 = cdna_water_nochim,
									  repeats = "error", 
									  orderBy = "abundance")
  
# assign taxonomy with specially prepared silva database
# NOTE: it should be possible to use the mothur-like amplicon-trimmed database
# to reduce memory amount. Not tested

ref_fasta <- "/data/db/silva_nr_v132_train_set.fa.gz"
taxtab <- assignTaxonomy(water_seqtable, 
						 refFasta = ref_fasta, 
						 multithread = TRUE)
colnames(taxtab) <- c("Kingdom", 
					  "Phylum", 
					  "Class", 
					  "Order", 
					  "Family", 
					  "Genus")
seqs <- getSequences(water_seqtable)

# turn dada2-output into phyloseq object
# the sample data can be added later 
# the tree can later be calculated based on the ps object, after filtering the
# low abundant OTUs
ps_dada <- phyloseq(tax_table(taxtab), 
					otu_table(water_seqtable, 
							  taxa_are_rows = FALSE), 
					#phy_tree(fitGTR$tree),
					refseq(seqs))
					
### dada2 uses actual sequences as OTU headers

# this changes the header from the actual sequence to Seq_001, Seq_002 etc
taxa_names(ps_new)
n_seqs <- seq(ntaxa(ps_dada))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps_dada) <- paste("Seq", 
							formatC(n_seqs, 
									width = len_n_seqs, 
									flag = "0"), 
							sep = "_")
taxa_names(ps_dada)