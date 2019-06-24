#!/usr/bin/env Rscript
# load libraries
library(phyloseq)
library(Biostrings)
library(dplyr)
library(tidyr)
library(DECIPHER)
library(phangorn)

# Set the working dir
setwd("/data/projects/glyphosate/reads/mothur_processed/")

# set path for plots, create directory if not existant
plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

# mothur files and additional files that need to be imported from the working dir
# shared file is the OTU count table
# constaxonomy contains the taxonomy of the OTUs
our_shared_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.shared"
our_cons_taxonomy_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.0.02.cons.taxonomy"

# sample data is in the metafile
metafile <- read.delim("/data/projects/glyphosate/analysis/metadata/metafile.tsv", row.names = 1, header = TRUE, na.strings = "")
metafile <- sample_data(metafile)

# read OTU representative sequences
OTU_seqs <- readDNAStringSet(file = "OTU_reps_fasta_002.fasta",
  format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

# for testing purposes
# on command line in working dir: 
# tail OTU_reps_fasta_002.fasta | sed 's/Otu/Ref/g' > example_ref.fasta
reference_seqs <- readDNAStringSet(file = "example_ref.fasta",
  format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)
rs <- reference_seqs

# import mothur output into phyloseq
mothur_phyloseq_object <- import_mothur(mothur_list_file = NULL,
  mothur_group_file = NULL, mothur_tree_file = NULL, cutoff = NULL,
  mothur_shared_file = our_shared_file,
  mothur_constaxonomy_file = our_cons_taxonomy_file,
  parseFunction = parse_taxonomy_default)
                                    
# add meta data and OTU representative seqs to phyloseq_object, reduce for testing purposes: 10 otus, 10 samples
mothur_full <- merge_phyloseq(mothur_phyloseq_object, metafile, refseq(OTU_seqs))
OTU <- otu_table(mothur_full)[1:10]
OTU <- OTU[, c(1:10)]
otu_table(mothur_full) <- OTU
ps0 <- mothur_full

# create an otu table out filled with 1 based on a 
# reference sample with the number of reference seqs
otumat <- matrix(1, nrow = 1, ncol = length(rs))
colnames(otumat) <- names(rs)
rownames(otumat) <- "Reference"
OTU <- otu_table(otumat, taxa_are_rows = FALSE)

# use the original sample data to create additional entry with NAs for reference samples
SAM <- sample_data(ps0)[1,]
SAM[,] <- NA
sample_names(SAM) <- "Reference"

# generate taxonomy, split strings into columns
tax_df <- data.frame(wholetax = c("level1.level2.level3.level4.level5.level6",
  "peter.paul.mary.joe.chandler.rachel",
  "one.two.three.four.five.six",
  "and.another.line.and.another.line",
  "now.it.is.enough.with.examples"))
rownames(tax_df) <- names(rs)  
# split strings into tax level columns, headers same as in to-be-merged-with ps object
tax_df <- tax_df %>% separate(wholetax, colnames(tax_table(ps0)))
ref_taxa <- tax_table(as.matrix(tax_df))
ps.ref <- phyloseq(OTU, SAM, ref_taxa, rs)

# now merge own data and external reference data
ps.merged <- merge_phyloseq(ps0, ps.ref)

# calculate tree including reference seqs
seqs <- refseq(ps.merged)
alignment <- AlignSeqs(seqs, anchor = NA)
phang.align <- as.phyDat(alignment, type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data = phang.align)
fitGTR <- update(fit, k = 4, inv = 0.2)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                       rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload = TRUE)
# merge ps object with tree and plot tree including reference sequences
ps.merged <- merge_phyloseq(ps.merged, phy_tree(fitGTR$tree))
plot_tree(ps.merged)
