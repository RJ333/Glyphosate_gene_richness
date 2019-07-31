#!/usr/bin/env Rscript
library(phyloseq)
library(Biostrings)
library(tidyverse)

setwd("/data/projects/glyphosate/reads/mothur_processed/")

plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

# import files
our_shared_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.shared"
our_cons_taxonomy_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.0.02.cons.taxonomy"

metafile2 <- read.delim("/data/projects/glyphosate/analysis/metadata/metafile2.tsv", 
  row.names = 1, header = TRUE, na.strings = "")
metafile2 <- sample_data(metafile2)

OTU_seqs <- readDNAStringSet(file = "OTU_reps_fasta_002.fasta", format = "fasta", 
  nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

# import mothur output into phyloseq
mothur_phyloseq_object <- import_mothur(mothur_list_file = NULL,
  mothur_group_file = NULL, mothur_tree_file = NULL, cutoff = NULL,
  mothur_shared_file = our_shared_file,
  mothur_constaxonomy_file = our_cons_taxonomy_file,
  parseFunction = parse_taxonomy_default)

# now all files are imported, we can adjust them to our needs
# add further taxonomy columns "OTU" and "wholetax" and adjust column names
wholetax <- do.call(paste, c(as.data.frame(tax_table(mothur_phyloseq_object))
  [c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6")], sep = "_"))
tax_table(mothur_phyloseq_object) <- cbind(tax_table(mothur_phyloseq_object),
  rownames(tax_table(mothur_phyloseq_object)), wholetax)
colnames(tax_table(mothur_phyloseq_object)) <- c(
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "otu_id",
  "wholetax")
                                    
# phyloseq_object variations for different analyses
mothur_full <- merge_phyloseq(mothur_phyloseq_object, metafile2, refseq(OTU_seqs))
mothur_1 <- filter_taxa(mothur_full, function (x) {sum(x > 1) >= 1}, prune = TRUE)
mothur_1_melt <- psmelt(mothur_1)
mothur_full_ra <- transform_sample_counts(mothur_full, function(x) {(x / sum(x)) * 100})
mothur_ra_0.01 <- filter_taxa(mothur_full_ra, function (x) {sum(x > 0.01) >= 1}, prune = TRUE)
mothur_ra_melt <- psmelt(mothur_ra_0.01)
mothur_ra_melt_mean <- aggregate(Abundance ~ OTU + time + days + new_day
  + treatment + nucleic_acid + habitat + disturbance + glyphosate + glyphosate_gone 
  + condition_diversity + kingdom + phylum + class + order + family + genus + otu_id,
  data = mothur_ra_melt, mean)

# store phyloseq objects for other scripts
save.image("glyphosate_mothur_in_phyloseq.RData")

sessionInfo()