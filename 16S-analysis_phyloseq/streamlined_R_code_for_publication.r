# Set the working dir with shared file, constaxonomy, sample and OTU_rep files in it
setwd("D:/Arbeit/denbi/chandler/")

# mkdir for plots
plot_path <- "./plots/"

# load libraries
library(scales)
library(ggplot2)
library(phyloseq)


# import mothur output into phyloseq
mothur_ps <- import_mothur(mothur_list_file = NULL, 
						   mothur_group_file = NULL,
						   mothur_tree_file = NULL, 
						   cutoff = NULL, 
						   mothur_shared_file = "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.shared",
						   mothur_constaxonomy_file = "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.0.02.cons.taxonomy", 
						   parseFunction = parse_taxonomy_default)
						   
# add taxonomy columns and adjust header
wholetax <- do.call(paste, c(as.data.frame(tax_table(mothur_ps))
                  [c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6")], 
				  sep = "_"))

tax_table(mothur_ps) <- cbind(tax_table(mothur_ps), 
							  rownames(tax_table(mothur_ps)), 
							  wholetax)
							  
colnames(tax_table(mothur_ps)) <- c("kingdom", 
									"phylum", 
									"class", 
									"order", 
									"family", 
									"genus", 
									"otu_id",
									"wholetax")
											
# read meta data, turn into phyloseq object, merge with existing ps object									
metafile <- read.delim("all_samples_with_metacond3.tsv", 
						row.names = 1, 
						header = TRUE,
						na.strings = "")
metafile <- sample_data(metafile)

# read OTU representative sequences

### for generation of file check gitlab #59

OTU_seqs <- readDNAStringSet(file = "OTU_reps_fasta_002.fasta", 
							  format = "fasta",
							  nrec = -1L, 
							  skip = 0L, 
							  seek.first.rec = FALSE, 
							  use.names = TRUE)
# add meta data and OTU representative seqs to phyloseq object
mothur_full <- merge_phyloseq(mothur_ps, metafile, refseq(OTU_seqs))

# remove OTUs with less than 2 reads in at least 1 sample
mothur_greater_1 <- filter_taxa(mothur_full, function (x) {sum(x > 1) >= 1}, prune = TRUE)
# transform counts into relative abundances, displayed as percentage
mothur_relative <- transform_sample_counts(mothur_greater_1, function(x){(x / sum(x)) * 100})
mothur_ra_melt <- psmelt(mothur_relative)
# add absolute abundance (rel abundance * total cell counts)
mothur_ra_melt$abs_Abundance <- (mothur_ra_melt$Abundance * mothur_ra_melt$cell_counts)/100
# factorize OTUs
mothur_ra_melt$OTU <- as.factor(mothur_ra_melt$OTU)

### depending on the plot, we need the parallels separately or averaged
# need to remove "Sample" to average the parallels
mothur_ra_melt_mean <- aggregate(cbind(Abundance, abs_Abundance) ~ OTU + time + days + new_day
								+ treatment + nucleic_acid + habitat + disturbance 
								+ cell_counts + glyphosate + glyphosate_gone + condition_diversity +
								+ kingdom + phylum + class + order + family + genus + otu_id + wholetax, 
								data = mothur_ra_melt, 
								mean)

save.image("mothur_glyph.RData")