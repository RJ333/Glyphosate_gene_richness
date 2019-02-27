# Set the working dir with shared file, constaxonomy, sample and OTU_rep files in it
setwd("/data/projects/glyphosate/analysis_16S/dada2/")

# mkdir for plots
plot_path <- "./plots/"

# load libraries
library(scales)
library(ggplot2)
library(phyloseq)
library(DESeq2)

## TO DO: copy data files from cloud for OTU reps and tree

# import mothur output into phyloseq
mothur_ps <- import_mothur(mothur_list_file = NULL, 
						   mothur_group_file = NULL,
						   mothur_tree_file = NULL, 
						   cutoff = NULL, 
						   mothur_shared_file = "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.shared",
						   mothur_constaxonomy_file = "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.0.02.cons.taxonomy", 
						   parseFunction = parse_taxonomy_default)
# or load from workspace
load("mothur_glyph_002.RData")   
						   
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
metafile <- read.delim("metafile.tsv", 
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

############## Alpha diversity (including singletons)
erich_mothur <- estimate_richness(mothur_full, measures = c("Observed", 
															"Chao1", 
															"ACE", 
															"Shannon", 
															"Simpson", 
															"InvSimpson", 
															"Fisher"))
erich_mothur_meta <- cbind(erich_mothur, sample_data(mothur_full)[,c(1:7)])

# reorder and rename factor levels for plotting, adjust lab names
erich_mothur_meta$habitat <- relevel(erich_mothur_meta$habitat, "water")
erich_mothur_meta$nucleic_acid <- relevel(erich_mothur_meta$nucleic_acid, "dna")
labs_nucleic_acid <- c("DNA", "RNA")
labs_habitat <- c("Free-living", "Biofilm")
levels(erich_mothur_meta$habitat) <- labs_habitat
levels(erich_mothur_meta$nucleic_acid) <- labs_nucleic_acid

# plotting Alpha diversity

shannon_plot <- ggplot(erich_mothur_meta, aes(x = new_day, 
											  y = Shannon, 
											  colour = treatment)) + 
	geom_point(alpha = 0.8, size = 4) +
	geom_vline(aes(xintercept = 1), 
			   linetype = "dashed", 
			   size = 1.2) +
	stat_summary(aes(colour = treatment), 
				 fun.y = "mean",  
				 geom = "line",
				 alpha = 0.75,
				 size = 2) +
	scale_colour_manual(values = c("glyph" = "black", 
								   "control" = "grey50"), 
						name = "Microcosm  ", 
						breaks = c("glyph", 
								   "control"), 
						labels = c("Treatment", 
								   "Control")) +
	coord_cartesian(ylim = c(1, 3)) +
	theme_bw() +
	theme(axis.text = element_text(size = 18),
		  axis.title = element_text(size = 20, face = "bold"),
		  legend.title = element_text(size = 15, face = "bold"), 
		  legend.text = element_text(size = 13),
		  panel.grid.major = element_line(colour = NA, size = 0.2),
		  panel.grid.minor = element_line(colour = NA, size = 0.5),
		  strip.text.x = element_text(size = 15, face = "bold")
		  ) +
  	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	labs(x = "Days", y = "Shannon index") +
	facet_wrap(~ habitat + nucleic_acid)

ggsave(shannon_plot, file = paste(plot_path, "Figure_4_Shannon_DNA_RNA.png", 
								  sep = ""),
								  height = 10,
								  width = 14)

############## community overview bar plot

# transform into relative abundance, displayed in percentage!
mothur_full_ra <- transform_sample_counts(mothur_full, function(x){(x / sum(x)) * 100})
# remove low abundant OTUs
mothur_ra_0.01 <- filter_taxa(mothur_full_ra, function (x) {sum(x > 0.01) >= 1}, prune = TRUE)
# melt into long format for plotting
mothur_ra_melt <- psmelt(mothur_ra_0.01)
# calculate mean of technical replicates
mothur_ra_melt_mean <- aggregate(Abundance ~ OTU + time + days + new_day
								+ treatment + nucleic_acid + habitat + disturbance 
								+ glyphosate + glyphosate_gone + condition_diversity +
								+ kingdom + phylum + class + order + family + genus + otu_id, 
								data = mothur_ra_melt, 
								mean)
# get data per OTU, setting threshold for samples and clusters
community_subset <- droplevels(subset(mothur_ra_melt_mean, days > 40 
							& Abundance > 0.15
							& habitat == "water" 
							& treatment == "glyph"))
# check required number of colours per order and number of classes
length(levels(droplevels(community_subset$class)))
length(levels(droplevels(community_subset$order))) 
# sort orders based on phylogenetic class for plotting
community_subset$order <- factor(community_subset$order, 
									   # alphaproteos
							levels = c("Caulobacterales",
									   "Rhizobiales",
									   "Rhodobacterales",
									   "Rhodospirillales",
									   "Sneathiellales",
									   "Sphingomonadales",
									   "Parvibaculales",
									   "Thalassobaculales",
									   # gammaproteos
									   "Alteromonadales",
									   "Betaproteobacteriales",
									   "Pseudomonadales",
									   # bacteroidetes/
									   "Chitinophagales",
									   "Sphingobacteriales",
									   "Flavobacteriales"))
# assign specific colour to make plot distuingishable
fill_values <- c("Alteromonadales" = "orange",
					"Betaproteobacteriales" = "pink",
					"Caulobacterales" = "black",
					"Chitinophagales" = "purple",
					"Flavobacteriales" = "green",
					"Sneathiellales" = "white",
					"Parvibaculales" = "green3",
					"Pseudomonadales" = "grey30",
					"Rhizobiales" = "red",
					"Rhodobacterales" = "lightblue",
					"Rhodospirillales" = "yellow",
					"Sphingobacteriales" = "darkred",
					"Sphingomonadales" = "grey",
					"Thalassobaculales" = "blue2")	

# plotting all selected clusters in bar plot ordered by class 
# and displaying orders over time for DNA and RNA
community_plot <- ggplot(community_subset, aes(x = new_day, group = order)) +
	scale_fill_manual(breaks = levels(community_subset$order), 
				      values = fill_values) +
	geom_bar(data = subset(community_subset, nucleic_acid == "dna" & 
                                             treatment == "glyph"),
			 aes(x = new_day - 0.5, 
				 y = Abundance), 
			 fill = "black", 
		     width = 0.9, 
		     stat = "sum") +
	geom_bar(data = subset(community_subset, nucleic_acid == "dna" & 
                                             treatment == "glyph"),
			 aes(x = new_day - 0.5, 
				 y = Abundance, 
				 fill = order), 
			 width = 0.6, 
			 stat = "identity") +
	geom_bar(data = subset(community_subset, nucleic_acid == "cdna" & 
                                             treatment == "glyph"),
			 aes(x = new_day + 0.5, 
				 y = Abundance), 
			 fill = "black", 
			 width = 0.9, 
			 stat = "sum") +
	geom_bar(data = subset(community_subset, nucleic_acid == "cdna" & 
                                             treatment == "glyph"),
			 aes(x = new_day + 0.5, 
				 y = Abundance, 
				 fill = order), 
			 width = 0.6, 
			 stat = "identity") +
	geom_vline(data = subset(community_subset, treatment == "glyph"),
			   aes(xintercept = 1.5),
			   linetype = "dashed", size = 1.2) +
	guides(colour = FALSE, 
		   size = FALSE, 
		   width = FALSE,
		   fill = guide_legend(ncol = 1,
							   keyheight = 1.5,
							   label.theme = element_text(size = 15,
														  face = "italic",
														  angle = 0),
											(title = NULL))) +
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	scale_y_continuous(expand = c(0,0)) +
	theme_bw() +
	theme(axis.text = element_text(size = 17)) +
	theme(axis.title = element_text(size = 20,
									face = "bold")) +
	theme(legend.background = element_rect(fill = "grey90", 
										   linetype = "solid")) +
	theme(panel.grid.major = element_line(colour = NA, 
										  size = 0.2)) +
	theme(panel.grid.minor = element_line(colour = NA, 
										  size = 0.5)) +
	labs(x = "Days", 
		 y = "Relative abundance [%]") +
  annotate("text", 
		   x = -27, 
		   y = 90, 
		   label = "a)", 
		   color = "black", 
		   size = 6, 
		   angle = 0, 
		   fontface = "bold") +
  annotate("text", 
		   x = -22.5, 
		   y = 90, 
		   label = "b)", 
		   color = "black", 
		   size = 6, 
		   angle = 0, 
		   fontface = "bold")

ggsave(community_plot, file = paste(plot_path, "Figure_4_relative_community_overview.png", 
                                    sep = ""),
                                    width = 16, 
                                    height = 8)

############# OTU abundance plot                                    
                                    
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