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
# remove low abundant OTUs (you may decrease the threshold, but it will increase melting time)
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

# define subset function for specific phyloseq-object
get_current_otu_data <- function(x) {
	subset(mothur_ra_melt, OTU == x)
}

# list of OTUs mentioned in paper and supplement 
OTU_list <- c("Otu000007",
              # "Otu000011",
              # "Otu000018",
              # "Otu000025",
              # "Otu000032",
              # "Otu000036",
              # "Otu000037",
              # "Otu000038",
              # "Otu000023",
              # "Otu000046",
              # "Otu000049",
              # "Otu000056",
              # "Otu000058",
              # "Otu000059",
              # "Otu000070",
              # "Otu000072",
              # "Otu000078",
              # "Otu000094",
              # "Otu000109",
              # "Otu000129",
              # "Otu000139",
              # "Otu000176",
              # "Otu000191",
              # "Otu000320",
              # "Otu000098",
              # "Otu000042",
              # "Otu000044",
              # "Otu000006",
              "Otu000001"
              )
              
# rename for plotting
labs_habitat <- c("Biofilm", "Free-living")

# run a for loop to plot each OTU in list with own title and file name
for (i in OTU_list){
    current_otu_data <- get_current_otu_data(i)
    print(paste("OTU is", i))

species_title <- unique(paste(current_otu_data$family, 
							  current_otu_data$genus, 
							  current_otu_data$OTU, 
							  sep = " "))
							  
current_otu_data$treatment2 <- factor(current_otu_data$treatment, 
									  labels = c("Control", "Treatment"))

levels(current_otu_data$habitat) <- labs_habitat
                                      
current_plot <- ggplot(data = current_otu_data, 
	                   aes(x = days - 69, 
						   y = Abundance, 
						   group = nucleic_acid, 
						   lty = nucleic_acid)) + 
	geom_vline(aes(xintercept = 0), 
			   linetype = "dashed", 
			   size = 1.2) +
	geom_point(data = subset(current_otu_data, treatment == "control"), 
		       aes(colour = treatment), 
			   alpha = 1) +
	stat_summary(data = subset(current_otu_data, treatment == "control"), 
	             aes(colour = treatment), 
				 fun.y = "mean",  
				 geom = "line", 
				 size = 2, 
				 alpha = 1) +
	stat_summary(data = subset(current_otu_data, treatment == "glyph"), 
	             aes(colour = treatment), 
				 fun.y = "mean",  
				 geom = "line", 
				 size = 2) +
	geom_point(data = subset(current_otu_data, treatment == "glyph"), 
	           aes(colour = treatment)) +
scale_linetype_manual(values = c("dna" = 1, 
									 "cdna" = 6), 
						name = "Nucleic acid  ", 
						breaks = c("cdna", 
								   "dna"), 
						labels = c("16S rRNA", 
								   "16S rRNA gene")) +
	scale_colour_manual(values = c("glyph" = "black", 
								   "control" = "grey50"), 
						name = "Microcosm  ", 
						breaks = c("glyph", 
								   "control"), 
						labels = c("Treatment", 
								   "Control")) +
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	theme_bw() +
	ggtitle(species_title) +
	theme(axis.text = element_text(size = 18),
		  axis.title = element_text(size = 20, face = "bold"),
          panel.grid.major = element_line(colour = NA, size = 0.2),
          panel.grid.minor = element_line(colour = NA, size = 0.5),
          strip.text.x = element_text(size = 15, face = "bold")) +
	labs(x = "Days", y = "Relative abundance [%]") +
	facet_wrap(~ habitat, scales = "free")
	# ggsave(current_plot, file = paste(plot_path, 
									  # species_title, 
									  # ".png", 
									  # sep = ""), 
						 # width = 13, 
						 # height = 7)
	print(current_plot)
}

# factorize OTUs
mothur_ra_melt$OTU <- as.factor(mothur_ra_melt$OTU)

save.image("mothur_glyph.RData")