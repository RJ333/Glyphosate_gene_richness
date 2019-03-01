#!/usr/bin/env Rscript

# load libraries
library(scales)
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(gridExtra)

# Set the working dir with shared file, constaxonomy, sample and OTU_rep files in it
setwd("/data/projects/glyphosate/analysis_16S/dada2/")

# set path for plots, create dir if not existant
plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

# mothur files that need to be imported
# shared file is the OTU count table
# constaxonomy contains the taxonomy of the OTUs
our_shared_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.shared"
our_cons_taxonomy_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.0.02.cons.taxonomy"

# import mothur output into phyloseq
mothur_ps <- import_mothur(mothur_list_file = NULL, 
						   mothur_group_file = NULL,
						   mothur_tree_file = NULL, 
						   cutoff = NULL, 
						   mothur_shared_file = our_shared_file,
						   mothur_constaxonomy_file = our_cons_taxonomy_file, 
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

# read table for cell count plot
cell_counts_glyph <- read.csv("cell_counts_glyph.csv", sep = ";")

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
# remove singletons
mothur_1 <- filter_taxa(mothur_full, function (x) {sum(x > 1) >= 1}, prune = TRUE)
# transform into relative abundance, displayed in percentage!
mothur_full_ra <- transform_sample_counts(mothur_full, function(x){(x / sum(x)) * 100})
# remove low abundant OTUs (you may decrease the threshold, but it will increase melting time)
mothur_ra_0.01 <- filter_taxa(mothur_full_ra, function (x) {sum(x > 0.01) >= 1}, prune = TRUE)
# melt into long format for plotting
mothur_ra_melt <- psmelt(mothur_ra_0.01)
mothur_1_melt <- psmelt(mothur_1)
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
strip_text_habitat <- c("Biofilm", "Free-living")

# run a for loop to plot each OTU in list with own title and file name
for (i in OTU_list){
    current_otu_data <- get_current_otu_data(i)
    print(paste("OTU is", i))

species_title <- unique(paste(current_otu_data$family, 
							  current_otu_data$genus, 
							  current_otu_data$OTU, 
							  sep = " "))

levels(current_otu_data$habitat) <- strip_text_habitat
                                      
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

######## Total cell counts and glyphosate concentration

cell_counts_glyph_0 <- subset(cell_counts_glyph, new_day >= -7)

cell_counts_glyph_plot <- ggplot(cell_counts_glyph_0, aes(x = new_day, 
                                colour = treatment, 
                                linetype = treatment, 
                                group = treatment)) +
	geom_errorbar(aes(ymin = cells_ml - cells_se, 
                      ymax = cells_ml + cells_se),
                  linetype = 1, width = 2, size = 1.2, alpha = 0.7) +
	geom_errorbar(aes(ymin = (glyph_mg_L - glyph_se) * 2200000, 
                      ymax = (glyph_mg_L + glyph_se) * 2200000), 
                  linetype = 1, width = 1, size = 1.0, alpha = 0.6) +
    geom_line(aes(y = cells_ml, 
                  group = treatment, 
                  colour = treatment), 
              linetype = 1, size = 2, alpha = 0.8) +
	geom_point(aes(y = glyph_theor * 2200000, 
                   shape = "glyph"),  
               alpha = 0.5, size = 5) +
	geom_point(aes(y = glyph_mg_L * 2200000, 
                   shape = "glyph_deg"), 
               alpha = 0.7, size = 4) +

	geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1.5, alpha = 0.5)+
		
	scale_y_continuous(label =  function(x) {ifelse(x == 0, "0", parse(text = gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))}, 
        sec.axis = sec_axis(~ . / 2200000, name = expression(paste("Glyphosate concentration  ", bgroup("[",mg~L^{-1},"]")))))+
	scale_colour_manual(values = c("control" = "grey70", "glyph" = "black"),
                        name = "Microcosm",
                        breaks = c("control", "glyph"),
                        labels = c("Control", "Treatment"))+
	scale_shape_manual(values = c("glyph" = 2, "glyph_deg" = 17),
                       name = "Glyphosate decrease by",
                       breaks = c("glyph", "glyph_deg"),
                       labels = c("calculated\ndilution", "measured\nconcentration"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
    theme_bw() +
	theme(panel.grid.major = element_line(colour = NA, size = 0.2),
          panel.grid.minor = element_line(colour = NA, size = 0.5),
          axis.title = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(angle = 90, vjust = 1),
          axis.text = element_text(size = 18),
          legend.title=element_text(size = 14), 
          legend.text=element_text(size = 12)) +
	labs(x = "Days",
		 y = expression(paste("Total cell count  ",bgroup("[",cells~mL^{-1},"]"))))
			 
ggsave(cell_counts_glyph_plot, file = paste(plot_path, "Figure_1_cellcounts_glyph.png", 
                                            sep = ""),
                               width = 14, 
                               height = 10)


################### NMDS plots

# exclude OTUs with less than 3 reads and transform to relative abundance
mothur_nmds <- filter_taxa(mothur_full, function (x) {sum(x > 2) >= 1}, prune = TRUE)
mothur_nmds_ra <- transform_sample_counts(mothur_nmds, function(x){(x / sum(x)) * 100})

ps <- mothur_nmds_ra
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
threshold <- 0


# define a function to obtain sample subsets per combination of habitat, nucleic acid, days and treatment
get_sample_subsets <- function(ps, nucleic_acid, habitat, days, threshold){
	sample_subset <- sample_data(ps)[ which(sample_data(ps)$nucleic_acid == nucleic_acid & 
											sample_data(ps)$habitat == habitat & 
											sample_data(ps)$days > days),]
	phy_subset <- merge_phyloseq(tax_table(ps), 
								 otu_table(ps),
								 refseq(ps),
								 sample_subset)
	phy_subset2 <- filter_taxa(phy_subset, function (x) {sum(x > threshold) >= 1 }, prune = TRUE)
	return(phy_subset2)
}

sample_subset_list <- list() 
if(length(sample_subset_list) == 0) {
		for (acid in acids) {
			for (habitat in habitats) {
				print(paste0("nucleic_acid is ", acid, " and habitat is ", 
							 habitat))
				tmp <-	get_sample_subsets(ps = ps, 
									   nucleic_acid = acid, 
									   habitat = habitat,  
									   threshold = threshold)
				sample_data(tmp)$days <- as.factor(sample_data(tmp)$days)					   
				sample_data(tmp)$new_day <- as.factor(sample_data(tmp)$new_day)
				sample_subset_list[[paste(habitat, 
										  acid, 
										  "min reads per OTU", 
										  threshold, 
										  sep = " ")]] <- tmp
			}
	}
print(sample_subset_list)
} else {
	print("list is not empty, abort to prevend appending...")
}

ordination_nmds <- list()
ordination_nmds <- lapply(sample_subset_list, ordinate, 
											  method = "NMDS", 
											  dist = "bray", 
											  try = 100, 
											  autotransform = TRUE)
                                              
# NMDS function
nmds_ordination_plots <- list()
counter <- 0
if(length(nmds_ordination_plots) == 0 & 
	all.equal(counter, 0)) {
nmds_ordination_plots <- mapply(function(x,y) {
						 counter <<- counter + 1 
						 plot_ordination(x, y, 
										 type = "sample",
										 color = "days",
										 shape = "treatment") + 
							geom_polygon(aes(fill = disturbance), alpha = 0.5, size = 0.01) + 
							geom_point(aes(colour = treatment), colour = "black", size = 4.5, alpha = 0.7) +
							scale_shape_manual(values = c("glyph" = 16, 
                                                          "control" = 17), 
                            name = "Microcosm  ", 
                            breaks = c("glyph", 
                                       "control"), 
                            labels = c("Treatment", 
                                       "Control")) +
							guides(color = FALSE, fill = FALSE, shape = FALSE) +
							coord_cartesian(ylim = c(-0.77, 0.95), xlim = c(-0.9, 0.7)) +
							#ggtitle(names(sample_subset_list)[counter]) +
							geom_text(aes(label = new_day), 
									  colour = "white", 
									  size = 2.5) +
							theme_bw() +
							theme(panel.grid.major = element_line(colour = NA, size = 0.2),
                                  panel.grid.minor = element_line(colour = NA, size = 0.5),
                                  axis.text = element_text(size = 18),
								  #axis.title = element_text(size = 20, face = "bold"),
								  axis.title = element_blank(),
                                  legend.title = element_text(size = 15, face = "bold"), 
								  legend.text = element_text(size = 13))
}, x = sample_subset_list, 
   y = ordination_nmds, 
   SIMPLIFY = FALSE)
} else {
	print(paste("list is not empty, or counter not 0 (counter is", counter, 
				"), abort to prevend appending..."))
}

do.call("grid.arrange", c(nmds_ordination_plots[c(1, 3, 2, 4)], nrow = 2))

g1 <- do.call("arrangeGrob", c(nmds_ordination_plots[c(1, 3, 2, 4)], nrow = 2))
ggsave(g1, file = paste(plot_path, "Figure_3_NMDS.png", 
								  sep = ""),
								  height = 10,
								  width = 10)
                
#################### DESeq2

# test variable is not allowed to contain NA
mothur_deseq <- subset_samples(mothur_full, !(is.na(condition)))

# these are the parameters passed to function
ps <- mothur_deseq
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
treatments <- c("glyph", "control")
threshold <- 0

# this is the function we call to split our data into different subsets
get_sample_subsets <- function(ps, nucleic_acid, habitat, treatment){
	sample_subset <- sample_data(ps)[ which(sample_data(ps)$nucleic_acid == nucleic_acid & 
											sample_data(ps)$habitat == habitat & 
											sample_data(ps)$treatment == treatment),]
	phy_subset <- merge_phyloseq(tax_table(ps), 
								 otu_table(ps),
								 sample_subset)
	phy_subset2 <- filter_taxa(phy_subset, function (x) {sum(x > 0) >= 1}, prune = TRUE)
	return(phy_subset2)
}

# this is the nested for loop which calls the subsetting function 
# for each combination of subsetting variables
deseq_subsets <- list() 
if(length(deseq_subsets) == 0) {
	for (treatment in treatments){
		for (acid in acids) {
			for (habitat in habitats) {
				print(paste0("nucleic_acid is ", acid, " and habitat is ", 
							 habitat, " and treatment is ", treatment))
				tmp <-	get_sample_subsets(ps = ps, 
									   nucleic_acid = acid, 
									   habitat = habitat, 
									   treatment = treatment)
				sample_data(tmp)$days <- as.factor(sample_data(tmp)$days)					   
				sample_data(tmp)$new_day <- as.factor(sample_data(tmp)$new_day)
				deseq_subsets[[paste(habitat,
									 treatment, 
									 acid, 
							   sep = "_")]] <- tmp
			}
		}
	}
print(deseq_subsets)
} else {
	print("list is not empty, abort to prevend appending...")
}

# we can now estimate or determine different diversity parameters on our subsets
deseq_tests <- list()						  
counter <- 0
if(length(deseq_tests) == 0 & 
	all.equal(counter, 0)) {						  
deseq_tests <- lapply(deseq_subsets, 
					   function(deseqs) {
								counter  <<- counter + 1
                                tmp = phyloseq_to_deseq2(deseqs, ~ condition)
								tmp$condition <- relevel(tmp$condition, ref = "untreated")
								tmp_dds = DESeq(tmp, test = "Wald", fitType = "parametric")								
})
} else {
	print(paste("list is not empty, or counter not 0 (counter is", counter, 
				"), abort to prevend appending..."))
}

# try cbind with mapply
sigtabs_list <- list()
if(length(sigtabs_list) == 0) {
sigtabs_list <- mapply(function(dds, ps) {res = results(dds, cooksCutoff = FALSE)
									alpha = 0.01
									sigtab = res[which(res$padj < alpha), ]
									sigtab = cbind(as(sigtab, "data.frame"), 
                                        as(tax_table(ps)[rownames(sigtab), ], "matrix"))
									print(head(sigtab))
									return(sigtab)
}, dds = deseq_tests, ps = deseq_subsets, SIMPLIFY = FALSE)
} else {
	print(paste("list is not empty, abort to prevent appending..."))
}

# all significant changes are now in the list
sigtabs_list
# sort by log fold change
sigs_ordered <- lapply(sigtabs_list, function(x) x[order(x$log2FoldChange),])

# create a vector of all identified OTUs in sigtabs_list (can be used for plotting abundances!)
deseq_otus <- row.names(sigtabs_list[[1]])
for (i in 2:8) {deseq_otus <- unique(append(deseq_otus, row.names(sigtabs_list[[i]])))}

               
################### Venn Diagrams
               
# define function to plot Venn diagram with 4 categories, here biofilm vs water column
fourway.Venn <- function(A,B,C,D,cat.names = c("Water\nDNA",
											   "Biofilm\nDNA",
											   "Water\nRNA",
											   "Biofilm\nRNA")){
  grid.newpage()
  area1 <- length(A)
  area2 <- length(B)
  area3 <- length(C)
  area4 <- length(D)
  n12<-length(Reduce(intersect, list(A,B)))
  n13<-length(Reduce(intersect, list(A,C)))
  n14<-length(Reduce(intersect, list(A,D)))
  n23<-length(Reduce(intersect, list(B,C)))
  n24<-length(Reduce(intersect, list(B,D)))
  n34<-length(Reduce(intersect, list(C,D)))
  n123<-length(Reduce(intersect, list(A,B,C)))
  n124<-length(Reduce(intersect, list(A,B,D)))
  n134<-length(Reduce(intersect, list(A,C,D)))
  n234<-length(Reduce(intersect, list(B,C,D)))
  n1234<-length(Reduce(intersect, list(A,B,C,D)))
  
venn.plot <- draw.quad.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  area4 = area4,
  n12 = n12,
  n13 = n13,
  n14 = n14,
  n23 = n23,
  n24 = n24,
  n34 = n34,
  n123 = n123,
  n124 = n124,
  n134 = n134,
  n234 = n234,
  n1234 = n1234,
  category = cat.names,
  cat.pos = c(0,180,0,200),
  fill = c("blue", "red", "green", "yellow"),
  alpha = .3,
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "green", "black")
)
grid.draw(venn.plot)
}

# factorize OTUs to count them
mothur_ra_melt$OTU <- as.factor(mothur_ra_melt$OTU)

# generate subsets to count OTUs per subset
# habitat * nucleic_acid
water_dna_genera <- subset(mothur_ra_melt, habitat == "water" & nucleic_acid == "dna" & Abundance > 0.05)
water_dna_unique_genera<-water_dna_genera[which(!duplicated(water_dna_genera[,"genus"])),]
nrow(water_dna_unique_genera)

water_cdna_genera <- subset(mothur_ra_melt, habitat == "water" & nucleic_acid == "cdna" & Abundance > 0.05)
water_cdna_unique_genera<-water_cdna_genera[which(!duplicated(water_cdna_genera[,"genus"])),]
nrow(water_cdna_unique_genera)

biofilm_dna_genera <- subset(mothur_ra_melt, habitat == "biofilm" & nucleic_acid == "dna" & Abundance > 0.05)
biofilm_dna_unique_genera<-biofilm_dna_genera[which(!duplicated(biofilm_dna_genera[,"genus"])),]
nrow(biofilm_dna_unique_genera)

biofilm_cdna_genera <- subset(mothur_ra_melt, habitat == "biofilm" & nucleic_acid == "cdna" & Abundance > 0.05)
biofilm_cdna_unique_genera<-biofilm_cdna_genera[which(!duplicated(biofilm_cdna_genera[,"genus"])),]
nrow(biofilm_cdna_unique_genera)

# plot Venn diagram
fourway.Venn(water_dna_unique_otus$genus,
			 biofilm_dna_unique_otus$genus,
			 water_cdna_unique_otus$genus,
			 biofilm_cdna_unique_otus$genus)
dev.copy(png, paste(plot_path, "Supplement_4wayVenn_nucleic_acids_genus_0.05.png"))
dev.off()

############################ OTU, sequence length and genera distribution

# library sizes are returned using
sample_sums(mothur_full)

# how many OTUs belong to which genus?
genus_distribution <- aggregate(Abundance ~ OTU + genus, 
								data = mothur_1_melt, 
								max)
# separated by habitat and nucleic acid								
genus_distribution2 <- aggregate(Abundance ~ OTU + habitat + genus + nucleic_acid, 
								data = mothur_1_melt, 
								max)

# OTUs per genus ordered by amount of OTUs								
otu_per_genus <- as.data.frame(table(genus_distribution$genus))
otu_per_genus[order(otu_per_genus$Freq),]
nrow(otu_per_genus)

# these are the arguments for the subsetting function
ps <- mothur_full
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
threshold <- 1
after_day <- 43

# in this list we store the different sample subsets, generated by the for loops
sample_subset_list <- list() 
if(length(sample_subset_list) == 0) {
	for (each_day in after_day){
		for (acid in acids) {
			for (habitat in habitats) {
				print(paste0("nucleic_acid is ", acid, " and habitat is ", 
							 habitat, " and first day is ", each_day))
				tmp <-	get_sample_subsets(ps = ps, 
									   nucleic_acid = acid, 
									   habitat = habitat, 
									   days = each_day, 
									   threshold = threshold)
				sample_data(tmp)$days <- as.factor(sample_data(tmp)$days)					   
				sample_data(tmp)$new_day <- as.factor(sample_data(tmp)$new_day)
				sample_subset_list[[paste(habitat, 
										  "after day", 
										  each_day, 
										  acid, 
										  "min reads per OTU", 
										  threshold, 
										  sep = " ")]] <- tmp
			}
		}
	}
print(sample_subset_list)
} else {
	print("list is not empty, abort to prevend appending...")
}

# the distribution of sequence length can be addressed using
table(width(refseq(sample_subset_list[[1]])))

sessionInfo()