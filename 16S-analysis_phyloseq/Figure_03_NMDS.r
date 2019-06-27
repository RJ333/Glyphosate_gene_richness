#!/usr/bin/env Rscript
library(ggplot2)
library(phyloseq)
library(vegan)
library(gridExtra)

setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("glyphosate_mothur_in_phyloseq.RData")

phyloseq_for_nmds <- filter_taxa(mothur_full, function (x) {sum(x > 2) >= 1}, prune = TRUE)
phyloseq_for_nmds_rel_abund <- transform_sample_counts(phyloseq_for_nmds, function(x) {(x / sum(x)) * 100})

# arguments for the subsetting function
phyloseq_object <- phyloseq_for_nmds_rel_abund
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
threshold <- 0
days <- 40

# define a function to obtain sample subsets from the phyloseq object 
# per combination of habitat, nucleic acid, days and minimum required reads per OTU
get_sample_subsets <- function(phyloseq_object, nucleic_acid, habitat, days, threshold) {
  sample_subset <- sample_data(phyloseq_object)[ which(sample_data(phyloseq_object)$nucleic_acid == nucleic_acid &
    sample_data(phyloseq_object)$habitat == habitat & sample_data(phyloseq_object)$days > days),]
  phyloseq_subset <- merge_phyloseq(tax_table(phyloseq_object),
    otu_table(phyloseq_object),
    refseq(phyloseq_object),
    sample_subset)
  phyloseq_subset2 <- filter_taxa(phyloseq_subset, function (x) {sum(x > threshold) >= 1 }, prune = TRUE)
  return(phyloseq_subset2)
}

# here we pass the arguments for subsetting over two for loops
# to create all possible combinations of habitat, nucleic acid etc.
# the subsets are stored within a list, which has to be empty before running the loops 
sample_subset_list <- list()
if(length(sample_subset_list) == 0) {
  for (acid in acids) {
	for (habitat in habitats) {
	  print(paste0("nucleic_acid is ", acid, " and habitat is ", habitat))
	  tmp <- get_sample_subsets(phyloseq_object = phyloseq_object,
		nucleic_acid = acid, habitat = habitat, threshold = threshold, days = days)
	  sample_data(tmp)$days <- as.factor(sample_data(tmp)$days)
	  sample_data(tmp)$new_day <- as.factor(sample_data(tmp)$new_day)
	  sample_subset_list[[paste(habitat, acid, "min reads per OTU", threshold,
		sep = " ")]] <- tmp
    }
  }
  print(sample_subset_list)
} else {
  print("list is not empty, abort to prevend appending...")
}

# create a list where the distance metrics for the sample subsets are stored
ordination_nmds <- list()
ordination_nmds <- lapply(sample_subset_list, ordinate, method = "NMDS",
  dist = "bray", try = 100, autotransform = TRUE)

# test hypothesis that microcosms are different
set.seed(1)

# Calculate bray curtis distance matrix and collect sample data in lists
sample_distances <- list()
sample_data_list <- list()
sample_distances <- lapply(sample_subset_list, phyloseq::distance, method = "bray")
sample_data_list <- lapply(sample_subset_list, function(x) data.frame(sample_data(x)))

# Adonis test against hardcoded "treatment", instead passing "treatment" as argument did not work
adonis_results <- list()
adonis_results <- mapply(function(distance_matrix, sample_data) {
  adonis(distance_matrix ~ treatment, data = data.frame(sample_data))
}, distance_matrix = sample_distances, sample_data = sample_data_list, SIMPLIFY = FALSE)

# Homogeneity of dispersion test 
# tests if dispersion might be reason for adonis results, should confirm null hypothesis
beta_list <- list()
beta_list <- mapply(function(distance_matrix, sample_data) {
  betadisper(distance_matrix, sample_data$treatment)
}, distance_matrix = sample_distances, sample_data = sample_data_list, SIMPLIFY = FALSE)
dispersion_list <- list()
dispersion_list <- lapply(beta_list, permutest)

# generate NMDS plots with distance information from the ordination list
nmds_ordination_plots <- list()
counter <- 0
if(length(nmds_ordination_plots) == 0 & all.equal(counter, 0)) {
  nmds_ordination_plots <- mapply(function(x,y) {
    counter <<- counter + 1
    plot_ordination(x, y, type = "sample", color = "days", shape = "treatment") +
      geom_polygon(aes(fill = disturbance), alpha = 0.3, size = 0.01) +
      geom_point(aes(colour = treatment), colour = "black", size = 4, alpha = 0.7) +
      scale_shape_manual(values = c("glyph" = 16, "control" = 17), name = "Microcosm  ",
        breaks = c("glyph", "control"), labels = c("Treatment", "Control")) +
      scale_fill_manual(values = c("high" = "red3", "low" = "darkorange2", "none" = "forestgreen"),
        name = "Glyphosate\nconcentration", breaks = c("high", "low", "none"),
        labels = c("> 5 µM", "< 5 µM", "0 µM")) +
      guides(color = FALSE,
        shape = FALSE, fill = FALSE) +
      coord_cartesian(ylim = c(-0.77, 0.95), xlim = c(-0.9, 0.7)) +
      #ggtitle(names(sample_subset_list)[counter]) +
      geom_text(aes(label = new_day), colour = "white", size = 2.5) +
      theme_bw() +
      theme(panel.grid.major = element_line(colour = NA, size = 0.2),
        panel.grid.minor = element_line(colour = NA, size = 0.5),
        axis.text = element_text(size = 18),
        #axis.title = element_text(size = 20, face = "bold"),
        axis.title = element_blank(),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 13))
}, x = sample_subset_list, y = ordination_nmds, SIMPLIFY = FALSE)
} else {
  print(paste("list is not empty, or counter not 0 (counter is", counter,"), 
    abort to prevend appending..."))
}

# rotate the last plot to match the orientation of the others
nmds_ordination_plots[[4]] <- nmds_ordination_plots[[4]] + 
  scale_x_reverse() +
  scale_y_reverse()
  
# plot the NMDS ordinations in a 2 x 2 array
do.call("grid.arrange", c(nmds_ordination_plots[c(1, 3, 2, 4)], nrow = 2))

# write the NMDS ordinations array to an object
NMDS_array <- do.call("arrangeGrob", c(nmds_ordination_plots[c(1, 3, 2, 4)], nrow = 2))
ggsave(NMDS_array, file = paste(plot_path, "Figure_3_NMDS_size4.png", sep = ""),
  height = 9, width = 9)
scp -r -i /drives/d/ssh/denbi.key centos@193.196.20.93:/data/projects/glyphosate/reads/mothur_processed/plots/Figure_3_NMDS*.png /drives/d/denbi/chandler/ordination
