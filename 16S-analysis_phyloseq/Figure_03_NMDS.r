#!/usr/bin/env Rscript
library(ggplot2)
library(phyloseq)
library(vegan)
library(gridExtra)
library(grid)
library(dplyr)

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
# extract the information from the list for the subsets into one long dataframe for plotting
nmds_1 <- merge(sample_data(sample_subset_list[[1]]), ordination_nmds[[1]]$points, by = "row.names")
nmds_2 <- merge(sample_data(sample_subset_list[[2]]), ordination_nmds[[2]]$points, by = "row.names")
nmds_3 <- merge(sample_data(sample_subset_list[[3]]), ordination_nmds[[3]]$points, by = "row.names")
nmds_4 <- merge(sample_data(sample_subset_list[[4]]), ordination_nmds[[4]]$points, by = "row.names")
nmds_df <- rbind(nmds_1, nmds_2, nmds_3, nmds_4)
row.names(nmds_df) <- nmds_df$Row.names
nmds_df <- nmds_df[, -c(1:3, 10:14)]

# calculate center of polygon per triplicate
nmds_mean <- nmds_df %>%
  group_by(treatment, habitat, nucleic_acid, new_day) %>%
  mutate(MDS1_mean = mean(MDS1)) %>% mutate(MDS2_mean = mean(MDS2))

levels(nmds_mean$nucleic_acid) <- c("16S rRNA gene", "16S rRNA") 
levels(nmds_mean$habitat) <- c("Free-living", "Biofilm")

all_NMDS <- ggplot(nmds_mean, aes(x = MDS1, y = MDS2, shape = treatment, colour = disturbance)) +
  geom_point(aes(shape = treatment), colour = "black", size = 0.5, alpha = 0.7) +
  geom_polygon(aes(fill = disturbance, group = interaction(new_day, treatment)), 
    alpha = 0.5, size = 0.3) +
  scale_shape_manual(values = c("glyph" = 16, "control" = 17), name = "Microcosm  ",
    breaks = c("glyph", "control"), labels = c("Treatment", "Control")) +
  scale_fill_manual(values = c("high" = "red3", "low" = "darkorange2", "none" = "forestgreen"),
    name = "Glyphosate", breaks = c("high", "low", "none"),
    labels = c("> 5 µM", "< 5 µM", "0 µM")) +
  scale_colour_manual(values = c("high" = "red3", "low" = "darkorange2", "none" = "forestgreen"),
    name = "Glyphosate", breaks = c("high", "low", "none"),
    labels = c("> 5 µM", "< 5 µM", "0 µM")) +
  guides(color = FALSE) +
  guides(shape = guide_legend(override.aes = list(size = 4))) +
  coord_cartesian(ylim = c(-0.77, 0.95), xlim = c(-0.9, 0.7)) +
  geom_text(aes(x = MDS1_mean, y = MDS2_mean, label = new_day), colour = "black", size = 2, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 14, face = "bold")) +
    labs(x = "NMDS1", y = "NMDS2") +
    facet_wrap(~habitat + nucleic_acid)

ggsave(all_NMDS, file = paste(plot_path, "Figure_3_NMDS_facet.pdf", sep = ""),
  device = "pdf", width = 18.0, height = 16, dpi = 300, unit = "cm")

# biofilm plots need to be reversed to be in similar direction as water plots
# copy the code above and paste lower row on old plot (e.g. with Inkscape)

# all_NMDS_reversed <-    
  # scale_x_reverse() +
  # scale_y_reverse() +

# ggsave(all_NMDS_reversed, file = paste(plot_path, "Figure_3_NMDS_facet_reversed.pdf", sep = ""),
   # device = "pdf", width = 18.0, height = 16, dpi = 300, unit = "cm")