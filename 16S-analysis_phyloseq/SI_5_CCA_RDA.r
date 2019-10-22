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

#### CCA

# create a list where the distance metrics for the sample subsets are stored
ordination_cca <- list()
ordination_cca <- lapply(sample_subset_list, ordinate, method = "CCA",
  formula = sample_subset_list[[]] ~ glyphosate, try = 100, autotransform = TRUE)
  
# generate  plots with distance information from the ordination list
# extract the information from the list for the subsets into one long dataframe for plotting
# cca scores found with help from  http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html
cca_1 <- merge(sample_data(sample_subset_list[[1]]), scores(ordination_cca[[1]], display = "wa"), by = "row.names")
cca_2 <- merge(sample_data(sample_subset_list[[2]]), scores(ordination_cca[[2]], display = "wa"), by = "row.names")
cca_3 <- merge(sample_data(sample_subset_list[[3]]), scores(ordination_cca[[3]], display = "wa"), by = "row.names")
cca_4 <- merge(sample_data(sample_subset_list[[4]]), scores(ordination_cca[[4]], display = "wa"), by = "row.names")
cca_df <- rbind(cca_1, cca_2, cca_3, cca_4)
row.names(cca_df) <- cca_df$Row.names
cca_df <- cca_df[, -c(1:3, 10:14)]

# calculate center of polygon per triplicate
cca_mean <- cca_df %>%
  group_by(treatment, habitat, nucleic_acid, new_day) %>%
  mutate(CCA_mean = mean(CCA1)) %>% mutate(CA_mean = mean(CA1))

levels(cca_mean$nucleic_acid) <- c("16S rRNA gene", "16S rRNA") 
levels(cca_mean$habitat) <- c("Free-living", "Biofilm")

all_CCA <- ggplot(cca_mean, aes(x = CCA1, y = CA1, shape = treatment, colour = disturbance)) +
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
  #coord_cartesian(ylim = c(-0.77, 0.95), xlim = c(-0.9, 0.7)) +
  geom_text(aes(x = CCA_mean - 0.11, y = CA_mean, label = new_day), colour = "black", size = 2.5, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 14, face = "bold")) +
    labs(x = "CCA1", y = "CA1") +
    facet_wrap(~habitat + nucleic_acid, scales = "free")

ggsave(all_CCA, file = paste(plot_path, "Figure_x_CCA_facet.svg", sep = ""),
  device = "svg", width = 18.0, height = 16, dpi = 300, unit = "cm")
  
#### RDA
  
# create a list where the distance metrics for the sample subsets are stored
ordination_rda <- list()
ordination_rda <- lapply(sample_subset_list, ordinate, method = "RDA",
  formula = sample_subset_list[[]] ~ glyphosate, try = 100, autotransform = TRUE)
  
# generate  plots with distance information from the ordination list
# extract the information from the list for the subsets into one long dataframe for plotting
# rda scores found with help from  http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html
rda_1 <- merge(sample_data(sample_subset_list[[1]]), scores(ordination_rda[[1]], display = "wa"), by = "row.names")
rda_2 <- merge(sample_data(sample_subset_list[[2]]), scores(ordination_rda[[2]], display = "wa"), by = "row.names")
rda_3 <- merge(sample_data(sample_subset_list[[3]]), scores(ordination_rda[[3]], display = "wa"), by = "row.names")
rda_4 <- merge(sample_data(sample_subset_list[[4]]), scores(ordination_rda[[4]], display = "wa"), by = "row.names")
rda_df <- rbind(rda_1, rda_2, rda_3, rda_4)
row.names(rda_df) <- rda_df$Row.names
rda_df <- rda_df[, -c(1:3, 10:14)]

# calculate center of polygon per triplicate
rda_mean <- rda_df %>%
  group_by(treatment, habitat, nucleic_acid, new_day) %>%
  mutate(RDA_mean = mean(RDA1)) %>% mutate(CA_mean = mean(PC1))

levels(rda_mean$nucleic_acid) <- c("16S rRNA gene", "16S rRNA") 
levels(rda_mean$habitat) <- c("Free-living", "Biofilm")

all_RDA <- ggplot(rda_mean, aes(x = RDA1, y = PC1, shape = treatment, colour = disturbance)) +
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
  #coord_cartesian(ylim = c(-0.77, 0.95), xlim = c(-0.9, 0.7)) +
  geom_text(aes(x = RDA_mean - 0.11, y = CA_mean, label = new_day), colour = "black", size = 2.5, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 14, face = "bold")) +
    labs(x = "RDA1", y = "PC1") +
    facet_wrap(~habitat + nucleic_acid, scales = "free")

ggsave(all_RDA, file = paste(plot_path, "Figure_x_RDA_facet.svg", sep = ""),
  device = "svg", width = 18.0, height = 16, dpi = 300, unit = "cm")