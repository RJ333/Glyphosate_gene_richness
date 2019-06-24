#!/usr/bin/env Rscript
# load libraries
library(scales)
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(gridExtra)
library(vegan)
library(Biostrings)
library(tidyverse)

# Set the working dir
# it should contain shared file, constaxonomy file,
# meta file, cell count file and OTU_rep fasta file in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")

# set path for plots, create directory if not existant
plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

load("full_test.RData")

# add meta data and OTU representative seqs to phyloseq_object
mothur_full <- merge_phyloseq(mothur_phyloseq_object, metafile2, refseq(OTU_seqs))

# remove singletons
mothur_1 <- filter_taxa(mothur_full, function (x) {sum(x > 1) >= 1}, prune = TRUE)

# transform into relative abundance, displayed in percentage!
mothur_full_ra <- transform_sample_counts(mothur_full, function(x) {(x / sum(x)) * 100})

# remove low abundant OTUs (you may decrease the threshold, but it will increase melting time)
mothur_ra_0.01 <- filter_taxa(mothur_full_ra, function (x) {sum(x > 0.01) >= 1}, prune = TRUE)

# melt into long format for plotting
mothur_ra_melt <- psmelt(mothur_ra_0.01)
mothur_1_melt <- psmelt(mothur_1)

# aggregate all columns except for "parallels" to 
# calculate the mean abundance of technical replicates
mothur_ra_melt_mean <- aggregate(Abundance ~ OTU + time + days + new_day
  + treatment + nucleic_acid + habitat + disturbance + glyphosate + glyphosate_gone 
  + condition_diversity + kingdom + phylum + class + order + family + genus + otu_id,
  data = mothur_ra_melt, mean)
  
# get data per water OTU, setting threshold for samples and clusters
community_subset_water <- droplevels(subset(mothur_ra_melt_mean, days > 40
  & Abundance > 0.15 & habitat == "water"))

community_subset_water %>%
  group_by(new_day, order, nucleic_acid, treatment) %>%
  summarise(Abundance2 = sum(Abundance)) %>%
  ungroup() %>%
  # Replace missing values by 0
  spread(key = order, value = Abundance2) %>%
  replace(., is.na(.), 0) %>%
  gather(key = order, value = Abundance2, -c(new_day, nucleic_acid, treatment)) -> 
    order_sums_water
  
order_sums_water$order <- factor(order_sums_water$order, 
  levels = c(
    # alphaproteos
    "Caulobacterales", #
    "Rhizobiales", #
    "Rhodobacterales", #
    "Rhodospirillales", #
    "Sneathiellales", #
    "Sphingomonadales", #
    "Parvibaculales",
    "Thalassobaculales",
    # gammaproteos
    "Aeromonadales", #
    "Alteromonadales", #
    "Betaproteobacteriales", #
    "Gammaproteobacteria_Incertae_Sedis",
    "Oceanospirillales", #
    "Pseudomonadales", #
    "Xanthomonadales", #
    # Actinobacteria
    "Corynebacteriales",
    # Bacteroidia
    "Bacillales", #
    "Bacteroidia_unclassified",
    "Chitinophagales", #
    "Flavobacteriales", #
    "Sphingobacteriales",
    # Planctomycetacia
    "Planctomycetales", #
    "OM190_or", #
    # Verrucomicrobia
    "Verrucomicrobiales" #
    ))
    
# assign specific colour to make plot distuingishable
fill_values_water <- c("Aeromonadales" = "green",
  "Alteromonadales" = "#e6194B",
  "Bacillales" = "red",
  "Bacteroidia_unclassified" = "maroon2",
  "Betaproteobacteriales" = "#3cb44b",
  "Caulobacterales" = "#ffe119",
  "Chitinophagales" = "#4363d8",
  "Corynebacteriales" = "darkblue",
  "Cytophagales" = "blue",
  "Flavobacteriales" = "#f58231",
  "Gammaproteobacteria_Incertae_Sedis" = "black",
  "Gaiellales" = "black",
  "Oceanospirillales" = "maroon4",
  "OM190_or" = "grey80",
  "Opitutales" = "yellow",
  "Sneathiellales" = "#42d4f4",
  "Parvibaculales" = "#f032e6",
  "Planctomycetales" = "yellow",
  "Pseudomonadales" = "#fabebe",
  "Rhizobiales" = "#469990",
  "Rhodobacterales" = "#000000",
  "Rhodospirillales" = "#9A6324",
  "Sphingobacteriales" = "#fffac8",
  "Sphingomonadales" = "#800000",
  "Thalassobaculales" = "#a9a9a9",
  "Verrucomicrobiales" = "turquoise1",
  "Xanthomonadales" = "orange"
  )
  
# sort and rename factor levels
order_sums_water$treatment <- factor(order_sums_water$treatment, 
  levels = c("glyph", "control"))
levels(order_sums_water$treatment) <- c("Treatment", "Control")  
order_sums_water$nucleic_acid <- factor(order_sums_water$nucleic_acid, 
  levels = c("dna", "cdna"))
levels(order_sums_water$nucleic_acid) <- c("16S rRNA gene", "16S rRNA") 

# plot an array of 4 geom_areas
ggplot(order_sums_water, aes(x = new_day, y = Abundance2, fill = order)) +
  geom_area(stat = "identity") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1.2) +
  scale_fill_manual(breaks = levels(order_sums_water$order), values = fill_values_water) +
  guides(colour = FALSE, size = FALSE, width = FALSE, fill = guide_legend(ncol = 1,
    keyheight = 1.5, label.theme = element_text(size = 15, face = "italic",
    angle = 0), (title = NULL))) +
  scale_x_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size = 17),
    axis.title = element_text(size = 20, face = "bold"),
    legend.background = element_rect(fill = "grey90", linetype = "solid"),
    panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    strip.text.x = element_text(size = 15, face = "bold")) +
  labs(x = "Days", y = "Relative abundance [%]") +
  facet_wrap(~ treatment + nucleic_acid, nrow = 2)
  
  
# get data per OTU, setting threshold for samples and clusters
community_subset_biofilm <- droplevels(subset(mothur_ra_melt_mean, days > 40
  & Abundance > 0.15 & habitat == "biofilm"))

community_subset_biofilm %>%
  group_by(new_day, order, nucleic_acid, treatment) %>%
  summarise(Abundance2 = sum(Abundance)) %>%
  ungroup() %>%
  # Replace missing values by 0
  spread(key = order, value = Abundance2) %>%
  replace(., is.na(.), 0) %>%
  gather(key = order, value = Abundance2, -c(new_day, nucleic_acid, treatment)) -> 
    order_sums_biofilm
    
# recycle values from water plot
order_sums_biofilm$order <- factor(order_sums_biofilm$order, 
  levels = c(
    # alphaproteos
    "Caulobacterales", #
    "Rhizobiales", #
    "Rhodobacterales", #
    "Rhodospirillales", #
    "Sneathiellales", #
    "Sphingomonadales", #
    "Parvibaculales",
    "Thalassobaculales",
    # gammaproteos
    "Aeromonadales", #
    "Alteromonadales", #
    "Betaproteobacteriales", #
    "Gammaproteobacteria_Incertae_Sedis",
    "Oceanospirillales", #
    "Pseudomonadales", #
    "Xanthomonadales", #
    # Actinobacteria
    "Corynebacteriales",
    # Bacteroidia
    "Bacillales", #
    "Bacteroidia_unclassified",
    "Chitinophagales", #
    "Flavobacteriales", #
    "Sphingobacteriales",
    # Planctomycetacia
    "Planctomycetales", #
    "OM190_or", #
    # Verrucomicrobia
    "Verrucomicrobiales" #
    ))
    
fill_values_biofilm <- fill_values_water

# sort and rename factor levels
order_sums_biofilm$treatment <- factor(order_sums_biofilm$treatment, 
  levels = c("glyph", "control"))
levels(order_sums_biofilm$treatment) <- c("Treatment", "Control")  
order_sums_biofilm$nucleic_acid <- factor(order_sums_biofilm$nucleic_acid, 
  levels = c("dna", "cdna"))
levels(order_sums_biofilm$nucleic_acid) <- c("16S rRNA gene", "16S rRNA") 

# biofilm are plot for SI
ggplot(order_sums_biofilm, aes(x = new_day, y = Abundance2, fill = order)) +
  geom_area(stat = "identity") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1.2) +
  scale_fill_manual(breaks = levels(order_sums_biofilm$order), values = fill_values_biofilm) +
  guides(colour = FALSE, size = FALSE, width = FALSE, fill = guide_legend(ncol = 1,
    keyheight = 1.5, label.theme = element_text(size = 15, face = "italic",
    angle = 0), (title = NULL))) +
  scale_x_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text = element_text(size = 17),
    axis.title = element_text(size = 20, face = "bold"),
    legend.background = element_rect(fill = "grey90", linetype = "solid"),
    panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    strip.text.x = element_text(size = 15, face = "bold")) +
  labs(x = "Days", y = "Relative abundance [%]") +
  facet_wrap(~ treatment + nucleic_acid, nrow = 2)