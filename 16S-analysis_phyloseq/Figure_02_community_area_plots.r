#!/usr/bin/env Rscript
#library(scales)
library(ggplot2)
library(phyloseq)
library(tidyverse)

setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("glyphosate_mothur_in_phyloseq.RData")

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
water_areas <- ggplot(order_sums_water, aes(x = new_day, y = Abundance2, fill = order)) +
  geom_area(stat = "identity") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1) +
  scale_fill_manual(breaks = levels(order_sums_water$order), values = fill_values_water) +
  guides(colour = FALSE, size = FALSE, width = FALSE, fill = guide_legend(ncol = 1,
    keyheight = 1.2, label.theme = element_text(size = 12, face = "italic",
    angle = 0), (title = NULL))) +
  scale_x_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(angle = 90, vjust = 1),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 14, colour = "black", face = "bold"),
    legend.background = element_rect(fill = "grey90", linetype = "solid")) +
  labs(x = "Days", y = "Relative abundance [%]") +
 # theme(legend.position = "bottom", legend.direction = "horizontal") +
  facet_wrap(~ treatment + nucleic_acid, nrow = 2)
  
ggsave(water_areas, file = paste(plot_path, "Figure_2_water_communities.pdf", sep = ""),
  device = "pdf", width = 26.0, height = 18, dpi = 300, unit = "cm")
  
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
biofilm_areas <- ggplot(order_sums_biofilm, aes(x = new_day, y = Abundance2, fill = order)) +
  geom_area(stat = "identity") +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1) +
  scale_fill_manual(breaks = levels(order_sums_biofilm$order), values = fill_values_biofilm) +
  guides(color = FALSE) +
  guides(size = FALSE) +
  guides(width = FALSE) +
  guides(fill = guide_legend(label.theme = element_text(size = 12, face = "italic",
    angle = 0), ncol = 1, keyheight = 1.2, (title = NULL))) +
  scale_x_continuous(expand = c(0, 0), breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(angle = 90, vjust = 1),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 14, colour = "black", face = "bold"),
    legend.background = element_rect(fill = "grey90", linetype = "solid")) +
  labs(x = "Days", y = "Relative abundance [%]") +
  facet_wrap(~ treatment + nucleic_acid, nrow = 2)
  
ggsave(biofilm_areas, file = paste(plot_path, "SI_4_biofilm_communities.pdf", sep = ""),
  device = "pdf", width = 26.0, height = 18, dpi = 300, unit = "cm")