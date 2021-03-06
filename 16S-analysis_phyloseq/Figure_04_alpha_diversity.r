#!/usr/bin/env Rscript
library(ggplot2)
library(phyloseq)

setwd("/data/projects/glyphosate/reads/mothur_processed/")

plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

load("glyphosate_mothur_in_phyloseq.RData")

erich_mothur <- estimate_richness(mothur_full, measures = c(
  "Observed",
  "Chao1",
  "ACE",
  "Shannon",
  "Simpson",
  "InvSimpson",
  "Fisher"))
erich_mothur_meta <- cbind(erich_mothur, sample_data(mothur_full)[,c(1:7)])

# reorder and rename factor levels for plotting, adjust displayed label names
erich_mothur_meta$habitat <- relevel(erich_mothur_meta$habitat, "water")
erich_mothur_meta$nucleic_acid <- relevel(erich_mothur_meta$nucleic_acid, "dna")
labels_nucleic_acid <- c("16S rRNA gene", "16S rRNA")
labels_habitat <- c("Free-living", "Biofilm")
levels(erich_mothur_meta$habitat) <- labels_habitat
levels(erich_mothur_meta$nucleic_acid) <- labels_nucleic_acid

# plotting Alpha diversity
shannon_plot <- ggplot(erich_mothur_meta, aes(x = new_day, y = Shannon, colour = treatment)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1, alpha = 0.5) +
  stat_summary(aes(colour = treatment), fun.y = "mean", geom = "line", 
    alpha = 0.7, size = 1) +
  scale_colour_manual(values = c("glyph" = "black", "control" = "grey50"),
	name = "Microcosm  ", breaks = c("glyph", "control"), 
    labels = c("Treatment", "Control")) +
  coord_cartesian(ylim = c(1, 3)) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 14, face = "bold")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(x = "Days", y = "Shannon index") +
  facet_wrap(~ habitat + nucleic_acid)

ggsave(shannon_plot, file = paste(plot_path, "Figure_4_Shannon_DNA_RNA.pdf", sep = ""),
  device = "pdf", width = 18.0, height = 15, dpi = 300, unit = "cm")