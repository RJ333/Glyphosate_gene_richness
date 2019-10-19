#!/usr/bin/env Rscript
library(scales)
library(ggplot2)
library(cowplot)

setwd("/data/projects/glyphosate/reads/mothur_processed/")

# set path for plots, create directory if not existant
plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

cell_counts_glyph <- read.csv("cell_counts_glyph.csv", sep = ";")
cell_counts_glyph_0 <- subset(cell_counts_glyph, new_day >= -7)

# Figure 1 A: cell counts
cell_counts_glyph_plot <- ggplot(cell_counts_glyph_0, aes(x = new_day,
  colour = treatment, linetype = treatment, group = treatment)) +
  geom_errorbar(aes(ymin = cells_ml - cells_se, ymax = cells_ml + cells_se),
    linetype = 1, width = 2, size = 1, alpha = 0.7) +
  geom_line(aes(y = cells_ml, group = treatment, colour = treatment), 
    linetype = 1, size = 1, alpha = 0.8) + 
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1.5, alpha = 0.5) +
  scale_colour_manual(values = c("control" = "grey70", "glyph" = "black"),
    name = "Microcosm", breaks = c("control", "glyph"), 
    labels = c("Control", "Treatment"),
      guide = guide_legend(override.aes = list(
        colour = c("grey70", "black"),
        shape = c(NA, NA),
        linetype = c("solid", "solid")))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(x = "Days", y = expression(bold(paste("Total cell counts  ",
    bgroup("[",cells~mL^{-1},"]"))))) +
  theme_bw() +
    theme(panel.grid.major = element_line(colour = NA, size = 0.2),
      panel.grid.minor = element_line(colour = NA, size = 0.5),
      axis.title = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(angle = 90, vjust = 1, face = "bold"),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      strip.text.x = element_text(size = 16, colour = "black", face = "bold")) +
    #theme(legend.position = "bottom", legend.direction = "horizontal")
    theme(legend.position = "none") 
  
# Figure 1 B: glyphosate and AMPA concentrations
glyphosate_plot <- ggplot(cell_counts_glyph_0, aes(x = new_day,
  colour = treatment, linetype = treatment, group = treatment)) +
  geom_errorbar(aes(ymin = (glyph_micromol - glyph_se_micromol),
    ymax = (glyph_micromol + glyph_se_micromol)), linetype = 1, 
    width = 2, size = 0.7, alpha = 0.8) +
  geom_point(aes(y = glyph_theor_micromol, shape = "glyph"), 
    alpha = 0.5, size = 3) +
  geom_point(aes(y = glyph_micromol, shape = "glyph_deg"), 
    alpha = 0.7, size = 2) +
  geom_point(aes(y = ampa_micromol * 340, shape = "AMPA"), alpha = 1, size = 3) +
  #geom_point(aes(y = sarc_ala_micromol, shape = "sarc"), alpha = 1, size = 4) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1.5, alpha = 0.5) +
  scale_y_continuous(sec.axis = sec_axis(~ . / 340, name = expression(bold(paste(
      "AMPA  ", bgroup("[",µM,"]")))))) +
  scale_colour_manual(values = c("control" = "grey70", "glyph" = "black"),
    name = "Microcosm", breaks = c("control", "glyph"), 
    labels = c("Control", "Treatment"),
      guide = guide_legend(override.aes = list(
        colour = c("grey70", "black"),
        shape = c(NA, NA),
        linetype = c("solid", "solid")))) +
  scale_shape_manual(values = c("glyph" = 2, "glyph_deg" = 17, "AMPA" = 1, "sarc" = 13), 
    name = "Parameter", breaks = c("glyph", "glyph_deg", "AMPA", "sarc"),
    labels = c("Glyphosate calculated dilution", "Glyphosate measured concentration", "AMPA measured concentration", "Sarcosine/L-Alanine"),
      guide = guide_legend(override.aes = list(
        shape = c(2, 17, 1),
        alpha = c(0.4, NA, NA),
        size = c(3, 3, 3)))
        ) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(x = "Days", y = expression(bold(paste("Glyphosate  ",
    bgroup("[",µM,"]"))))) +
  theme_bw() +
    theme(panel.grid.major = element_line(colour = NA, size = 0.2),
      panel.grid.minor = element_line(colour = NA, size = 0.5),
      axis.title = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(angle = 90, vjust = 1, face = "bold"),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      strip.text.x = element_text(size = 16, colour = "black", face = "bold")) +
    theme(legend.position = "bottom", legend.direction = "vertical")

# save plot, add asterisk and labels A/B in paint
Figure_1_combined <- plot_grid(cell_counts_glyph_plot, glyphosate_plot, align = "v", nrow = 2, rel_heights = c(4/10, 6/10))
save_plot(filename = paste(plot_path, "Figure_1_combined.tiff", sep = ""), 
  Figure_1_combined, base_height = 10, base_width = 7.0866141732, dpi = 450)

# scp -r -i /drives/d/ssh/denbi.key centos@193.196.20.93:/data/projects/glyphosate/reads/mothur_processed/plots/Figure_1_combined.pdf /mnt/d/denbi/chandler/cell_counts/Figure_1_combined.pdf

