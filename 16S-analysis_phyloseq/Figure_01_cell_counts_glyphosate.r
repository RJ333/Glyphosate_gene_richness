#!/usr/bin/env Rscript
library(scales)
library(ggplot2)

setwd("/data/projects/glyphosate/reads/mothur_processed/")

# set path for plots, create directory if not existant
plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

cell_counts_glyph <- read.csv("cell_counts_glyph.csv", sep = ";")
cell_counts_glyph_0 <- subset(cell_counts_glyph, new_day >= -7)

cell_counts_glyph_plot <- ggplot(cell_counts_glyph_0, aes(x = new_day,
  colour = treatment, linetype = treatment, group = treatment)) +
  geom_errorbar(aes(ymin = cells_ml - cells_se, ymax = cells_ml + cells_se),
    linetype = 1, width = 2, size = 1, alpha = 0.7) +
  geom_errorbar(aes(ymin = (glyph_micromol - glyph_se_micromol) * 400000,
    ymax = (glyph_micromol + glyph_se_micromol) * 400000), linetype = 1, 
    width = 2, size = 0.7, alpha = 0.8) +
  geom_line(aes(y = cells_ml, group = treatment, colour = treatment), 
    linetype = 1, size = 1, alpha = 0.8) +
  geom_point(aes(y = glyph_theor_micromol * 400000, shape = "glyph"), 
    alpha = 0.5, size = 3) +
  geom_point(aes(y = glyph_micromol * 400000, shape = "glyph_deg"), 
    alpha = 0.7, size = 2) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1.5, alpha = 0.5) +
  scale_y_continuous(label =  function(x) {ifelse(x == 0, "0", parse(text = 
    gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))},
    sec.axis = sec_axis(~ . / 400000, name = expression(bold(paste(
      "Glyphosate  ", bgroup("[",µM,"]")))))) +
  scale_colour_manual(values = c("control" = "grey70", "glyph" = "black"),
    name = "Microcosm", breaks = c("control", "glyph"), 
    labels = c("Control", "Treatment")) +
  scale_shape_manual(values = c("glyph" = 2, "glyph_deg" = 17), 
    name = "Glyphosate decrease by", breaks = c("glyph", "glyph_deg"),
    labels = c("calculated\ndilution", "measured\nconcentration")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(x = "Days", y = expression(bold(paste("Total cell count  ",
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
    theme(legend.position = "bottom", legend.direction = "vertical")

ggsave(cell_counts_glyph_plot, file = paste(plot_path, 
  "Figure_1_cellcounts_glyph.pdf", sep = ""),
  device = "pdf", width = 18, height = 14, dpi = 300, unit = "cm")