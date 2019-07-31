#!/usr/bin/env Rscript
library(scales)
library(ggplot2)

setwd("/data/projects/glyphosate/reads/mothur_processed/")

# set path for plots, create directory if not existant
plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

adsorption <- read.csv("glyph_adsorption_data.csv", sep = ";")

strip_text_bottle <- c("Glass", "Polypropylene")
levels(adsorption$bottle) <- strip_text_bottle

ggplot(adsorption, aes(x = hours, y = relative, colour = setup)) + 
  geom_line(alpha = 0.8, size = 1) +
  geom_errorbar(aes(ymin = relative - relative_SD, ymax = relative + relative_SD), 
    size = 1, width = 3, alpha = 0.5, lty = 1) +
  #coord_cartesian(ylim = c(0, 140)) +
  facet_wrap(~ bottle) +
  scale_y_continuous(limits = c(0, 140), breaks = seq(0, 140, by = 20)) +
  scale_colour_manual(values = c("liquid" = "blue", "liquid_sediment" = "black", "liquid_sediment_incubated" = "red"), 
    name = "Bottle setups  ", breaks = c("liquid", "liquid_sediment", "liquid_sediment_incubated"), 
    labels = c("Liquid", "Liquid + sediment", "Incubated liquid + \nsediment")) +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(angle = 90, vjust = 1),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 16, colour = "black", face = "bold")) +
  labs(x = "Hours", y = "Glyphosate in aqueous phase (%)")
  
ggsave(file = paste(plot_path, "SI_1_glyphosate_adsorption.pdf", sep = ""),
  device = "pdf", width = 18.0, height = 12, dpi = 300, unit = "cm")