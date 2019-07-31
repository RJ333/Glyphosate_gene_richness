#!/usr/bin/env Rscript
library(ggplot2)
library(phyloseq)

setwd("/data/projects/glyphosate/reads/mothur_processed/")

plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

load("glyphosate_mothur_in_phyloseq.RData")

# Figure 5: Gallaecimonas OTU abundance plot
# define subset function

Gallaeci <-	subset(mothur_ra_melt, OTU == "Otu000011")
Gallaeci_title <- unique(paste(Gallaeci$family, Gallaeci$genus,
  Gallaeci$OTU, sep = "_"))
  
strip_text_habitat <- c("Biofilm", "Free-living")
levels(Gallaeci$habitat) <- strip_text_habitat

Gallaeci_plot <- ggplot(data = Gallaeci, aes(x = days - 69, y = Abundance, 
  group = nucleic_acid, lty = nucleic_acid)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1) +
  geom_point(data = subset(Gallaeci, treatment == "control"),
    aes(colour = treatment), alpha = 1, size = 1) +
  geom_point(data = subset(Gallaeci, treatment == "glyph"),
    aes(colour = treatment), size = 1) +
  stat_summary(data = subset(Gallaeci, treatment == "glyph"),
    aes(colour = treatment), fun.y = "mean", geom = "line", size = 1.5) +
  stat_summary(data = subset(Gallaeci, treatment == "control"),
    aes(colour = treatment), fun.y = "mean", geom = "line", size = 1.5, alpha = 1) +
  scale_linetype_manual(values = c("dna" = 1, "cdna" = 6), 
    name = "Nucleic acid  ", breaks = c("cdna", "dna"), labels = c("16S rRNA", "16S rRNA gene")) +
  scale_colour_manual(values = c("glyph" = "black", "control" = "grey50"),
    name = "Microcosm  ", breaks = c("glyph", "control"), labels = c("Treatment",
    "Control")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  #ggtitle(Gallaeci_title) +
  theme(panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    axis.text = element_text(size = 11.8),
    axis.title = element_text(size = 14, face = "bold"),
    #legend.title = element_text(size = 13, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    strip.text.x = element_text(size = 14, face = "bold")) +
  theme(legend.position = "bottom", legend.direction = "vertical", 
    legend.box = "horizontal") +
  theme(legend.margin = margin(t = 0, unit='cm')) +
  guides(colour = guide_legend(keywidth = 1, keyheight = 0.1, default.unit = "cm")) +
  guides(linetype = guide_legend(keywidth = 1, keyheight = 0.1, default.unit = "cm")) +
  labs(x = "Days", y = "Relative abundance [%]") +
  facet_wrap(~ habitat, scales = "free")
  print(Gallaeci_plot)
  
ggsave(Gallaeci_plot, file = paste(plot_path, "Figure_5_", Gallaeci_title,".pdf", sep = ""),
  device = "pdf", width = 18.0, height = 12, dpi = 300, unit = "cm")

