#!/usr/bin/env Rscript
library(ggplot2)
library(phyloseq)
library(gridExtra)

setwd("/data/projects/glyphosate/reads/mothur_processed/")

plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

load("glyphosate_mothur_in_phyloseq.RData")

# Supplement 6: OTU abundance plots

# define subset function
get_current_otu_data <- function(x) {
	subset(mothur_ra_melt, OTU == x)
}

# list of OTUs mentioned in paper and supplement
OTU_list <- c(
  "Otu000001", # SI 6 a) most abundant Rhizobiales
  "Otu000006", # Pseudomonas control increasing
  "Otu000007", # SI 6 b)
  "Otu000009",
  "Otu000011", # SI 6 c)
  "Otu000018",
  "Otu000019",
  "Otu000023", # SI 6 g) 1. Methylotenera
  "Otu000024",
  "Otu000025", # SI 6 f)
  "Otu000028",
  "Otu000029", # SI 6 k)
  "Otu000032",
  "Otu000034",
  "Otu000035",
  "Otu000036", # SI 6 i)
  "Otu000037",
  "Otu000038", # SI 6 n) Pseudolabrys
  "Otu000042", # SI 6 l) Brevundimonas
  "Otu000043",
  "Otu000044", # 2. Methylotenera 
  "Otu000046", # SI 6 h)
  "Otu000049",
  "Otu000050",
  "Otu000056",
  "Otu000058",
  "Otu000059", # SI 6 p)
  "Otu000069",
  "Otu000070",
  "Otu000072",
  "Otu000078", # SI 6 j)
  "Otu000086",
  "Otu000094",
  "Otu000098", # SI 6 m) Defluviimonas
  "Otu000109", # SI 6 d)
  "Otu000129", # SI 6 e)
  "Otu000139",
  "Otu000176",
  "Otu000191",
  "Otu000320",
  "Otu000098",
  "Otu000042"
  )

# run a for loop to store each OTU in list with own title and file name
current_plot_list <- list()
strip_text_habitat <- c("Biofilm", "Free-living")

for (i in OTU_list) {
    current_otu_data <- get_current_otu_data(i)
    print(paste("OTU is", i))

species_title <- unique(paste(current_otu_data$family, current_otu_data$genus,
  current_otu_data$OTU, sep = " "))

levels(current_otu_data$habitat) <- strip_text_habitat

current_plot_list[[i]] <- ggplot(data = current_otu_data, aes(x = days - 69,
  y = Abundance, group = nucleic_acid, lty = nucleic_acid)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1) +
  geom_point(data = subset(current_otu_data, treatment == "control"),
    aes(colour = treatment), alpha = 1, size = 1) +
  geom_point(data = subset(current_otu_data, treatment == "glyph"),
    aes(colour = treatment), size = 1) +
  stat_summary(data = subset(current_otu_data, treatment == "glyph"),
    aes(colour = treatment), fun.y = "mean", geom = "line", size = 1.5) +
  stat_summary(data = subset(current_otu_data, treatment == "control"),
    aes(colour = treatment), fun.y = "mean", geom = "line", size = 1.5, alpha = 1) +
  scale_linetype_manual(values = c("dna" = 1, "cdna" = 6), 
    name = "Nucleic acid  ", breaks = c("cdna", "dna"), labels = c("16S rRNA", "16S rRNA gene")) +
  scale_colour_manual(values = c("glyph" = "black", "control" = "grey50"),
    name = "Microcosm  ", breaks = c("glyph", "control"), labels = c("Treatment",
    "Control")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  ggtitle(species_title) +
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
}

pdf(paste(plot_path, "allplots.pdf", sep = ""), onefile = TRUE, height = 6, width = 9)
for (plott in seq_along(current_plot_list)) {
  print(current_plot_list[[plott]])  
}
dev.off()