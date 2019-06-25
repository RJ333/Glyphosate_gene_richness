#!/usr/bin/env Rscript
library(ggplot2)
library(phyloseq)

setwd("/data/projects/glyphosate/reads/mothur_processed/")

plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

load("glyphosate_mothur_in_phyloseq.RData")

# Figure 5 and Supplement 5: OTU abundance plots
# define subset function
get_current_otu_data <- function(x) {
	subset(mothur_ra_melt, OTU == x)
}

# list of OTUs mentioned in paper and supplement
OTU_list <- c(
  #DESeq2 OTUs
  # "Otu000007",
  # "Otu000011",
  # "Otu000018",
  # "Otu000025",
  # "Otu000032",
  # "Otu000036",
  # "Otu000037",
  # "Otu000038",
  # "Otu000023",
  # "Otu000046",
  # "Otu000049",
  # "Otu000056",
  # "Otu000058",
  # "Otu000059",
  # "Otu000070",
  # "Otu000072",
  # "Otu000078",
  # "Otu000094",
  # "Otu000109",
  # "Otu000129",
  # "Otu000139",
  # "Otu000176",
  # "Otu000191",
  # "Otu000320",
  # "Otu000098",
  # "Otu000042",
  # 2. Methylotenera 
  "Otu000044",
  # Pseudomonas control increasing
  "Otu000006",
  # most abundant Rhizobiales
  "Otu000001")

Pseudomonas_OTUs <- c(
  # outcommented included above
  #"Otu000006",
  "Otu000007",
  "Otu000009",
  "Otu000019",
  "Otu000024",
  "Otu000028",
  "Otu000029",
  "Otu000034",
  "Otu000035",
  #"Otu000036",
  "Otu000043",
  "Otu000050",
  "Otu000069",
  #"Otu000078",
  "Otu000086"
  )

strip_text_habitat <- c("Biofilm", "Free-living")

# run a for loop to plot each OTU in list with own title and file name
for (i in OTU_list) {
    current_otu_data <- get_current_otu_data(i)
    print(paste("OTU is", i))

species_title <- unique(paste(current_otu_data$family, current_otu_data$genus,
  current_otu_data$OTU, sep = " "))

levels(current_otu_data$habitat) <- strip_text_habitat

current_plot <- ggplot(data = current_otu_data, aes(x = days - 69,
  y = Abundance, group = nucleic_acid, lty = nucleic_acid)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1.2) +
  geom_point(data = subset(current_otu_data, treatment == "control"),
    aes(colour = treatment), alpha = 1) +
  stat_summary(data = subset(current_otu_data, treatment == "control"),
    aes(colour = treatment), fun.y = "mean", geom = "line", size = 2, alpha = 1) +
  stat_summary(data = subset(current_otu_data, treatment == "glyph"),
    aes(colour = treatment), fun.y = "mean", geom = "line", size = 2) +
  geom_point(data = subset(current_otu_data, treatment == "glyph"),
    aes(colour = treatment)) +
  scale_linetype_manual(values = c("dna" = 1, "cdna" = 6), name = "Nucleic acid  ",
    breaks = c("cdna", "dna"), labels = c("16S rRNA", "16S rRNA gene")) +
  scale_colour_manual(values = c("glyph" = "black", "control" = "grey50"),
    name = "Microcosm  ", breaks = c("glyph", "control"), labels = c("Treatment",
    "Control")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  ggtitle(species_title) +
  theme(axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    strip.text.x = element_text(size = 15, face = "bold")) +
  labs(x = "Days", y = "Relative abundance [%]") +
  facet_wrap(~ habitat, scales = "free")
  # ggsave(current_plot, file = paste(plot_path, species_title,".png", sep = ""),
    # width = 13, height = 7)
  print(current_plot)
}
