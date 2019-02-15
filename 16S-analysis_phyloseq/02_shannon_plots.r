# Set the working dir with mothur files in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("mothur_glyph_002.RData")

library("ggplot2")
library("gridExtra")
library("phyloseq")

# include singletons
mothur_ps2
mothur_div <- mothur_ps2
									   
erich_mothur <- estimate_richness(mothur_div, measures = c("Observed", 
																	 "Chao1", 
																	 "ACE", 
																	 "Shannon", 
																	 "Simpson", 
																	 "InvSimpson", 
																	 "Fisher"))
																	 
erich_mothur_meta <- cbind(erich_mothur, sample_data(mothur_div)[,c(1:7)])

# reorder and rename factor levels for plot
erich_mothur_meta$habitat <- relevel(erich_mothur_meta$habitat, "water")
erich_mothur_meta$nucleic_acid <- relevel(erich_mothur_meta$nucleic_acid, "dna")
labs_nucleic_acid <- c("DNA", "RNA")
labs_habitat <- c("Water column", "Biofilm")
levels(erich_mothur_meta$habitat) <- labs_habitat
levels(erich_mothur_meta$nucleic_acid) <-labs_nucleic_acid
	
shannon_plot <- ggplot(erich_mothur_meta, aes(x = new_day, y = Shannon, colour = treatment)) + 
	geom_point(alpha = 0.8, size = 4) +
	geom_vline(aes(xintercept = 1), 
			   linetype = "dashed", 
			   size = 1.2) +
	stat_summary(aes(colour = treatment), 
				 fun.y = "mean",  
				 geom = "line",
				 alpha = 0.75,
				 size = 2) +
	scale_colour_manual(values = c("glyph" = "black", 
								   "control" = "grey50"), 
						name = "Microcosm  ", 
						breaks = c("glyph", 
								   "control"), 
						labels = c("Treatment", 
								   "Control")) +
	#ggtitle("Shannon index") +
	coord_cartesian(ylim = c(1, 3)) +
	theme_bw() +
	theme(axis.text = element_text(size = 18),
		  axis.title = element_text(size = 20, face = "bold"),
		  legend.title = element_text(size = 15, face = "bold"), 
		  legend.text = element_text(size = 13),
		  panel.grid.major = element_line(colour = NA, size = 0.2),
		  panel.grid.minor = element_line(colour = NA, size = 0.5),
		  #strip.background = element_blank(),
		  strip.text.x = element_text(size = 15, face = "bold")
		  ) +
  	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	labs(x = "Days", y = "Shannon index") +
	#theme(legend.position = "none")+
	#theme(axis.title = element_blank()) +
	facet_wrap(~ habitat + nucleic_acid)

plot_folder <- "/data/projects/glyphosate/plots/R/diversity/"

ggsave(shannon_plot, file = paste(plot_folder, "Shannon_DNA_RNA.png", 
								  sep = ""),
								  height = 10,
								  width = 14)