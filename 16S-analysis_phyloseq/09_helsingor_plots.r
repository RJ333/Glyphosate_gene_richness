# Set the working dir with mothur files in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("mothur_glyph_002.RData")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
library(gridExtra)

mothur_ps_bla <- filter_taxa(mothur_ps2, function (x) {sum(x > 2) >= 1}, prune = TRUE)
mothur_ps2_ra <- transform_sample_counts(mothur_ps_bla, function(x){(x / sum(x)) * 100})

mothur_water_dna <- subset_samples(mothur_ps2_ra, habitat == "water" &
												  nucleic_acid == "dna" &
												  days > 43,
												  prune = TRUE)

sample_data(mothur_water_dna)$days <- as.factor(sample_data(mothur_water_dna)$days)					   
sample_data(mothur_water_dna)$new_day <- as.factor(sample_data(mothur_water_dna)$new_day)
											  
water_ordi <- ordinate(mothur_water_dna, "NMDS", "bray", try = 100, autotransform = TRUE)

plot_ordination(mothur_water_dna, water_ordi, type = "samples", shape = "treatment", color = "new_day") + 
	 geom_polygon(aes(fill = disturbance), alpha = 0.5, size = 0.01) + 
	 geom_point(size=5, colour = "grey60") +
	 ggtitle("Water column, DNA") +
	 geom_text(aes(label = new_day), colour = "black", size = 4) +
	 guides(color = FALSE) +
	 labs(shape = "Treatment", fill = "Disturbance") +
	 theme_bw()+
	 theme(plot.title = element_text(color = "black", size = 20),
		   axis.text = element_text(size = 18),
		   axis.title = element_text(size = 20, face = "bold"),
		   legend.title = element_text(size = 15, face = "bold"), 
		   legend.text = element_text(size = 13),
		   panel.grid.major = element_line(colour = NA, size = 0.2),
		   panel.grid.minor = element_line(colour = NA, size = 0.5))



		   
sample_data(mothur_ps2)$condition_diversity[sample_data(mothur_ps2)$condition_diversity == "start"] <- "untreated"
# define order of factor levels
condition_order <- c("untreated", "treated", "22 to 36", "43 to 71")

waterdnaglyph <- subset_samples(mothur_ps2, habitat == "water" & 
									   nucleic_acid == "dna" &
									   treatment == "glyph" &
									   days > 43,
									   prune = TRUE)
									   
richness_plot_dna_water_control <- plot_richness(waterdnacontrol, x = "condition_diversity", 
												 color = "as.factor(new_day)",
												 measures = c("Shannon" )) + 
									 geom_point(alpha = 0.8, size = 4) +
									 coord_cartesian(ylim = c(1.8, 2.7)) +
									 ggtitle("Shannon index of water column, DNA, control") +
									 labs(x = "Time period") +
									 guides(color = FALSE) +
									 theme_bw()+
									 theme(plot.title = element_text(color = "black", size = 20),
										 axis.text = element_text(size = 18),
										 axis.title = element_text(size = 20, face = "bold"),
										 legend.title = element_text(size = 15, face = "bold"), 
										 legend.text = element_text(size = 13),
										 panel.grid.major = element_line(colour = NA, size = 0.2),
										 panel.grid.minor = element_line(colour = NA, size = 0.5))

richness_plot_dna_water_control$data$condition_diversity <- as.character(richness_plot_dna_water_control$data$condition_diversity)
richness_plot_dna_water_control$data$condition_diversity <- factor(richness_plot_dna_water_control$data$condition_diversity, levels=condition_order)
print(richness_plot_dna_water_control)

plots_div <- list()
plots_div <- list(richness_plot_dna_water_glyph, richness_plot_dna_water_control)							  

div_plot_object <- do.call("grid.arrange", c(plots_div, nrow = 1))
