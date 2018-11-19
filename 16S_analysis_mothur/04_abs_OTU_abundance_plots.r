# this interactive script subsets a specific OTU and plots it absolute
# abundance in the water column over days per nucleic acid and treatment
# absolute abundance is relative abundance * total cell count!
require(scales)
require(ggplot2)

# where the plots should be stored
plot_folder <- "/data/projects/glyphosate/plots/R/OTU_abundance/"

abs_OTU_abundance <- subset(mothur_ra_melt, OTU == "Otu000011" &
										    habitat == "water")
# add taxonomy										   
species_title <- unique(paste(abs_OTU_abundance$family, 
							  abs_OTU_abundance$genus, 
							  abs_OTU_abundance$OTU, 
							  sep = "_"))
							  
ggplot(data = abs_OTU_abundance, 
	   aes(x = days - 69, 
		   y = abs_Abundance, 
	       group = nucleic_acid, 
		   lty = nucleic_acid)) + 
	geom_vline(aes(xintercept = 0), 
			   linetype = "dashed", 
			   size = 1.2) +
	geom_point(data = subset(abs_OTU_abundance, treatment == "control"), 
		       aes(colour = treatment), 
			   alpha = 1) +
	stat_summary(data = subset(abs_OTU_abundance, treatment == "control"), 
	             aes(colour = treatment), 
				 fun.y = "mean",  
				 geom = "line", 
				 size = 2, 
				 alpha = 1) +
	stat_summary(data = subset(abs_OTU_abundance, treatment == "glyph"), 
	             aes(colour = treatment), 
				 fun.y = "mean",  
				 geom = "line", 
				 size = 2) +
	geom_point(data = subset(abs_OTU_abundance, treatment == "glyph"), 
	           aes(colour = treatment)) +
	scale_linetype_manual(values = c("dna" = 1, 
									 "cdna" = 6), 
						name = "Nucleic acid  ", 
						breaks = c("cdna", 
								   "dna"), 
						labels = c("16S rRNA", 
								   "16S rRNA gene")) +
	scale_colour_manual(values = c("glyph" = "black", 
								   "control" = "grey50"), 
						name = "Microcosm  ", 
						breaks = c("glyph", 
								   "control"), 
						labels = c("Treatment", 
								   "Control")) +
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	scale_y_continuous(label= function(x) {ifelse(x == 0, "0", 
										   parse(text = gsub("[+]", "", 
										   gsub("e", " %*% 10^", 
										     scientific_format()(x)))))}) +
	theme_bw() +
	ggtitle(species_title) +
	theme(axis.text = element_text(size = 18))+
	theme(panel.grid.major = element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor = element_line(colour = NA, size = 0.5))+
	#theme(legend.position = "none")+
	theme(axis.title = element_blank()) #+
	ggsave(file = paste(plot_folder, 
									  species_title, 
									  "abs.png", 
									  sep = ""), 
		   width = 13, 
		   height = 7)