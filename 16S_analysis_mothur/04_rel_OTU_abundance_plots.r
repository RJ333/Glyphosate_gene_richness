# this interactive script subsets a specific OTU and plots it relative
# abundance in water column and biofilm over days per nucleic acid and treatment
# the samples include parallels
# Input is a custom list of OTUs 


# first a vector of OTUs to be plotted with abundance threshold
OTU_list <- subset(aggregate(Abundance ~ OTU, 
							 data = mothur_ra_melt, 
							 max), 
				   Abundance > 0.05)[,1]
				   
# define subset function 
get_current_otu_data <- function(x) {
	subset(mothur_ra_melt, OTU == x)
}
# where the plots should be stored
plot_folder <- "/data/projects/glyphosate/plots/R/OTU_abundance/"

# run a for loop to ggplot each OTU in list with own title and file name
for (i in OTU_list){
current_otu_data <- get_current_otu_data(i)
print(paste("OTU is", i))

species_title <- unique(paste(current_otu_data$family, 
							  current_otu_data$genus, 
							  current_otu_data$OTU, 
							  sep = "_"))

current_otu_data$treatment2 <- factor(current_otu_data$treatment, 
											  labels = c("Control", "Treatment"))
		
current_plot <- ggplot(data = current_otu_data, 
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
	theme_bw() +
	ggtitle(species_title) +
	theme(axis.text = element_text(size = 18))+
	theme(panel.grid.major = element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor = element_line(colour = NA, size = 0.5))+
	#theme(legend.position = "none")+
	theme(axis.title = element_blank()) +
	facet_wrap(~ habitat, scales = "free")
	ggsave(current_plot, file = paste(plot_folder, 
									  species_title, 
									  ".png", 
									  sep = ""), 
						 width = 13, 
						 height = 7)
}