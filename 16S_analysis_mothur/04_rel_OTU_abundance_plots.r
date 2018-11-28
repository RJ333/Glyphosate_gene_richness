# this interactive script subsets a specific OTU and plots it relative
# abundance in water column and biofilm over days per nucleic acid and treatment
# the samples include parallels
# Input is a custom list of OTUs 


# first a vector of OTUs to be plotted with abundance threshold
# OTU_list <- subset(aggregate(Abundance ~ OTU, 
							 # data = mothur_ra_melt, 
							 # max), 
				   # Abundance > 0.05)[,1]
				

OTU_list <- c("Otu000001", "Otu000003", "Otu000004", "Otu000007", "Otu000008",
				"Otu000009", "Otu000010", "Otu000011", "Otu000012", "Otu000013",
				"Otu000014", "Otu000015", "Otu000025", "Otu000030", "Otu000032",
				"Otu000034", "Otu000036", "Otu000037", "Otu000041", "Otu000044",
				"Otu000046", "Otu000049", "Otu000050", "Otu000056", "Otu000058",
				"Otu000059", "Otu000070", "Otu000072", "Otu000078", "Otu000094",
				"Otu000228", "Otu000401", "Otu000038", "Otu000047", "Otu000048",
				"Otu000051", "Otu000062", "Otu000065", "Otu000066", "Otu000081",
				"Otu000098", "Otu000121", "Otu000135", "Otu000180", "Otu000204",
				"Otu000005", "Otu000006", "Otu000017", "Otu000018", "Otu000019",
				"Otu000020", "Otu000109", "Otu000112", "Otu000129", "Otu000139",
				"Otu000176", "Otu000191", "Otu000214", "Otu000320", "Otu000016",
				"Otu000033", "Otu000042", "Otu000002", "Otu000021", "Otu000023",
				"Otu000054", "Otu000064", "Otu000087", "Otu000123", "Otu000272",
				"Otu000039", "Otu000111", "Otu000210", "Otu000024", "Otu000028",
				"Otu000096", "Otu000097", "Otu000103", "Otu000149", "Otu000181",
				"Otu000186")
			
				
# define subset function for specific phyloseq-object
# get_current_otu_data <- function(x) {
	# subset(mothur_ra_melt, OTU == x)
#}

get_current_otu_data <- function(x) {
	subset(deseq_melt, OTU == x)
}

# where the plots should be stored
#plot_folder <- "/data/projects/glyphosate/plots/R/OTU_abundance/"
plot_folder <- "/data/projects/glyphosate/plots/R/deseq/"
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
						   y = Abundance, 
						   group = nucleic_acid, 
						   lty = nucleic_acid)) + 
	geom_vline(aes(xintercept = 0), 
			   linetype = "dashed", 
			   size = 1.2) +
	geom_point(data = subset(current_otu_data, treatment == "control"), 
		       aes(colour = treatment), 
			   alpha = 1) +
	stat_summary(data = subset(current_otu_data, treatment == "control"), 
	             aes(colour = treatment), 
				 fun.y = "mean",  
				 geom = "line", 
				 size = 2, 
				 alpha = 1) +
	stat_summary(data = subset(current_otu_data, treatment == "glyph"), 
	             aes(colour = treatment), 
				 fun.y = "mean",  
				 geom = "line", 
				 size = 2) +
	geom_point(data = subset(current_otu_data, treatment == "glyph"), 
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