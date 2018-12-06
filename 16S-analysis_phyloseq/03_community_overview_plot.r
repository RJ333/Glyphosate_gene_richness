require(ggplot2)
require(scales)

# get data per OTU, setting threshold for samples and clusters
sample_subset <- droplevels(subset(mothur_ra_melt_mean, days > 40 
							& Abundance > 0.15
							& habitat == "water" 
							& treatment == "glyph"))
# check required number of colours per order and number of classes
length(levels(droplevels(sample_subset$class)))
length(levels(droplevels(sample_subset$order))) 
# sort orders based on class for plotting
sample_subset$order <- factor(sample_subset$order, 
									   # alphaproteos
							levels = c("Caulobacterales",
									   "Rhizobiales",
									   "Rhodobacterales",
									   "Rhodospirillales",
									   "Sneathiellales",
									   "Sphingomonadales",
									   "Parvibaculales",
									   "Thalassobaculales",
									   # gammaproteos
									   "Alteromonadales",
									   "Betaproteobacteriales",
									   "Pseudomonadales",
									   # bacteroidetes/
									   "Chitinophagales",
									   "Sphingobacteriales",
									   "Flavobacteriales"))
# assign specific colour to make plot distuingishable
fill_values <- c("Alteromonadales" = "orange",
					"Betaproteobacteriales" = "pink",
					"Caulobacterales" = "black",
					"Chitinophagales" = "purple",
					"Flavobacteriales" = "green",
					"Sneathiellales" = "white",
					"Parvibaculales" = "green3",
					"Pseudomonadales" = "grey30",
					"Rhizobiales" = "red",
					"Rhodobacterales" = "lightblue",
					"Rhodospirillales" = "yellow",
					"Sphingobacteriales" = "darkred",
					"Sphingomonadales" = "grey",
					"Thalassobaculales" = "blue2")	

#plotting all selected clusters in bar plot ordered by class 
# and displaying orders over time for DNA and RNA
ggplot(sample_subset, aes(x = new_day, group = order)) +
	scale_fill_manual(breaks = levels(sample_subset$order), 
				      values = fill_values) +
	geom_bar(data = subset(sample_subset, 
						   nucleic_acid == "dna" & 
						   treatment == "glyph"),
			 aes(x = new_day - 0.5, 
				 y = Abundance), 
			 fill = "black", 
		     width = 0.9, 
		     stat = "sum") +
	geom_bar(data = subset(sample_subset, 
						   nucleic_acid == "dna" & 
						   treatment == "glyph"),
			 aes(x = new_day - 0.5, 
				 y = Abundance, 
				 fill = order), 
			 width = 0.6, 
			 stat = "identity") +
	geom_bar(data = subset(sample_subset, 
						   nucleic_acid == "cdna" & 
						   treatment == "glyph"),
			 aes(x = new_day + 0.5, 
				 y = Abundance), 
			 fill = "black", 
			 width = 0.9, 
			 stat = "sum") +
	geom_bar(data = subset(sample_subset, 
						   nucleic_acid == "cdna" & 
						   treatment == "glyph"),
			 aes(x = new_day + 0.5, 
				 y = Abundance, 
				 fill = order), 
			 width = 0.6, 
			 stat = "identity") +
	geom_vline(data = subset(sample_subset, 
							 treatment == "glyph"),
			   aes(xintercept = 1.5),
			   linetype = "dashed", size = 1.2) +
	guides(colour = FALSE, 
		   size = FALSE, 
		   width = FALSE,
		   fill = guide_legend(ncol = 1,
							   keyheight = 1.5,
							   label.theme = element_text(size = 15,
														  face = "italic",
														  angle = 0),
											(title = NULL))) +
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
	scale_y_continuous(expand = c(0,0)) +
	theme_bw() +
	theme(axis.text = element_text(size = 17)) +
	theme(axis.title = element_text(size = 20,
									face = "bold")) +
	theme(legend.background = element_rect(fill = "grey90", 
										   linetype = "solid")) +
	theme(panel.grid.major = element_line(colour = NA, 
										  size = 0.2)) +
	theme(panel.grid.minor = element_line(colour = NA, 
										  size = 0.5)) +
	labs(x = "Days", 
		 y = "Relative abundance [%]") +
  annotate("text", 
		   x = -27, 
		   y = 90, 
		   label = "a)", 
		   color = "black", 
		   size = 6, 
		   angle = 0, 
		   fontface = "bold") +
  annotate("text", 
		   x = -22.5, 
		   y = 90, 
		   label = "b)", 
		   color = "black", 
		   size = 6, 
		   angle = 0, 
		   fontface = "bold")

ggsave(file = "fig_02_relative_community_order_level_sorted_for_manuscript.png", 
	   width = 16, 
	   height = 8)