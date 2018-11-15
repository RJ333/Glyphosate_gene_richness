require(ggplot2)
require(scales)
############# clean final version ################################################################
#data preparation, setting threshold for samples and clusters which are included in plot
sample_subset<-subset(mothur_ra_melt_mean, days > 40 
						& Abundance > 0.15
						& habitat == "water" 
						& treatment == "glyph")
sample_subset<-droplevels(sample_subset)

length(levels(sample_subset$class)) # number of colors needed for all classes
length(levels(sample_subset$order)) # number of colors needed for all orders


#sort orders based on class, first alpha, then beta, then gammaproteos
sample_subset$order<-factor(sample_subset$order, levels=c("Caulobacterales","Rhizobiales","Rhodobacterales","Rhodospirillales","Sphingomonadales",
"Burkholderiales","Hot Creek 32","Methylophilales","Alteromonadales","Gammaproteobacteria Incertae Sedis","Pseudomonadales","Sphingobacteriales","ARKICE-90","Flavobacteriales"))
#define fill colours for all orders
# fill_values<-c("Flavobacteriales"="green",
	# "Rhizobiales"="red",
	# "Pseudomonadales"="grey30",
	# "Rhodobacterales"="black",
	# "Caulobacterales"="lightblue",
	# "Gammaproteobacteria Incertae Sedis"="white",
	# "Burkholderiales"="pink",
	# "Sphingobacteriales"="grey75",
	# "Hot Creek 32"="grey",
	# "Alteromonadales"="orange",
	# "Sphingomonadales"="purple",
	# "Rhodospirillales"="yellow",
	# "ARKICE-90"="blue2",
	# "Methylophilales"="green3")

fill_values <- c("Alteromonadales" = "orange",
	"Betaproteobacteriales" = "pink",
	"Caulobacterales" = "lightblue",
	"Chitinophagales" = "purple",
	"Flavobacteriales" = "green",
	"Parvibaculales" = "green3",
	"Pseudomonadales" = "grey30",
	"Rhizobiales" = "red",
	"Rhodobacterales" = "black",
	"Rhodospirillales" = "yellow",
	"Sneathielles" = "white",
	"Sphingobacteriales" = "grey75",
	"Sphingomonadales" = "grey",
	"Thalassobaculales" = "blue2")	
	
#plotting all selected clusters in bar plot ordered by class 
# and displaying orders over time for DNA and RNA
test_groesser_0.5_class <- ggplot(sample_subset, aes(x = new_day, group = order))+
	scale_fill_manual(breaks=levels(sample_subset$order), values = fill_values)+
	geom_bar(data=subset(sample_subset, nucleic_acid == "dna" & treatment == "glyph"),
	  aes(x = new_day-0.5, y = Abundance), fill = "black", width = 0.9, stat = "sum")+
	geom_bar(data=subset(sample_subset, nucleic_acid == "dna" & treatment == "glyph"),
	  aes(x = new_day-0.5, y=Abundance, fill = order), width = .6, stat = "identity")+
	geom_bar(data=subset(sample_subset, nucleic_acid == "cdna" & treatment == "glyph"),
	  aes(x = new_day+0.5, y = Abundance), fill = "black", width = 0.9, stat = "sum")+
	geom_bar(data=subset(sample_subset, nucleic_acid == "cdna" & treatment == "glyph"),
	  aes(x = new_day+0.5, y = Abundance, fill = order), width = 0.6, stat = "identity")+
	geom_vline(data=subset(sample_subset, treatment == "glyph"),aes(xintercept=1.5),
	  linetype = "dashed", size = 1.2)+
	guides(colour =FALSE, size = FALSE, width = FALSE,
		   fill = guide_legend(ncol = 1,
						keyheight = 1.5,
						label.theme = element_text(size = 15,
												face = "italic",
												angle = 0),
						(title = NULL)))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	scale_y_continuous(expand = c(0,0))+
	theme_bw()+
	theme(axis.text = element_text(size = 17))+
	theme(axis.title = element_text(size = 20,face = "bold"))+
	theme(legend.background  =  element_rect(fill = "grey90",linetype = "solid"))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	labs(x="Days", y="Relative abundance [%]")+
  annotate("text", x = -27, y = 90, label = "a)" , color="black", size=6 , 
  angle=0, fontface="bold")+
  annotate("text", x = -22.5, y = 90, label = "b)" , color="black", size=6 , 
  angle=0, fontface="bold")
test_groesser_0.5_class
ggsave(file="fig_02_relative_community_order_level_sorted_for_manuscript.png", 
width=16, height=8)
