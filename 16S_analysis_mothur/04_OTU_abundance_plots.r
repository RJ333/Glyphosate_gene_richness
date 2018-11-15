Pseudomonas_water_plot <- subset(mothur_ra_melt, 
								   grepl('Pseudomonas', genus) 
								  # & OTU == "Otu000011"
								   & days > 40)
								   
b <- aggregate(Abundance ~ OTU + nucleic_acid + habitat, 
								data = Pseudomonas_water_plot, 
								max)
# number of OTUs per genus								   
length(table(Pseudomonas_water_plot$OTU))
								
# find abundant OTUs within genus
head(b[order(-b$Abundance),], 30)

14425 Otu000006         cdna   water 41.1363379
21639 Otu000009          dna   water 26.2518968
14426 Otu000007         cdna   water 24.5058737
21637 Otu000006          dna   water 23.5841785
14428 Otu000019         cdna   water 23.3428603
21638 Otu000007          dna   water 17.5272880
14427 Otu000009         cdna   water 16.9166368
21640 Otu000019          dna   water 10.3588078
14430 Otu000028         cdna   water  4.6808622
14429 Otu000024         cdna   water  4.1757185
14436 Otu000050         cdna   water  3.0617846
14432 Otu000034         cdna   water  2.5944663
21641 Otu000024          dna   water  2.5719397
21642 Otu000028          dna   water  2.5626722
21644 Otu000034          dna   water  2.3585519
21645 Otu000035          dna   water  1.3154524
21648 Otu000050          dna   water  0.9649769
14434 Otu000036         cdna   water  0.9285809


# with parallels
Pseudomonas_water_plot <- subset(mothur_ra_melt,
								   habitat == "water" 
								  # & grepl('Pseudomonas', genus) 
								   & OTU == "Otu000006"
								   & days > 40)
							   
species_title_Pseudomonas<-expression(paste(,italic("Pseudomonas"), " sp."))
Pseudomonas_water_plot$treatment2 <- factor(Pseudomonas_water_plot$treatment, 
											  labels = c("Control", "Treatment"))
		
ggplot(data = Pseudomonas_water_plot, aes(x = days - 69, y = Abundance, 
			group = nucleic_acid, lty = nucleic_acid))+ 
	#coord_cartesian(ylim = c(0, 1))+
	#annotate("segment", x = -8.0, xend = 18, y = 4300000, yend = 4300000, 
	#		 colour = "grey50", size = 2, alpha = 1,lineend = "butt")+
	#annotate("segment", x = 18.0, xend = 18, y = 4440000, yend = 4160000, 
	#		 colour = "grey50", size = 2, alpha = 1,lineend = "butt")+
	#annotate("segment", x = -8.0, xend = -8, y = 4440000, yend = 4160000, 
	#		 colour = "grey50", size = 2, alpha = 1,lineend = "butt")+
	geom_vline(aes(xintercept = 0), linetype = "dashed", size=1.2)+
	geom_point(data = subset(Pseudomonas_water_plot, treatment == "control"), 
		       aes(colour = treatment), alpha = 1)+
	stat_summary(data = subset(Pseudomonas_water_plot, treatment == "control"), 
	             aes(colour = treatment), fun.y = "mean",  geom = "line",  size = 2, alpha = 1)+
	stat_summary(data = subset(Pseudomonas_water_plot, treatment == "glyph"), 
	             aes(colour = treatment), fun.y = "mean",  geom = "line",  size = 2)+
	geom_point(data = subset(Pseudomonas_water_plot, treatment == "glyph"), 
	             aes(colour = treatment))+
	scale_linetype_manual(values = c("dna" = 1, "cdna" = 6), 
						name = "Nucleic acid  ", 
						breaks = c("cdna", "dna"), 
						labels = c("16S rRNA", "16S rRNA gene"))+
	scale_colour_manual(values = c("glyph" = "black", "control" = "grey50"), 
						name = "Microcosm  ", 
						breaks = c("glyph", "control"), 
						labels = c("Treatment", "Control"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	theme_bw()+
	#scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", 
	#					gsub("e", " %*% 10^", scientific_format()(x)))))})+
	theme(axis.text = element_text(size = 18))+
	theme(panel.grid.major = element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor = element_line(colour = NA, size = 0.5))+
	theme(legend.position = "none")+
	theme(axis.title = element_blank())

	
	
### without parallels, does not work yet, grouping is wrong	
Gallaecimonas_water_plot <- subset(mothur_ra_melt_mean,
								   habitat == "water" 
								   & grepl('Gallaecimonas', genus) 
								   & days > 40)

# number of OTUs per genus								   
length(table(Gallaecimonas_water_plot$OTU))
								   
species_title_Gallaecimonas<-expression(paste(,italic("Gallaecimonas"), " sp."))
Gallaecimonas_water_plot$treatment2 <- factor(Gallaecimonas_water_plot$treatment, 
											  labels = c("Control", "Treatment"))

a <- ggplot(data = Gallaecimonas_water_plot, aes(x = days - 69, y = Abundance, colour = OTU))+ 
	#coord_cartesian(ylim = c(0, 1))+
	#annotate("segment", x = -8.0, xend = 18, y = 4300000, yend = 4300000, 
	#		 colour = "grey50", size = 2, alpha = 1,lineend = "butt")+
	#annotate("segment", x = 18.0, xend = 18, y = 4440000, yend = 4160000, 
	#		 colour = "grey50", size = 2, alpha = 1,lineend = "butt")+
	#annotate("segment", x = -8.0, xend = -8, y = 4440000, yend = 4160000, 
	#		 colour = "grey50", size = 2, alpha = 1,lineend = "butt")+
	geom_vline(aes(xintercept = 0), linetype = "dashed", size=1.2)+
	geom_line(data = subset(Gallaecimonas_water_plot, treatment == "control"), 
		       aes(colour = treatment), alpha = 1)+
	geom_line(data = subset(Gallaecimonas_water_plot, treatment == "glyph"), 
	             aes(colour = treatment))+
	scale_linetype_manual(values = c("dna" = 1, "cdna" = 6), 
						name = "Nucleic acid  ", 
						breaks = c("cdna", "dna"), 
						labels = c("16S rRNA", "16S rRNA gene"))+
	scale_colour_manual(values = c("glyph" = "black", "control" = "grey50"), 
						name = "Microcosm  ", 
						breaks = c("glyph", "control"), 
						labels = c("Treatment", "Control"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	theme_bw()+
	#scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", 
	#					gsub("e", " %*% 10^", scientific_format()(x)))))})+
	theme(axis.text = element_text(size = 18))+
	theme(panel.grid.major = element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor = element_line(colour = NA, size = 0.5))+
	theme(legend.position = "none")+
	theme(axis.title = element_blank())