# how many genera and how many OTUs per genera, abundant OTUs?
str(mothur_ra_melt)

# highest genus abundance per habitat
genus_max <- aggregate(Abundance ~ genus + habitat, data = mothur_ra_melt, max)

# number of OTUs per genus
genus_distribution <- aggregate(Abundance ~ OTU + habitat + genus, 
								data = mothur_ra_melt, max)
# most abundant	OTU within genus per habitat							
genus_distribution_abu <- subset(genus_distribution, Abundance > 0.05)
OTUs_per_genus <- as.data.frame(table(genus_distribution$genus, 
										genus_distribution$habitat))
										
# highest genus abundance per habitat + number of OTU per genus										
OTUs_per_genus_max <- merge(genus_max, OTUs_per_genus, 
							by.x = c("genus", "habitat"), by.y = c("Var1", "Var2"))

OTUs_per_genus_max_freq <- subset(OTUs_per_genus_max, Freq > 2 & Abundance > 0.05)
OTUs_per_genus_max_abu <- droplevels(subset(OTUs_per_genus_max, Abundance > 0.05))

# plot OTUs per genera, ordered by most OTUs
ggplot(OTUs_per_genus_max_freq, aes(x = reorder(genus, -Freq), y = Freq)) +
	geom_bar(stat = "identity") +
	#coord_cartesian(ylim = c(0, 100)) +
	theme(axis.text = element_text(size = 17, angle = 90))

# plot abundance per genera, ordered by abundance
ggplot(OTUs_per_genus_max_abu, aes(x = reorder(genus, -Abundance), y = Abundance)) +
	geom_bar(stat = "identity") +
	#coord_cartesian(ylim = c(0, 1)) +
	theme(axis.text = element_text(size = 13, angle = 90))

# plot abundance per OTU, colored by genus
ggplot(genus_distribution_abu, aes(x = OTU, y = Abundance, 
	color = genus, fill = order)) +
	geom_bar(stat = "identity") +
	geom_text(aes(y = 40, label = genus), size = 3, angle = 85) +
	facet_wrap(~ habitat, nrow = 2) +
	theme(axis.text = element_text(size = 10, angle = 90)) +
	theme(legend.position = "none")
	
# distribution of OTUs per genus


# OTUs of certain genus
genus_distribution[grep("Pseudomonas", genus_distribution$genus),]

OTUs_per_genus_max[grep("Escherichia", OTUs_per_genus_max$genus),]

# use dplyr?
# max abundance of an OTU of a certain genus
max(subset(mothur_ra_melt, grepl('Pseudomonas', genus))$Abundance)

# generate subset of specific genus
Pseudomonas <- droplevels(subset(mothur_ra_melt, 
								   grepl('Pseudomonas', genus)))							   
str(Pseudomonas)

								   
length(table(Pseudomonas_water_plot$OTU))
								
# find abundant OTUs within genus
head(b[order(-b$Abundance),], 30)