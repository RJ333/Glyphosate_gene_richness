# this is a collection of useful commands to explore your phyloseq dataset
# in terms of diversity and richness. When I've figured out better, what values
# are ultimately needed for comparison, I can adjust the script


# let's generate an own ps-object for the diversity tests and
# remove OTUs with less than 2 reads for diversity and richness analysis
# mothur_div <- filter_taxa(mothur_ps2, function (x) {sum(x > 0) >= 1}, prune = TRUE)
# Set the working dir with mothur files in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("mothur_glyph_002.RData")

.cran_packages <- c("ggplot2", 
					"gridExtra")
.bioc_packages <- c("dada2", 
					"phyloseq", 
					"DECIPHER", 
					"phangorn")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# include singletons
mothur_ps2
mothur_div <- mothur_ps2
# this is the function we call to split our data into different subsets
get_sample_subsets <- function(ps, nucleic_acid, habitat, days, threshold){
	sample_subset <- sample_data(ps)[ which(sample_data(ps)$nucleic_acid == nucleic_acid & 
											sample_data(ps)$habitat == habitat & 
											sample_data(ps)$days > days),]
	phy_subset <- merge_phyloseq(tax_table(ps), 
								 otu_table(ps),
								 #phy_tree(ps),
								 refseq(ps),
								 sample_subset)
	phy_subset2 <- filter_taxa(phy_subset, function (x) {sum(x > threshold) >= 1}, prune = TRUE)
	return(phy_subset2)
}

# these are the parameters passed to function
ps <- mothur_div
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
threshold <- 0
after_day <- 43

# this is the nested for loop which calls the subsetting function 
# for each combination of subsetting variables
div_subset_list <- list() 
if(length(div_subset_list) == 0) {
	for (each_day in after_day){
		for (acid in acids) {
			for (habitat in habitats) {
				print(paste0("nucleic_acid is ", acid, " and habitat is ", 
							 habitat, " and first day is ", each_day))
				tmp <-	get_sample_subsets(ps = ps, 
									   nucleic_acid = acid, 
									   habitat = habitat, 
									   days = each_day, 
									   threshold = threshold)
				sample_data(tmp)$days <- as.factor(sample_data(tmp)$days)					   
				sample_data(tmp)$new_day <- as.factor(sample_data(tmp)$new_day)
				div_subset_list[[paste(habitat, 
										  "after day", 
										  each_day, 
										  acid, 
										  "min reads per OTU", 
										  threshold, 
										  sep = " ")]] <- tmp
			}
		}
	}
print(div_subset_list)
} else {
	print("list is not empty, abort to prevend appending...")
}

# we can now estimate or determine different diversity parameters on our subsets
all_measures <- list()						  
counter <- 0
if(length(all_measures) == 0 & 
	all.equal(counter, 0)) {						  
all_measures <- lapply(div_subset_list, 
					   function(ps) {
									 counter  <<- counter + 1
									 plot_richness(ps,
										x = "new_day",
										title = names(div_subset_list)[counter],
										color = "treatment",
										#scales = "fixed",
										measures = c("Observed", 
													 "Chao1", 
													 "ACE", 
													 "Shannon", 
													 "Simpson", 
													 "InvSimpson", 
													 "Fisher"))
})
} else {
	print(paste("list is not empty, or counter not 0 (counter is", counter, 
				"), abort to prevend appending..."))
}
# if the list was empty and the counter 0 (import to get the matching plot title)
# we now have a plot with all measure per subset in the list

# we can plot those individually using 
all_measures[1]
# or altogether (not a good idea here)
require(gridExtra)
do.call("grid.arrange", c(all_measures, 
						  nrow = 4, 
						  top = "Chao1 est richness and Shannon Index"))

# this function allows us to plot one diversity measure for all samples	  
shannon <- list()						  
counter <- 0
if(length(shannon) == 0 & 
	all.equal(counter, 0)) {						  
shannon <- lapply(div_subset_list, 
				  function(ps) {
								counter  <<- counter + 1
								plot_richness(ps,
									 x = "new_day",
									 title = names(div_subset_list)[counter],
									 color = "treatment",
									 scales = "fixed",
									 measures = "Shannon")
})
} else {
	print(paste("list is not empty, or counter not 0 (counter is", counter, 
				"), abort to prevend appending..."))
}

do.call("grid.arrange", c(shannon, 
						  nrow = 2, 
						  top = "Chao1 est richness and Shannon Index"))


# we can perform a t test to check for significantly different measures, 
# based on a factor in the sample data
# ttest https://github.com/joey711/phyloseq/issues/704 
# don't try with list objects
# set day 62 and 69 to untreated and 72, 76, 79, 83 ,86 to treated
# sampleData(GP)$human <- getVariable(GP, "SampleType") %in% c("Feces", "Mock", "Skin", "Tongue")

ttest_list <- list()
richness_subset_list <- list()
richness_subset_list[["water_cdna_control"]] <- watercdnacontrol2
richness_subset_list[["water_cdna_glyph"]] <- watercdnaglyph2
richness_subset_list[["water_dna_control"]] <- waterdnacontrol2
richness_subset_list[["water_dna_glyph"]] <- waterdnaglyph2


# set factor levels for tests, first for plotting, second for ttest
sample_data(mothur_div)$condition_diversity[sample_data(mothur_div)$condition_diversity == "start"] <- "untreated"
sample_data(mothur_div)$condition[sample_data(mothur_div)$new_day == -25] <- "untreated"

# define order of factor levels
condition_order <- c("untreated", "treated", "22 to 36", "43 to 71")

waterdnacontrol <- subset_samples(mothur_div, habitat == "water" & 
									   nucleic_acid == "dna" &
									   treatment == "control" &
									   days > 43,
									   prune = TRUE)

# calculate measures
erich_waterdnacontrol <- estimate_richness(waterdnacontrol, measures = c("Observed", 
																	 "Chao1", 
																	 "ACE", 
																	 "Shannon", 
																	 "Simpson", 
																	 "InvSimpson", 
																	 "Fisher"))
# perform t test with measures												   
ttest_waterdnacontrol <- t(sapply(erich_waterdnacontrol, function(x) unlist(t.test(x~sample_data(waterdnacontrol)$condition)[c("estimate",
																														 "p.value",
																														 "statistic",
																														 "conf.int")])))																					 

ttest_list[["dnawatercontrol"]] <- ttest_waterdnacontrol	
#richness_subset_list[["biofilm_dna_control"]] <- biofilmdnacontrol2




# plot measures only divided by condition. NA marks samples neither treated or untreated
richness_plot_dna_water_control <- plot_richness(waterdnacontrol, x = "condition_diversity", 
						 color = "as.factor(new_day)",
						 measures = c(
									 #"Observed", 
									 #"Chao1", 
									  #"ACE", 
									  "Shannon" 
									  #"Simpson", 
									  #"InvSimpson", 
									 # "Fisher"
									  )) + geom_point(alpha = 0.8, size = 4) +
									  ggtitle("richness_plot_dna_water_control") + #+ geom_violin(alpha = 0.5)
									      theme(axis.text = element_text(size = 18),
								  axis.title = element_text(size = 20, face = "bold"),
								  legend.title = element_text(size = 15, face = "bold"), 
								  legend.text = element_text(size = 13))

richness_plot_dna_water_control$data$condition_diversity <- as.character(richness_plot_dna_water_control$data$condition_diversity)
richness_plot_dna_water_control$data$condition_diversity <- factor(richness_plot_dna_water_control$data$condition_diversity, levels = condition_order)
print(richness_plot_dna_water_control)	


# put into list and plot together
plots_div <- list()
plots_div <- list(richness_plot_cdna_water_glyph, richness_plot_dna_water_glyph,  richness_plot_cdna_water_control, richness_plot_dna_water_control)							  
plot_folder <- "/data/projects/glyphosate/plots/R/diversity/"
div_plot_object <- do.call("arrangeGrob", c(plots_div, nrow = 2, top = "Shannon Diversity"))
ggsave(div_plot_object, file = paste(plot_folder, "Shannon_water_after43_larger_font.png", 
								  sep = ""),
								  height = 13,
								  width = 20)
								  
#scp -r -i /drives/d/ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/plots/R/diversity/* /mnt/d/denbi/chandler/diversity/
########################### general commands to get information on genera and OTUs and their distribution


# for biofilm
genus_max <- aggregate(Abundance ~ OTU + genus + family + order + habitat + nucleic_acid, data = mothur_ra_melt, max)
genus_max_biofilm <- subset(genus_max, habitat == "biofilm" & order == "Rhizobiales")
max(genus_max_biofilm$Abundance)
genus_max_biofilm <- genus_max_biofilm[order(-genus_max_biofilm$Abundance),]
										  								  
# how many genera and how many OTUs per genera, abundant OTUs?
str(mothur_ra_melt)

# highest genus abundance per habitat
genus_max <- aggregate(Abundance ~ OTU + genus + habitat, data = mothur_ra_melt, max)
#genus_max <- aggregate(Abundance ~ genus + family + habitat, data = deseq_melt, max)
genus_max$Abundance <- round(genus_max$Abundance, 2)
genus_max[order(-genus_max$Abundance),]

# number of OTUs per genus
genus_distribution <- aggregate(Abundance ~ OTU + habitat + genus, 
								data = mothur_ra_melt, 
								max)
# most abundant	OTU within genus per habitat							
genus_distribution_abu <- subset(genus_distribution, Abundance > 0.05)
OTUs_per_genus <- as.data.frame(table(genus_distribution$genus, 
										genus_distribution$habitat))
										
# highest genus abundance per habitat + number of OTU per genus										
OTUs_per_genus_max <- merge(genus_max, 
							OTUs_per_genus, 
							by.x = c("genus", 
									 "habitat"), 
							by.y = c("Var1", 
									 "Var2"))
# subsets genera with more than 2 OTUs
OTUs_per_genus_max_freq <- subset(OTUs_per_genus_max, Freq > 2 & Abundance > 0.05)
# or only above a threshold
OTUs_per_genus_max_abu <- droplevels(subset(OTUs_per_genus_max, Abundance > 0.05))

# plot OTUs per genera, ordered by most OTUs
ggplot(OTUs_per_genus_max_freq, 
	   aes(x = reorder(genus, -Freq), 
		   y = Freq)) +
	   geom_bar(stat = "identity") +
	   #coord_cartesian(ylim = c(0, 100)) +
	   theme(axis.text = element_text(size = 17, 
								      angle = 90))

# plot abundance per genera, ordered by abundance
ggplot(OTUs_per_genus_max_abu, 
	   aes(x = reorder(genus, -Abundance), 
		   y = Abundance)) +
	   geom_bar(stat = "identity") +
	   #coord_cartesian(ylim = c(0, 1)) +
	   theme(axis.text = element_text(size = 13, 
									  angle = 90))

# plot abundance per OTU, colored and named by genus
ggplot(genus_distribution_abu, 
	   aes(x = OTU, 
	       y = Abundance, 
		   color = genus, 
		   fill = order)) +
	   geom_bar(stat = "identity") +
	   geom_text(aes(y = 40, 
					 label = genus), 
				 size = 3, 
				 angle = 85) +
	   facet_wrap(~ habitat, 
			      nrow = 2) +
	   theme(axis.text = element_text(size = 10, 
									  angle = 90)) +
	   theme(legend.position = "none")
	
# these commands show you the OTUs per specific genus
nrow(genus_distribution[grep("Pseudomonas", genus_distribution$genus),])
OTUs_per_genus_max[grep("Escherichia", OTUs_per_genus_max$genus),]



# max abundance of an OTU of a certain genus
max(subset(mothur_ra_melt, grepl('Pseudomonas', genus))$Abundance)