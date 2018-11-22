# this is a first test of how to create ordination plots 
# for glyphosate data in R phyloseq

# use the relative abundance phyloseq object, no transformation needed

# Set the working dir with mothur files in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("mothur_glyph.RData")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)



# this is copied from a tutorial from:
# https://joey711.github.io/phyloseq/plot_ordination-examples.html



# careful! subset_samples does not work inside functions 
# https://github.com/joey711/phyloseq/issues/487


# set variables
ps <- mothur_ps4_ra
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
min_OTU_counts <- 5
after_day <- 44
sample_subset_list <- list() 

ordinations <- c("CCA", "NMDS")

# for testing the function
ps = mothur_ps4_ra
nucleic_acid = "cdna"
habitat = "water"
threshold = 5


get_sample_subsets <- function(ps, nucleic_acid, habitat, days, threshold){
	sample_subset <- sample_data(ps)[ which(sample_data(ps2)$nucleic_acid == nucleic_acid & 
											sample_data(ps2)$habitat == habitat & 
											sample_data(ps2)$days > days),]
	phy_subset <- merge_phyloseq(tax_table(ps), 
								 otu_table(ps),
								 phy_tree(ps),
								 refseq(ps),
								 sample_subset)
	phy_subset2 <- filter_taxa(phy_subset, function (x) {sum(x > 0) > threshold}, prune = TRUE)
	return(phy_subset2)
}

if(length(sample_subset_list) == 0) {
	for (acid in acids){
		for (habitat in habitats) {
			print(paste0("nucleic_acid is ", acid, " and habitat is ", habitat))
			sample_subset_list[[length(sample_subset_list)+1]] <- 
				get_sample_subsets(ps = ps, 
								   nucleic_acid = acid, 
								   habitat = habitat, 
								   days = after_day, 
								   threshold = threshold)
		}
	}
print(sample_subset_list)
} else {
	print("list is not empty, abort to prevend appending...")
}






GP.ord <- ordinate(dna_water_5_44, "NMDS", "bray", autotransform = FALSE)
plot_ordination(dna_water_5_44, GP.ord, type="sample", color="days", shape="treatment") + 
	 geom_polygon(aes(fill=disturbance)) + 
	 geom_point(size=5) + 
	 ggtitle("dna_water_nmds") +
	 geom_text(aes(label = new_day), colour = "black", size = 2)


# Remove OTUs that do not show appear more than 5 times in more than half the samples
# wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
# GP1 = prune_taxa(wh0, GP)

# remove data from sample -25-days?? 

ps_dna_water <- subset_samples(mothur_ps4_ra, nucleic_acid == "dna" & 
													habitat == "water" &
													days > 44,
													prune = TRUE)
										  
ps_dna_water1 <- filter_taxa(ps_dna_water, function (x) {sum(x > 0) > 1}, prune = TRUE)
# turn days into factors for plotting
sample_data(ps_dna_water1)$days <- as.factor(sample_data(ps_dna_water1)$days)
GP1 = ps_dna_water1

# (1) Just OTUs

# Let’s start by plotting just the OTUs, and shading the points by Phylum. 
# Note that even in our “trimmed” dataset there are ntaxa(GP1)= 204 OTUs.

#GP.ord <- ordinate(GP1, "NMDS", "bray")
# p1 <- plot_ordination(GP1, GP.ord, type="taxa", color="phylum", title="taxa") +
			# facet_wrap(~phylum, 3)
		
# print(p1)


GP.ord <- ordinate(GP1, "NMDS", "bray", autotransform = FALSE)
plot_ordination(GP1, GP.ord, type="sample", color="days", shape="treatment") + 
	 geom_polygon(aes(fill=disturbance)) + 
	 geom_point(size=5) + 
	 ggtitle("dna_water_nmds") +
	 geom_text(aes(label = new_day), colour = "black", size = 2)
	 
	 
	 
# CCA plot	 
#cca.phyloseq() apparently does not exist anymore

ps_dna_water_44 <- subset_samples(mothur_ps4_ra, nucleic_acid == "dna" & 
													habitat == "water" &
													days > 44,
													prune = TRUE)
										  
ps_dna_water_44_1 <- filter_taxa(ps_dna_water_44, function (x) {sum(x > 0) > 1}, prune = TRUE)
sample_data(ps_dna_water_44_1)$days <- as.factor(sample_data(ps_dna_water_44_1)$days)
test <- ordinate(ps_dna_water_44_1, "CCA", ~ glyphosate, try = 100, autotransform = FALSE)

plot_ordination(ps_dna_water_44_1, test, type="sample", color="days", shape="treatment") + 
	 geom_polygon(aes(fill=disturbance)) + 
	 geom_point(size=5) + 
	 ggtitle("dna_water_cca") +
	 geom_text(aes(label = new_day), colour = "black", size = 3)


## dna_biofilm nmds

ps_dna_biofilm <- subset_samples(mothur_ps4_ra, nucleic_acid == "dna" & 
													habitat == "biofilm" &
													days > 44,
													prune = TRUE)
										  
ps_dna_biofilm1 <- filter_taxa(ps_dna_biofilm, function (x) {sum(x > 0) > 1}, prune = TRUE)
# turn days into factors for plotting
sample_data(ps_dna_biofilm1)$days <- as.factor(sample_data(ps_dna_biofilm1)$days)
GP1 = ps_dna_biofilm1
GP.ord <- ordinate(GP1, "NMDS", "bray", autotransform = FALSE)
plot_ordination(GP1, GP.ord, type="sample", color="days", shape="treatment") + 
	 geom_polygon(aes(fill=disturbance)) + 
	 geom_point(size=5) + 
	 ggtitle("dna_biofilm_nmds") +
	 geom_text(aes(label = new_day), colour = "black", size = 3)
	 
## cdna_biofilm
	 
ps_cdna_biofilm <- subset_samples(mothur_ps4_ra, nucleic_acid == "cdna" & 
													habitat == "biofilm" &
													days > 44,
													prune = TRUE)
										  
ps_cdna_biofilm1 <- filter_taxa(ps_cdna_biofilm, function (x) {sum(x > 0) > 1}, prune = TRUE)
# turn days into factors for plotting
sample_data(ps_cdna_biofilm1)$days <- as.factor(sample_data(ps_cdna_biofilm1)$days)
GP1 = ps_cdna_biofilm1
GP.ord <- ordinate(GP1, "NMDS", "bray", autotransform = FALSE)
plot_ordination(GP1, GP.ord, type="sample", color="days", shape="treatment") + 
	 geom_polygon(aes(fill=disturbance)) + 
	 geom_point(size=5) + 
	 ggtitle("cdna_biofilm_nmds") +
	 geom_text(aes(label = new_day), colour = "black", size = 3)