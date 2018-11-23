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
ordinations <- c("CCA", "NMDS")

ps <- mothur_ps4_ra
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
threshold = 5
after_day <- c(43, 44)








get_sample_subsets <- function(ps, nucleic_acid, habitat, days, threshold){
	sample_subset <- sample_data(ps)[ which(sample_data(ps)$nucleic_acid == nucleic_acid & 
											sample_data(ps)$habitat == habitat & 
											sample_data(ps)$days > days),]
	phy_subset <- merge_phyloseq(tax_table(ps), 
								 otu_table(ps),
								 phy_tree(ps),
								 refseq(ps),
								 sample_subset)
	phy_subset2 <- filter_taxa(phy_subset, function (x) {sum(x > 0) > threshold}, prune = TRUE)
	return(phy_subset2)
}

sample_subset_list <- list() 
if(length(sample_subset_list) == 0) {
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
				sample_subset_list[[paste(habitat, 
										  "after day", 
										  each_day, 
										  acid, 
										  "min reads per OTU", 
										  threshold, 
										  sep = " ")]] <- tmp
			}
		}
	}
print(sample_subset_list)
} else {
	print("list is not empty, abort to prevend appending...")
}

# TO DO:
# integrate checks for: counter = 0, list is empty
# I should turn this very similar functions into a bigger function or loop
ordination_nmds <- list()
ordination_nmds <- lapply(sample_subset_list, ordinate, method = "NMDS", dist = "bray", try = 100, autotransform = FALSE)

ordination_cca <- list()
ordination_cca <- lapply(sample_subset_list, ordinate, method = "CCA", try = 100, formula = ~ glyphosate, autotransform = FALSE)

ordination_rda <- list()
ordination_rda <- lapply(sample_subset_list, ordinate, method = "RDA", try = 100, formula = ~ glyphosate, autotransform = FALSE)

# this security check works
nmds_ordination_plots <- list()
counter <- 0
if(length(nmds_ordination_plots) == 0 & 
	identical(counter, integer(0))) {
nmds_ordination_plots <- mapply(function(x,y) {
						 counter <<- counter + 1 
						 plot_ordination(x, y, 
										 type = "sample", 
										 color = "new_day", 
										 shape = "treatment") + 
							geom_polygon(aes(fill = disturbance)) + 
							geom_point(size = 5) + 
							guides(color = FALSE) +
							ggtitle(names(sample_subset_list)[counter]) +
							geom_text(aes(label = new_day), 
									  colour = "black", 
									  size = 2)
}, x = sample_subset_list, 
   y = ordination_nmds, 
   SIMPLIFY = FALSE)
} else {
	print(paste("list is not empty, or counter not 0 (counter is", counter, 
				"), abort to prevend appending..."))
}

# not implemented in the others so far
cca_ordination_plots <- list()
counter <- 0
cca_ordination_plots <- mapply(function(x,y) {
						print(paste("counter")
						counter <<- counter + 1 
						plot_ordination(x, y, type = "sample", color = "new_day", shape="treatment") + 
							geom_polygon(aes(fill = disturbance)) + 
							geom_point(size = 5) + 
							guides(color = FALSE) +
							ggtitle(names(sample_subset_list)[counter]) +
							geom_text(aes(label = new_day), colour = "black", size = 2)
}, x = sample_subset_list, y = ordination_cca, SIMPLIFY = FALSE)


rda_ordination_plots <- list()
counter <- 0
rda_ordination_plots <- mapply(function(x,y) {
						counter <<- counter + 1 
						plot_ordination(x, y, type = "sample", color = "new_day", shape="treatment") + 
							geom_polygon(aes(fill = disturbance)) + 
							geom_point(size = 5) + 
							guides(color = FALSE) +
							ggtitle(names(sample_subset_list)[counter]) +
							geom_text(aes(label = new_day), colour = "black", size = 2)
}, x = sample_subset_list, y = ordination_rda, SIMPLIFY = FALSE)

# samples 5:8 are those without the very early and different sample from day 44

# NMDS plots
require(gridExtra)
do.call("grid.arrange", c(nmds_ordination_plots[c(1:4)], nrow = 2, top = "NMDS"))	
do.call("grid.arrange", c(nmds_ordination_plots[c(5:8)], nrow = 2, top = "NMDS"))



# RDA plots
rda_with_constraint <- rda_ordination_plots
rda_no_constraint <- rda_ordination_plots

do.call("grid.arrange", c(rda_with_constraint[c(5:8)], nrow = 2, top = "rda_with_constraint"))
do.call("grid.arrange", c(rda_no_constraint[c(5:8)], nrow = 2, top = "rda_no_constraint"))	


# CCA plots
with_constraint <- cca_ordination_plots
no_constraint <- cca_ordination_plots

do.call("grid.arrange", c(with_constraint[c(5:8)], nrow = 2, top = "cca_with_constraint"))
do.call("grid.arrange", c(no_constraint[c(5:8)], nrow = 2, top = "cca_no_constraint"))	


################# single plots for testing
GP.ord <- ordinate(dna_water_5_44, "NMDS", "bray", autotransform = FALSE)
plot_ordination(dna_water_5_44, GP.ord, type="sample", color="days", shape="treatment") + 
	 geom_polygon(aes(fill=disturbance)) + 
	 geom_point(size=5) + 
	 ggtitle("dna_water_nmds") +
	 geom_text(aes(label = new_day), colour = "black", size = 2)


ps_dna_water <- subset_samples(mothur_ps4_ra, nucleic_acid == "dna" & 
													habitat == "water" &
													days > 44,
													prune = TRUE)
										  
ps_dna_water1 <- filter_taxa(ps_dna_water, function (x) {sum(x > 0) > 1}, prune = TRUE)

GP2 = GP1 = ps_dna_water1
# turn days into factors for plotting
sample_data(GP2)$days <- as.factor(sample_data(GP1)$days)
sample_data(GP2)$new_day <- as.factor(sample_data(GP1)$new_day)

GP.ord <- ordinate(GP1, "NMDS", "bray", autotransform = FALSE)
plot_ordination(GP2, GP.ord, type="sample", color="new_day", shape="treatment") + 
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