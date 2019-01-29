# this is a first test of how to create ordination plots 
# for glyphosate data in R phyloseq

# use the relative abundance phyloseq object, no transformation needed

# Set the working dir with mothur files in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("mothur_glyph_002.RData")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)



# this is copied from a tutorial from:
# https://joey711.github.io/phyloseq/plot_ordination-examples.html



# careful! subset_samples does not work inside functions 
# https://github.com/joey711/phyloseq/issues/487


# set variables
#ordinations <- c("NMDS", "CCA")
mothur_ps_bla <- filter_taxa(mothur_ps2, function (x) {sum(x > 2) >= 1}, prune = TRUE)
mothur_ps2_ra <- transform_sample_counts(mothur_ps_bla, function(x){(x / sum(x)) * 100})

ps <- mothur_ps2_ra 
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
threshold <- 0
after_day <- 43



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


# formula <- list("~ glyphosate", NULL)

# ordination_list <- list()
	# for (current_formula in formula) {
							# tmp <- lapply(sample_subset_list, 
										  # ordinate, 
										  # method = "CCA",
										  # formula = as.formula(current_formula))
							# ordination_list[[paste(current_formula)]] <- tmp
# }



# ordination_list <- list()
# for (current_formula in formula) {
  # tmp <- lapply(sample_subset_list, 
                # ordinate, 
                # method = "CCA",
                # formula = as.formula(current_formula)) 
                # ordination_list[[length(ordination_list) + 1]] <- tmp
# }



# ordination_list <- list()
# for (current_formula in formula) {
  # tmp <- lapply(sample_subset_list, 
                # ordinate, 
                # method = "CCA",
                # formula = if(is.null(current_formula)) NULL else as.formula(current_formula))
                # ordination_list[[length(ordination_list) + 1]] <- tmp
# }


							# tmp <- lapply(sample_subset_list, 
										  # ordinate, 
										  # method = "CCA",
										  # formula = as.formula(formula[[2]]))
							# ordination_list[[paste(current_formula)]] <- tmp

# formula <- list("~ Al", NULL)
							
# for (current_formula in formula) {
# cca(paste(varespec, as.formula(current_formula), sep = " "), data = varechem)
# }

# ps2 <- filter_taxa(ps, function (x) {sum(x > 0) > 50}, prune = TRUE)
# test_otu <- as.data.frame(otu_table(ps2))
# test_sample <- phyloseq_to_df(sample_data(ps2))
# str(ordination_list, max = 2)


# ordination_list <- list()
# counter <- 0
# for (method in ordinations) {
	# for (current_formula in formula) {
							
							# tmp <- lapply(sample_subset_list, ordinate, 
											  # method = method, 
											  # dist = "bray", 
											  # try = 2,
											  # maxtry = 5,
											  # formula = current_formula,
											  # autotransform = FALSE)
							# ordination_list[[paste(method, current_formula)]] <- tmp
# }
# }
# str(ordination_list, max = 2)
 
				


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
# I should turn this very similar functions into a bigger function or loop
ordination_nmds <- list()
ordination_nmds <- lapply(sample_subset_list, ordinate, 
											  method = "NMDS", 
											  dist = "bray", 
											  try = 100, 
											  autotransform = TRUE)
											  
# ordination_nmds <- lapply(sample_subset_list, ordinate, 
											  # method = "NMDS", 
											  # dist = "unifrac",
											  # weighted = TRUE)
											  

# ordination_cca <- list()
# ordination_cca <- lapply(sample_subset_list, ordinate, 
											 # method = "CCA", 
											 # formula = ~ glyphosate + Condition(days),
											 # try = 100, 
											 # autotransform = FALSE)

# ordination_rda <- list()
# ordination_rda <- lapply(sample_subset_list, ordinate, 
											 # method = "RDA", 
											 # try = 100, 
											 # formula = ~ glyphosate, 
											 # autotransform = FALSE)

# NMDS function
nmds_ordination_plots <- list()
counter <- 0
if(length(nmds_ordination_plots) == 0 & 
	all.equal(counter, 0)) {
nmds_ordination_plots <- mapply(function(x,y) {
						 counter <<- counter + 1 
						 plot_ordination(x, y, 
										 type = "sample", 
										 color = "new_day", 
										 shape = "treatment") + 
							geom_polygon(aes(fill = disturbance), alpha = 0.5, size = 0.01) + 
							geom_point(size = 6) + 
							guides(color = FALSE) +
							ggtitle(names(sample_subset_list)[counter]) +
							geom_text(aes(label = new_day), 
									  colour = "black", 
									  size = 3.5) +
							theme(axis.text = element_text(size = 18),
								  axis.title = element_text(size = 20, face = "bold"),
								  legend.title = element_text(size = 15, face = "bold"), 
								  legend.text = element_text(size = 13))
}, x = sample_subset_list, 
   y = ordination_nmds, 
   SIMPLIFY = FALSE)
} else {
	print(paste("list is not empty, or counter not 0 (counter is", counter, 
				"), abort to prevend appending..."))
}
do.call("grid.arrange", c(nmds_ordination_plots[c(1, 3)], nrow = 1, top = "NMDS"))
do.call("grid.arrange", c(nmds_ordination_plots[c(2, 4)], nrow = 1, top = "NMDS"))

# CCA, RDA function
# cca_ordination_plots <- list()
# counter <- 0
# if(length(cca_ordination_plots) == 0 & 
	# all.equal(counter, 0)) {
# cca_ordination_plots <- mapply(function(x,y) {
						# counter <<- counter + 1 
						# plot_ordination(x, y, type = "sample", color = "new_day", shape="treatment") + 
							# geom_polygon(aes(fill = disturbance)) + 
							# geom_point(size = 5) + 
							# guides(color = FALSE) +
							# ggtitle(names(sample_subset_list)[counter]) +
							# geom_text(aes(label = new_day), colour = "black", size = 2)
# }, x = sample_subset_list, y = ordination_cca, SIMPLIFY = FALSE)} else {
	# print(paste("list is not empty, or counter not 0 (counter is", counter, 
				# "), abort to prevend appending..."))
# }

# samples 5:8 are those without the very early and different sample from day 44
plot_folder <- "/data/projects/glyphosate/plots/R/ordination/"
# NMDS plots
require(gridExtra)
# do.call("grid.arrange", c(nmds_ordination_plots[c(1:4)], nrow = 2, top = "NMDS"))	
# ggsave(file = paste(plot_folder, threshold,"_nmds_after_43.png", 
								  # sep = ""),
								  # height = 13,
								  # width = 20)
#do.call("grid.arrange", c(nmds_ordination_plots[c(1:2)], nrow = 1, top = "NMDS"))								  
g1 <- do.call("arrangeGrob", c(nmds_ordination_plots[c(1,3)], nrow = 1, top = "NMDS"))	
ggsave(g1, file = paste(plot_folder, threshold,"water_nmds_after_43_AT_larger_font.png", 
								  sep = ""),
								  height = 10,
								  width = 20)
#do.call("grid.arrange", c(nmds_ordination_plots[c(5:8)], nrow = 2, top = "NMDS"))								  
g2 <- do.call("arrangeGrob", c(nmds_ordination_plots[c(2,4)], nrow = 2, top = "NMDS"))
ggsave(g2, file = paste(plot_folder, threshold,"biofilm_nmds_after_43_AT_larger_font.png", 
								  sep = ""),
								  height = 10,
								  width = 20)
						  
# do.call("arrangeGrob", c(nmds_ordination_plots[c(5:8)], nrow = 2, top = "NMDS"))
# ggsave(file = paste(plot_folder, "nmds_after_44.png", 
								  # sep = ""),
								  # height = 13,
								  # width = 20)
# CCA plots
# g <- do.call("arrangeGrob", c(cca_ordination_plots[c(1:4)], nrow = 2, top = "cca_two_conditions"))
# ggsave(g, file = paste(plot_folder, "cca_after_43_two_conditions.png", 
								  # sep = ""),
								  # height = 13,
								  # width = 20)
# g <- do.call("arrangeGrob", c(cca_ordination_plots[c(5:8)], nrow = 2, top = "cca_two_conditions"))
# ggsave(g, file = paste(plot_folder, "cca_after_44_two_conditions.png", 
								  # sep = ""),
								  # height = 13,
								  # width = 20)
								  
# do.call("grid.arrange", c(no_constraint[c(5:8)], nrow = 2, top = "cca_no_constraint"))	

# RDA plots
# do.call("grid.arrange", c(rda_ordination_plots[c(1:4)], nrow = 2, top = "rda_with_constraint"))
# ggsave(file = paste(plot_folder, "rda_after_43_constrained.png", 
								  # sep = ""),
								  # height = 13,
								  # width = 20)
# do.call("grid.arrange", c(rda_ordination_plots[c(5:8)], nrow = 2, top = "rda_with_constraint"))
# ggsave(file = paste(plot_folder, "rda_after_44_constrained.png", 
								  # sep = ""),
								  # height = 13,
								  # width = 20)


#scp -r -i /drives/d/ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/plots/R/ordination/* /mnt/d/denbi/chandler/ordination/


dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
        ordi = ordinate(physeq, method=i, distance=dist)
        plot_ordination(physeq, ordi, "samples", color="SampleType")
}, GP1, dist)

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