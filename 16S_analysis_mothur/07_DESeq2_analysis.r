# running DESeq2 for before/after glyphosate incubation pairs. 
# This will be done for control and treatment samples individually
# http://joey711.github.io/phyloseq-extensions/DESeq2.html

# install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2", version = "3.8")



# Set the working dir with mothur files in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("mothur_deseq_002.RData")

save.image("mothur_deseq_002.RData")

library("phyloseq")
library("DESeq2")
packageVersion("DESeq2")
# generate a phyloseq object from scratch to avoid conflicts with packages

# import mothur output into phyloseq
deseq_ps <- import_mothur(mothur_list_file = NULL, 
						   mothur_group_file = NULL,
						   mothur_tree_file = NULL, 
						   cutoff = NULL, 
						   mothur_shared_file = "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.shared",
						   mothur_constaxonomy_file = "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.0.02.cons.taxonomy", 
						   parseFunction = parse_taxonomy_default)
  
# adjust taxonomy header
tax_table(deseq_ps) <- cbind(tax_table(deseq_ps), 
    rownames(tax_table(deseq_ps)))
colnames(tax_table(deseq_ps)) <- c("kingdom", 
									"phylum", 
									"class", 
									"order", 
									"family", 
									"genus", 
									"otu_id")
rank_names(deseq_ps)

# read meta data, turn into phyloseq object, merge with existing ps object									
metafile <- read.delim("/data/projects/glyphosate/analysis/metadata/all_samples_with_meta_cond2.tsv", 
						row.names = 1, 
						header = TRUE,
						na.strings = "")
metafile2 <- sample_data(metafile)

# add meta data and OTU representative seqs to phyloseq object
deseq_ps2 <- merge_phyloseq(deseq_ps, metafile2)
# test variable is not allowed to contain NA
deseq_ps2 <- subset_samples(deseq_ps2, !(is.na(condition)))
# check levels for test variable
head(sample_data(deseq_ps2)$condition, 50)

### for single tasks

deseq_ps2_glyph_water <- subset_samples(deseq_ps2, treatment == "glyph" & habitat == "water", prune = TRUE)
diagdds = phyloseq_to_deseq2(deseq_ps2_glyph_water, ~ condition)
diagdds$condition <- relevel(diagdds$condition, ref = "untreated")
diagdds = DESeq(diagdds, test = "Wald", fitType = "parametric")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
head(sigtab)




# this is the function we call to split our data into different subsets
get_sample_subsets <- function(ps, nucleic_acid, habitat, threshold, treatment){
	sample_subset <- sample_data(ps)[ which(sample_data(ps)$nucleic_acid == nucleic_acid & 
											sample_data(ps)$habitat == habitat & 
											sample_data(ps)$treatment == treatment),]
	phy_subset <- merge_phyloseq(tax_table(ps), 
								 otu_table(ps),
								 #phy_tree(ps),
								 #refseq(ps),
								 sample_subset)
	phy_subset2 <- filter_taxa(phy_subset, function (x) {sum(x > 0) > threshold}, prune = TRUE)
	return(phy_subset2)
}

# these are the parameters passed to function
ps <- deseq_ps2
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
treatments <- c("glyph", "control")
threshold <- 1


# this is the nested for loop which calls the subsetting function 
# for each combination of subsetting variables
deseq_subsets <- list() 
if(length(deseq_subsets) == 0) {
	for (treatment in treatments){
		for (acid in acids) {
			for (habitat in habitats) {
				print(paste0("nucleic_acid is ", acid, " and habitat is ", 
							 habitat, " and treatment is ", treatment))
				tmp <-	get_sample_subsets(ps = ps, 
									   nucleic_acid = acid, 
									   habitat = habitat, 
									   treatment = treatment, 
									   threshold = threshold)
				sample_data(tmp)$days <- as.factor(sample_data(tmp)$days)					   
				sample_data(tmp)$new_day <- as.factor(sample_data(tmp)$new_day)
				deseq_subsets[[paste(habitat,
										  treatment, 
										  acid, 
										  "min_reads_per_OTU", 
										  threshold, 
										  sep = "_")]] <- tmp
			}
		}
	}
print(deseq_subsets)
} else {
	print("list is not empty, abort to prevend appending...")
}

# if the list was empty, the function now wrote all subsets into the list


# we can now estimate or determine different diversity parameters on our subsets
deseq_tests <- list()						  
counter <- 0
if(length(deseq_tests) == 0 & 
	all.equal(counter, 0)) {						  
deseq_tests <- lapply(deseq_subsets, 
					   function(deseqs) {
									 counter  <<- counter + 1
									 tmp = phyloseq_to_deseq2(deseqs, ~ condition)
									 tmp$condition <- relevel(tmp$condition, ref = "untreated")
									 tmp_dds = DESeq(tmp, test = "Wald", fitType = "parametric")								
})
} else {
	print(paste("list is not empty, or counter not 0 (counter is", counter, 
				"), abort to prevend appending..."))
}
# to adress specific elements in list
results(deseq_tests[[8]])
results(deseq_tests$biofilm_control_cdna_min_reads_per_OTU_1)
results(deseq_tests[["biofilm_control_cdna_min_reads_per_OTU_1"]])

# try cbind with mapply
sigtabs_list <- list()
if(length(sigtabs_list) == 0) {
sigtabs_list <- mapply(function(dds, ps) {res = results(dds, cooksCutoff = FALSE)
									alpha = 0.01
									sigtab = res[which(res$padj < alpha), ]
									sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
									print(head(sigtab))
									return(sigtab)
}, dds = deseq_tests, ps = deseq_subsets, SIMPLIFY = FALSE)
} else {
	print(paste("list is not empty, abort to prevent appending..."))
}


# to adress specific elements in list
names(sigtabs_list)
sigtabs_list[[1]]
sigtabs_list$biofilm_control_cdna_min_reads_per_OTU_1

# extract column from list in R "otu_id" for VENN diagram, get counts with "length"
# https://stackoverflow.com/questions/28305685/get-column-from-list-of-dataframes-r

A <- row.names(sigtabs_list[[1]])
B <- row.names(sigtabs_list[[3]])
A %in% B

intersect(A,B)
setdiff(A,B)
setdiff(B,A)
intersect(intersect(A,B),C)
Reduce(intersect, list(A,B,C))

# plot from otu_id and ps_subset
# we could use the OTU_list function from 04_abs_OTU_abundance_plots
for_melt <- merge_phyloseq(deseq_ps, metafile2)
for_melt2 <- filter_taxa(for_melt, function (x) {sum(x > 0) > 2}, prune = TRUE)
for_melt3 <- transform_sample_counts(for_melt2, function(x){(x / sum(x)) * 100})
deseq_melt <- psmelt(for_melt3)

current_otu_data <- subset(deseq_melt, OTU == "Otu000001")

species_title <- unique(paste(current_otu_data$family, 
							  current_otu_data$genus, 
							  current_otu_data$OTU, 
							  sep = "_"))
							  
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
