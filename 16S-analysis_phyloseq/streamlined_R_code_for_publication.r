#!/usr/bin/env Rscript
# load libraries
library(scales)
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(gridExtra)
library(Biostrings)

# Set the working dir
# it should contain shared file, constaxonomy file,
# meta file, cell count file and OTU_rep fasta file in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")

# set path for plots, create directory if not existant
plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

# test for sourcing
source_test <- "/home/centos/scripts/Glyphosate_gene_richness/16S-analysis_phyloseq/source_test_read_mothur.r"
source(source_test)

# mothur files and additional files that need to be imported from the working dir
# shared file is the OTU count table
# constaxonomy contains the taxonomy of the OTUs
our_shared_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.shared"
our_cons_taxonomy_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.0.02.cons.taxonomy"

# sample data is in the metafile
metafile <- read.delim("/data/projects/glyphosate/analysis/metadata/metafile.tsv", row.names = 1, header = TRUE, na.strings = "")
metafile <- sample_data(metafile)

# read table for cell count plot
cell_counts_glyph <- read.csv("cell_counts_glyph.csv", sep = ";")

# read OTU representative sequences
OTU_seqs <- readDNAStringSet(file = "OTU_reps_fasta_002.fasta",
  format = "fasta", nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

# import mothur output into phyloseq
mothur_phyloseq_object <- import_mothur(mothur_list_file = NULL,
  mothur_group_file = NULL, mothur_tree_file = NULL, cutoff = NULL,
  mothur_shared_file = our_shared_file,
  mothur_constaxonomy_file = our_cons_taxonomy_file,
  parseFunction = parse_taxonomy_default)

# now all files are imported, we can adjust them to our needs
# add further taxonomy columns "OTU" and "wholetax" and adjust column names
wholetax <- do.call(paste, c(as.data.frame(tax_table(mothur_phyloseq_object))
  [c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6")], sep = "_"))
tax_table(mothur_phyloseq_object) <- cbind(tax_table(mothur_phyloseq_object),
  rownames(tax_table(mothur_phyloseq_object)), wholetax)
colnames(tax_table(mothur_phyloseq_object)) <- c(
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "otu_id",
  "wholetax")
                                    
# add meta data and OTU representative seqs to phyloseq_object
mothur_full <- merge_phyloseq(mothur_phyloseq_object, metafile, refseq(OTU_seqs))

# code for Figure 4: Alpha diversity (including singletons)
erich_mothur <- estimate_richness(mothur_full, measures = c(
  "Observed",
  "Chao1",
  "ACE",
  "Shannon",
  "Simpson",
  "InvSimpson",
  "Fisher"))
erich_mothur_meta <- cbind(erich_mothur, sample_data(mothur_full)[,c(1:7)])

# reorder and rename factor levels for plotting, adjust displayed label names
erich_mothur_meta$habitat <- relevel(erich_mothur_meta$habitat, "water")
erich_mothur_meta$nucleic_acid <- relevel(erich_mothur_meta$nucleic_acid, "dna")
labels_nucleic_acid <- c("16S rRNA gene", "16S rRNA")
labels_habitat <- c("Free-living", "Biofilm")
levels(erich_mothur_meta$habitat) <- labels_habitat
levels(erich_mothur_meta$nucleic_acid) <- labels_nucleic_acid

# plotting Alpha diversity
shannon_plot <- ggplot(erich_mothur_meta, aes(x = new_day, y = Shannon, colour = treatment)) +
  geom_point(alpha = 0.8, size = 4) +
  geom_vline(aes(xintercept = 1), linetype = "dashed", size = 1.2) +
  stat_summary(aes(colour = treatment), fun.y = "mean", geom = "line", 
    alpha = 0.75, size = 2) +
  scale_colour_manual(values = c("glyph" = "black", "control" = "grey50"),
	name = "Microcosm  ", breaks = c("glyph", "control"), 
    labels = c("Treatment", "Control")) +
  coord_cartesian(ylim = c(1, 3)) +
  theme_bw() +
  theme(axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 13),
    panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    strip.text.x = element_text(size = 15, face = "bold")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  labs(x = "Days", y = "Shannon index") +
  facet_wrap(~ habitat + nucleic_acid)

ggsave(shannon_plot, file = paste(plot_path, "Figure_4_Shannon_DNA_RNA.png",
  sep = ""), height = 10, width = 14)

# Figure 2: treatment community overview bar plot
# remove singletons
mothur_1 <- filter_taxa(mothur_full, function (x) {sum(x > 1) >= 1}, prune = TRUE)

# transform into relative abundance, displayed in percentage!
mothur_full_ra <- transform_sample_counts(mothur_full, function(x) {(x / sum(x)) * 100})

# remove low abundant OTUs (you may decrease the threshold, but it will increase melting time)
mothur_ra_0.01 <- filter_taxa(mothur_full_ra, function (x) {sum(x > 0.01) >= 1}, prune = TRUE)

# melt into long format for plotting
mothur_ra_melt <- psmelt(mothur_ra_0.01)
mothur_1_melt <- psmelt(mothur_1)

# aggregate all columns except for "parallels" to 
# calculate the mean abundance of technical replicates
mothur_ra_melt_mean <- aggregate(Abundance ~ OTU + time + days + new_day
  + treatment + nucleic_acid + habitat + disturbance + glyphosate + glyphosate_gone 
  + condition_diversity + kingdom + phylum + class + order + family + genus + otu_id,
  data = mothur_ra_melt, mean)
                                
# get data per OTU, setting threshold for samples and clusters
community_subset <- droplevels(subset(mothur_ra_melt_mean, days > 40
  & Abundance > 0.15 & habitat == "water" & treatment == "glyph"))
                            
# check required number of colours per order and number of classes
length(levels(droplevels(community_subset$class)))
length(levels(droplevels(community_subset$order)))

# sort orders for plotting based on phylogenetic 
# classes (Alphaproteos, Gammaproteos and Bacteriodetes) 
community_subset$order <- factor(community_subset$order, levels = c(
  # alphaproteos
  "Caulobacterales",
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
fill_values2 <- c("Alteromonadales" = "#e6194B",
  "Betaproteobacteriales" = "#3cb44b",
  "Caulobacterales" = "#ffe119",
  "Chitinophagales" = "#4363d8",
  "Flavobacteriales" = "#f58231",
  "Sneathiellales" = "#42d4f4",
  "Parvibaculales" = "#f032e6",
  "Pseudomonadales" = "#fabebe",
  "Rhizobiales" = "#469990",
  "Rhodobacterales" = "#000000",
  "Rhodospirillales" = "#9A6324",
  "Sphingobacteriales" = "#fffac8",
  "Sphingomonadales" = "#800000",
  "Thalassobaculales" = "#a9a9a9")

# plotting all selected clusters in bar plot ordered by class
# and displaying orders over time for DNA and RNA
community_plot <- ggplot(community_subset, aes(x = new_day, group = order)) +
  scale_fill_manual(breaks = levels(community_subset$order), values = fill_values2) +
  geom_bar(data = subset(community_subset, nucleic_acid == "dna" &
    treatment == "glyph"), aes(x = new_day - 0.5, y = Abundance), fill = "black",
    width = 0.9, stat = "sum") +
  geom_bar(data = subset(community_subset, nucleic_acid == "dna" & treatment == "glyph"),
    aes(x = new_day - 0.5, y = Abundance, fill = order), width = 0.6, stat = "identity") +
  geom_bar(data = subset(community_subset, nucleic_acid == "cdna" & treatment == "glyph"),
    aes(x = new_day + 0.5, y = Abundance), fill = "black", width = 0.9, stat = "sum") +
  geom_bar(data = subset(community_subset, nucleic_acid == "cdna" & treatment == "glyph"),
    aes(x = new_day + 0.5, y = Abundance, fill = order), width = 0.6, stat = "identity") +
  geom_vline(data = subset(community_subset, treatment == "glyph"), aes(xintercept = 1.5),
    linetype = "dashed", size = 1.2) +
  guides(colour = FALSE, size = FALSE, width = FALSE, fill = guide_legend(ncol = 1,
    keyheight = 1.5, label.theme = element_text(size = 15, face = "italic",
    angle = 0), (title = NULL))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text = element_text(size = 17),
    axis.title = element_text(size = 20, face = "bold"),
    legend.background = element_rect(fill = "grey90", linetype = "solid"),
    panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5)) +
  labs(x = "Days", y = "Relative abundance [%]") +
  annotate("text", x = -27, y = 90, label = "a)", color = "black", size = 6,
    angle = 0, fontface = "bold") +
  annotate("text", x = -22.5, y = 90, label = "b)", color = "black", size = 6,
    angle = 0, fontface = "bold")

ggsave(community_plot, file = paste(plot_path, 
  "Figure_4_relative_community_overview.png", sep = ""), width = 16, height = 8)

# Figure 5 and Supplement 5: OTU abundance plots
# define subset function for specific phyloseq-object
get_current_otu_data <- function(x) {
	subset(mothur_ra_melt, OTU == x)
}

# list of OTUs mentioned in paper and supplement
OTU_list <- c("Otu000007",
  "Otu000011",
  "Otu000018",
  "Otu000025",
  "Otu000032",
  "Otu000036",
  "Otu000037",
  "Otu000038",
  "Otu000023",
  "Otu000046",
  "Otu000049",
  "Otu000056",
  "Otu000058",
  "Otu000059",
  "Otu000070",
  "Otu000072",
  "Otu000078",
  "Otu000094",
  "Otu000109",
  "Otu000129",
  "Otu000139",
  "Otu000176",
  "Otu000191",
  "Otu000320",
  "Otu000098",
  "Otu000042",
  "Otu000044",
  "Otu000006",
  "Otu000001")

# rename for plotting
strip_text_habitat <- c("Biofilm", "Free-living")

# run a for loop to plot each OTU in list with own title and file name
for (i in OTU_list) {
    current_otu_data <- get_current_otu_data(i)
    print(paste("OTU is", i))

species_title <- unique(paste(current_otu_data$family, current_otu_data$genus,
  current_otu_data$OTU, sep = " "))

levels(current_otu_data$habitat) <- strip_text_habitat

current_plot <- ggplot(data = current_otu_data, aes(x = days - 69,
  y = Abundance, group = nucleic_acid, lty = nucleic_acid)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1.2) +
  geom_point(data = subset(current_otu_data, treatment == "control"),
    aes(colour = treatment), alpha = 1) +
  stat_summary(data = subset(current_otu_data, treatment == "control"),
    aes(colour = treatment), fun.y = "mean", geom = "line", size = 2, alpha = 1) +
  stat_summary(data = subset(current_otu_data, treatment == "glyph"),
    aes(colour = treatment), fun.y = "mean", geom = "line", size = 2) +
  geom_point(data = subset(current_otu_data, treatment == "glyph"),
    aes(colour = treatment)) +
  scale_linetype_manual(values = c("dna" = 1, "cdna" = 6), name = "Nucleic acid  ",
    breaks = c("cdna", "dna"), labels = c("16S rRNA", "16S rRNA gene")) +
  scale_colour_manual(values = c("glyph" = "black", "control" = "grey50"),
    name = "Microcosm  ", breaks = c("glyph", "control"), labels = c("Treatment",
    "Control")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_bw() +
  ggtitle(species_title) +
  theme(axis.text = element_text(size = 18),
    axis.title = element_text(size = 20, face = "bold"),
    panel.grid.major = element_line(colour = NA, size = 0.2),
    panel.grid.minor = element_line(colour = NA, size = 0.5),
    strip.text.x = element_text(size = 15, face = "bold")) +
  labs(x = "Days", y = "Relative abundance [%]") +
  facet_wrap(~ habitat, scales = "free")
  ggsave(current_plot, file = paste(plot_path, species_title,".png", sep = ""),
    width = 13, height = 7)
  print(current_plot)
}

# Figure 1: Total cell counts and glyphosate concentration
cell_counts_glyph_0 <- subset(cell_counts_glyph, new_day >= -7)

cell_counts_glyph_plot <- ggplot(cell_counts_glyph_0, aes(x = new_day,
  colour = treatment, linetype = treatment, group = treatment)) +
  geom_errorbar(aes(ymin = cells_ml - cells_se, ymax = cells_ml + cells_se),
    linetype = 1, width = 2, size = 1.2, alpha = 0.7) +
  geom_errorbar(aes(ymin = (glyph_micromol - glyph_se_micromol) * 400000,
    ymax = (glyph_micromol + glyph_se_micromol) * 400000), linetype = 1, 
    width = 1, size = 1.0, alpha = 0.6) +
  geom_line(aes(y = cells_ml, group = treatment, colour = treatment), linetype = 1, 
    size = 2, alpha = 0.8) +
  geom_point(aes(y = glyph_theor_micromol * 400000, shape = "glyph"), alpha = 0.5, 
    size = 5) +
  geom_point(aes(y = glyph_micromol * 400000, shape = "glyph_deg"), alpha = 0.7, 
    size = 4) +
  geom_vline(aes(xintercept = 0), linetype = "dashed", size = 1.5, alpha = 0.5) +
  scale_y_continuous(label =  function(x) {ifelse(x == 0, "0", parse(text = 
    gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))},
    sec.axis = sec_axis(~ . / 400000, name = expression(paste(
      "Glyphosate concentration  ", bgroup("[",µM,"]"))))) +
  scale_colour_manual(values = c("control" = "grey70", "glyph" = "black"),
    name = "Microcosm", breaks = c("control", "glyph"), 
    labels = c("Control", "Treatment")) +
  scale_shape_manual(values = c("glyph" = 2, "glyph_deg" = 17), 
    name = "Glyphosate decrease by", breaks = c("glyph", "glyph_deg"),
    labels = c("calculated\ndilution", "measured\nconcentration")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme_bw() +
	theme(panel.grid.major = element_line(colour = NA, size = 0.2),
      panel.grid.minor = element_line(colour = NA, size = 0.5),
      axis.title = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(angle = 90, vjust = 1),
      axis.text = element_text(size = 18),
      legend.title=element_text(size = 14),
      legend.text=element_text(size = 12)) +
  labs(x = "Days", y = expression(paste("Total cell count  ",
    bgroup("[",cells~mL^{-1},"]"))))

ggsave(cell_counts_glyph_plot, file = paste(plot_path, 
  "Figure_1_cellcounts_glyph.png", sep = ""), width = 14, height = 10)

# Figure 3: NMDS plots
# exclude OTUs with less than 3 reads and transform to relative abundance
phyloseq_for_nmds <- filter_taxa(mothur_full, function (x) {sum(x > 2) >= 1}, prune = TRUE)
phyloseq_for_nmds_rel_abund <- transform_sample_counts(phyloseq_for_nmds, function(x) {(x / sum(x)) * 100})

# arguments for the subsetting function
phyloseq_object <- phyloseq_for_nmds_rel_abund
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
threshold <- 0
days <- 40

# define a function to obtain sample subsets from the phyloseq object 
# per combination of habitat, nucleic acid, days and minimum required reads per OTU
get_sample_subsets <- function(phyloseq_object, nucleic_acid, habitat, days, threshold) {
  sample_subset <- sample_data(phyloseq_object)[ which(sample_data(phyloseq_object)$nucleic_acid == nucleic_acid &
    sample_data(phyloseq_object)$habitat == habitat & sample_data(phyloseq_object)$days > days),]
  phyloseq_subset <- merge_phyloseq(tax_table(phyloseq_object),
    otu_table(phyloseq_object),
    refseq(phyloseq_object),
    sample_subset)
  phyloseq_subset2 <- filter_taxa(phyloseq_subset, function (x) {sum(x > threshold) >= 1 }, prune = TRUE)
  return(phyloseq_subset2)
}

# here we pass the arguments for subsetting over two for loops
# to create all possible combinations of habitat, nucleic acid etc.
# the subsets are stored within a list, which has to be empty before running the loops 
sample_subset_list <- list()
if(length(sample_subset_list) == 0) {
  for (acid in acids) {
	for (habitat in habitats) {
	  print(paste0("nucleic_acid is ", acid, " and habitat is ", habitat))
	  tmp <- get_sample_subsets(phyloseq_object = phyloseq_object,
		nucleic_acid = acid, habitat = habitat, threshold = threshold, days = days)
	  sample_data(tmp)$days <- as.factor(sample_data(tmp)$days)
	  sample_data(tmp)$new_day <- as.factor(sample_data(tmp)$new_day)
	  sample_subset_list[[paste(habitat, acid, "min reads per OTU", threshold,
		sep = " ")]] <- tmp
    }
  }
  print(sample_subset_list)
} else {
  print("list is not empty, abort to prevend appending...")
}

# create a list where the distance metrics for the sample subsets are stored
ordination_nmds <- list()
ordination_nmds <- lapply(sample_subset_list, ordinate, method = "NMDS",
  dist = "bray", try = 100, autotransform = TRUE)

# test hypothesis that microcosms are different
set.seed(1)

# Calculate bray curtis distance matrix and collect sample data in lists
sample_distances <- list()
sample_data_list <- list()
sample_distances <- lapply(sample_subset_list, phyloseq::distance, method = "bray")
sample_data_list <- lapply(sample_subset_list, function(x) data.frame(sample_data(x)))

# Adonis test against hardcoded "treatment", instead passing "treatment" as argument did not work
adonis_results <- list()
adonis_results <- mapply(function(distance_matrix, sample_data) {
  adonis(distance_matrix ~ treatment, data = data.frame(sample_data))
}, distance_matrix = sample_distances, sample_data = sample_data_list, SIMPLIFY = FALSE)

# Homogeneity of dispersion test 
# tests if dispersion might be reason for adonis results, should confirm null hypothesis
beta_list <- list()
beta_list <- mapply(function(distance_matrix, sample_data) {
  betadisper(distance_matrix, sample_data$treatment)
}, distance_matrix = sample_distances, sample_data = sample_data_list, SIMPLIFY = FALSE)
dispersion_list <- list()
dispersion_list <- lapply(beta_list, permutest)

# generate NMDS plots with distance information from the ordination list
nmds_ordination_plots <- list()
counter <- 0
if(length(nmds_ordination_plots) == 0 & all.equal(counter, 0)) {
  nmds_ordination_plots <- mapply(function(x,y) {
    counter <<- counter + 1
    plot_ordination(x, y, type = "sample", color = "days", shape = "treatment") +
      geom_polygon(aes(fill = disturbance), alpha = 0.5, size = 0.01) +
      geom_point(aes(colour = treatment), colour = "black", size = 4.5, alpha = 0.7) +
      scale_shape_manual(values = c("glyph" = 16, "control" = 17), name = "Microcosm  ",
        breaks = c("glyph", "control"), labels = c("Treatment", "Control")) +
      scale_fill_manual(values = c("high" = "red", "low" = "orange", "none" = "green"),
        name = "Present glyphosate\nconcentration", breaks = c("high", "low", "none"),
        labels = c("High", "Low", "None")) +
      guides(color = FALSE) +
      coord_cartesian(ylim = c(-0.77, 0.95), xlim = c(-0.9, 0.7)) +
      #ggtitle(names(sample_subset_list)[counter]) +
      geom_text(aes(label = new_day), colour = "white", size = 2.5) +
      theme_bw() +
      theme(panel.grid.major = element_line(colour = NA, size = 0.2),
        panel.grid.minor = element_line(colour = NA, size = 0.5),
        axis.text = element_text(size = 18),
        #axis.title = element_text(size = 20, face = "bold"),
        axis.title = element_blank(),
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 13))
}, x = sample_subset_list, y = ordination_nmds, SIMPLIFY = FALSE)
} else {
  print(paste("list is not empty, or counter not 0 (counter is", counter,"), 
    abort to prevend appending..."))
}

# plot the NMDS ordinations in a 2 x 2 array
do.call("grid.arrange", c(nmds_ordination_plots[c(1, 3, 2, 4)], nrow = 2))

# write the NMDS ordinations array to an object
NMDS_array <- do.call("arrangeGrob", c(nmds_ordination_plots[c(1, 3, 2, 4)], nrow = 2))
ggsave(NMDS_array, file = paste(plot_path, "Figure_3_NMDS.png", sep = ""),
  height = 10, width = 10)

# Table 1: differentially abundant OTUs tested by DESeq2
# test variable is not allowed to contain NA
mothur_deseq <- subset_samples(mothur_full, !(is.na(condition)))

# these are the parameters passed to function
phyloseq_object<- mothur_deseq
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
treatments <- c("glyph", "control")
threshold <- 0

# this is the function we call to split our data into different subsets
get_sample_subsets_deseq <- function(phyloseq_object, nucleic_acid, habitat, treatment) {
  sample_subset <- sample_data(phyloseq_object)[ which(sample_data(phyloseq_object)$nucleic_acid == nucleic_acid &
    sample_data(phyloseq_object)$habitat == habitat &
    sample_data(phyloseq_object)$treatment == treatment),]
  phyloseq_subset <- merge_phyloseq(tax_table(phyloseq_object),
    otu_table(phyloseq_object), sample_subset)
  phyloseq_subset2 <- filter_taxa(phyloseq_subset, function (x) {sum(x > 0) >= 1}, prune = TRUE)
  return(phyloseq_subset2)
}

# here we pass the arguments for subsetting over three for loops
# to create all possible combinations of habitat, nucleic acid and microcosm 
deseq_subsets <- list()
if(length(deseq_subsets) == 0) {
  for (treatment in treatments) {
    for (acid in acids) {
      for (habitat in habitats) {
        print(paste0("nucleic_acid is ", acid, " and habitat is ",
          habitat, " and treatment is ", treatment))
        tmp <- get_sample_subsets_deseq(phyloseq_object = phyloseq_object,
          nucleic_acid = acid, habitat = habitat, treatment = treatment)
        sample_data(tmp)$days <- as.factor(sample_data(tmp)$days)
        sample_data(tmp)$new_day <- as.factor(sample_data(tmp)$new_day)
        deseq_subsets[[paste(habitat, treatment, acid, sep = "_")]] <- tmp
      }
    }
  }
print(deseq_subsets)
} else {
  print("list is not empty, abort to prevend appending...")
}

# run DESeq2 on the sample subsets to test for differentially abundant OTUs
deseq_tests <- list()
counter <- 0
if(length(deseq_tests) == 0 & all.equal(counter, 0)) {
  deseq_tests <- lapply(deseq_subsets, function(deseqs) {
    counter  <<- counter + 1
    tmp = phyloseq_to_deseq2(deseqs, ~ condition)
    tmp$condition <- relevel(tmp$condition, ref = "untreated")
    tmp_dds = DESeq(tmp, test = "Wald", fitType = "parametric")
  })
} else {
  print(paste("list is not empty, or counter not 0 (counter is", counter,"), 
    abort to prevend appending..."))
}

# combine the DESeq2-test results with the corresponding taxonomy
# and sort detected OTUs by log fold change
sigtabs_list <- list()
if(length(sigtabs_list) == 0) {
sigtabs_list <- mapply(function(dds, phyloseq_object) {res = results(dds, cooksCutoff = FALSE)
  alpha = 0.01
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"),
    as(tax_table(phyloseq_object)[rownames(sigtab), ], "matrix"))
  print(head(sigtab))
  return(sigtab)
  }, dds = deseq_tests, phyloseq_object = deseq_subsets, SIMPLIFY = FALSE)
} else {
  print(paste("list is not empty, abort to prevent appending..."))
}
sigs_ordered <- lapply(sigtabs_list, function(x) x[order(x$log2FoldChange),])

# create a vector of all identified OTUs in sigtabs_list (can be used for plotting abundances!)
deseq_otus <- row.names(sigtabs_list[[1]])
for (i in 2:8) {deseq_otus <- unique(append(deseq_otus, row.names(sigtabs_list[[i]])))}

# Supplement 3: Venn Diagrams for genera and OTU distribution
# define function to plot Venn diagram with 4 categories, here biofilm vs water column
fourway.Venn <- function(A, B, C, D, cat.names = c("Water\nDNA",
  "Biofilm\nDNA", "Water\nRNA", "Biofilm\nRNA")) {
    grid.newpage()
    # calculate the values for the different intersections and areas
    area1 <- length(A)
    area2 <- length(B)
    area3 <- length(C)
    area4 <- length(D)
    n12 <- length(Reduce(intersect, list(A, B)))
    n13 <- length(Reduce(intersect, list(A, C)))
    n14 <- length(Reduce(intersect, list(A, D)))
    n23 <- length(Reduce(intersect, list(B, C)))
    n24 <- length(Reduce(intersect, list(B, D)))
    n34 <- length(Reduce(intersect, list(C, D)))
    n123 <- length(Reduce(intersect, list(A, B, C)))
    n124 <- length(Reduce(intersect, list(A, B, D)))
    n134 <- length(Reduce(intersect, list(A, C, D)))
    n234 <- length(Reduce(intersect, list(B, C, D)))
    n1234 <- length(Reduce(intersect, list(A, B, C, D)))

  venn.plot <- draw.quad.venn(
    area1 = area1,
    area2 = area2,
    area3 = area3,
    area4 = area4,
    n12 = n12,
    n13 = n13,
    n14 = n14,
    n23 = n23,
    n24 = n24,
    n34 = n34,
    n123 = n123,
    n124 = n124,
    n134 = n134,
    n234 = n234,
    n1234 = n1234,
    category = cat.names,
    cat.pos = c(0, 180, 0, 200),
    fill = c("blue", "red", "green", "yellow"),
    alpha = 0.3,
    lty = "blank",
    cex = 2,
    cat.cex = 2,
    cat.col = c("blue", "red", "green", "black"))
  grid.draw(venn.plot)
}

# factorize OTUs to count them
mothur_ra_melt$OTU <- as.factor(mothur_ra_melt$OTU)

# these are the arguments for the Venn subsetting function
melt_object<- mothur_ra_melt
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
abundance <- 0.05	
    
# function to generate subsets from a long format dataframe
# this is the function we call to split our data into different subsets
subset_melt <- function(melt_object, nucleic_acid, habitat, abundance) {
  melt_subsetted <- melt_object[which(melt_object$nucleic_acid == nucleic_acid &
    melt_object$habitat == habitat &  melt_object$Abundance > abundance),]
    #return(melt_subsetted)
}

# here we pass the arguments for subsetting over three for loops
# to create all possible combinations of habitat, nucleic acid and microcosm 
Venn_subsets <- list()
if(length(Venn_subsets) == 0) {
  for (current_abundance in abundance) {
    for (acid in acids) {
      for (habitat in habitats) {
        print(paste0("nucleic_acid is ", acid, " and habitat is ", habitat, 
          " and threshold is ", current_abundance))
        tmp <- subset_melt(melt_object = melt_object, nucleic_acid = acid,
          habitat = habitat, abundance = current_abundance)
        Venn_subsets[[paste(habitat, current_abundance, acid, sep = "_")]] <- tmp
      }
    }
  }
} else {
  print("list is not empty, abort to prevend appending...")
}

# retrieve the unique genera per subset from the list
water_dna_unique_genera <- Venn_subsets[[1]][which(!duplicated(Venn_subsets[[1]]$genus)),]
biofilm_dna_unique_genera <- Venn_subsets[[2]][which(!duplicated(Venn_subsets[[2]]$genus)),]
water_cdna_unique_genera <- Venn_subsets[[3]][which(!duplicated(Venn_subsets[[3]]$genus)),]
biofilm_cdna_unique_genera <- Venn_subsets[[4]][which(!duplicated(Venn_subsets[[4]]$genus)),]

# plot Venn diagram
fourway.Venn(water_dna_unique_genera$genus, biofilm_dna_unique_genera$genus,
  water_cdna_unique_genera$genus, biofilm_cdna_unique_genera$genus)
dev.copy(png, paste(plot_path, "Supplement_4wayVenn_nucleic_acids_genus_0.05.png"))
dev.off()

# SI 3: values for sequence length and OTU/genera distribution
# library sizes are returned using
sample_sums(mothur_full)

# how many OTUs belong to which genus?
genus_distribution <- aggregate(Abundance ~ OTU + genus, data = mothur_1_melt, max)

# separated by habitat and nucleic acid
genus_distribution2 <- aggregate(Abundance ~ OTU + habitat + genus + nucleic_acid,
  data = mothur_1_melt, max)

# OTUs per genus ordered by amount of OTUs
otu_per_genus <- as.data.frame(table(genus_distribution$genus))
otu_per_genus[order(otu_per_genus$Freq),]
nrow(otu_per_genus)

# these are the arguments for the subsetting function
phyloseq_object <- mothur_full
acids <- c("dna", "cdna")
habitats <- c("water", "biofilm")
threshold <- 1
after_day <- 43

# in this list we store the different sample subsets, generated by the for loops
sample_subset_list <- list()
if(length(sample_subset_list) == 0) {
  for (each_day in after_day) {
    for (acid in acids) {
      for (habitat in habitats) {
        print(paste0("nucleic_acid is ", acid, " and habitat is ", habitat, 
          " and first day is ", each_day))
        tmp <- get_sample_subsets(phyloseq_object = phyloseq_object,
          nucleic_acid = acid, habitat = habitat, days = each_day, threshold = threshold)
        sample_data(tmp)$days <- as.factor(sample_data(tmp)$days)
        sample_data(tmp)$new_day <- as.factor(sample_data(tmp)$new_day)
        sample_subset_list[[paste(habitat, "after day", each_day, acid,
          "min reads per OTU", threshold, sep = " ")]] <- tmp
      }
    }
  }
print(sample_subset_list)
} else {
  print("list is not empty, abort to prevend appending...")
}

# the distribution of sequence length can be addressed using
table(width(refseq(sample_subset_list[[1]])))

sessionInfo()