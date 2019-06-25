#!/usr/bin/env Rscript
library(phyloseq)
library(Biostrings)
library(tidyverse)

setwd("/data/projects/glyphosate/reads/mothur_processed/")

plot_path <- "./plots/"
if(!dir.exists(plot_path)) {
  dir.create(plot_path)
}

# import files
our_shared_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.shared"
our_cons_taxonomy_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.0.02.cons.taxonomy"

metafile2 <- read.delim("/data/projects/glyphosate/analysis/metadata/metafile2.tsv", 
  row.names = 1, header = TRUE, na.strings = "")
metafile2 <- sample_data(metafile2)

OTU_seqs <- readDNAStringSet(file = "OTU_reps_fasta_002.fasta", format = "fasta", 
  nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

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
                                    
# phyloseq_object variatons for different analyses
mothur_full <- merge_phyloseq(mothur_phyloseq_object, metafile2, refseq(OTU_seqs))
mothur_1 <- filter_taxa(mothur_full, function (x) {sum(x > 1) >= 1}, prune = TRUE)
mothur_full_ra <- transform_sample_counts(mothur_full, function(x) {(x / sum(x)) * 100})
mothur_ra_0.01 <- filter_taxa(mothur_full_ra, function (x) {sum(x > 0.01) >= 1}, prune = TRUE)
mothur_ra_melt <- psmelt(mothur_ra_0.01)
mothur_ra_melt_mean <- aggregate(Abundance ~ OTU + time + days + new_day
  + treatment + nucleic_acid + habitat + disturbance + glyphosate + glyphosate_gone 
  + condition_diversity + kingdom + phylum + class + order + family + genus + otu_id,
  data = mothur_ra_melt, mean)

# store phyloseq objects for other scripts
save.image("glyphosate_mothur_in_phyloseq.RData")

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