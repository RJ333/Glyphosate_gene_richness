#!/usr/bin/env Rscript
library(phyloseq)
library(DESeq2)

setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("glyphosate_mothur_in_phyloseq.RData")

# Table 1: differentially abundant OTUs tested by DESeq2
# test variable is not allowed to contain NA
mothur_deseq <- subset_samples(mothur_full, !(is.na(condition)))

# these are the parameters passed to function
phyloseq_object <- mothur_deseq
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
