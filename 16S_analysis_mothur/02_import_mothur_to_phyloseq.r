# Set the working dir with mothur files in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")

# define required packages
library("knitr")
library("BiocStyle")

.cran_packages <- c("ggplot2", "gridExtra")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
   install.packages(.cran_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
   source("http://bioconductor.org/biocLite.R")
   biocLite(.bioc_packages[!.inst], ask = F)
}

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# when you combine different objects such as tax and otu table,
# phyloseq performs an inner joint!

# import mothur output into phyloseq
mothur_ps <- import_mothur(mothur_list_file = NULL, 
  mothur_group_file = NULL,
  mothur_tree_file = NULL, 
  cutoff = NULL, 
  mothur_shared_file = "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.shared",
  mothur_constaxonomy_file = "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.0.02.cons.taxonomy", 
  parseFunction = parse_taxonomy_default)

rank_names(mothur_ps)
  
# adjust taxonomy header
colnames(tax_table(mothur_ps)) <- c(k = "Kingdom", p = "Phylum", c= "Class", 
									o = "Order", f = "Family", g = "Genus" )
									
# when you combine different objects such as otu and meta table,
# phyloseq performs an inner joint!

# read meta data, turn into phyloseq object, merge with existing ps object									
metafile <- read.delim("/data/projects/glyphosate/analysis/metadata/all_samples_with_meta.tsv", row.names = 1, header = TRUE)
metafile2 <- sample_data(metafile)

# read OTU representative sequences, for generation of file check gitlab #59
OTU_seqs <- readDNAStringSet(file = "OTU_reps.fasta", format = "fasta",
               nrec = -1L, skip = 0L, seek.first.rec = FALSE, use.names = TRUE)

mothur_ps2 <- merge_phyloseq(mothur_ps, metafile2, refseq(OTU_seqs))

mothur_ps3 <- filter_taxa(mothur_ps, function (x) {sum(x > 0) > 1}, prune=TRUE)

# this changes the header from the actual sequence to Seq_001, Seq_002 etc
taxa_names(mothur_ps)

# taxa are already "Otuxxxxxx", no adjustment needed
# n_seqs <- seq(ntaxa(mothur_ps))
# len_n_seqs <- nchar(max(n_seqs))
# taxa_names(mothur_ps) <- paste("Seq", formatC(n_seqs, 
# 											width = len_n_seqs, 
# 											flag = "0"), sep = "_")
# taxa_names(mothur_ps)


# generate a tree
seqs <- getSequences(OTU_seqs)
#names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)

# phangorn package: we first construct a neighbor-joining tree,
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data = phang.align)

# then fit a GTR+G+I maximum likelihood tree using the neighbor-joining tree as a starting point.
fitGTR <- update(fit, k = 4, inv = 0.2)
# the following step takes longer
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

# generate taxonomy from tax table as df
wholetax <- do.call(paste, c(as.data.frame(tax_table(mothur_ps))
                  [c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 
				  sep = "__"))
				  
# generate otu table
otu_export <- as.data.frame(otu_table(mothur_ps))
tmp <- names(otu_export)

for(i in 1:length(tmp)){
names(tmp)[i] = paste(wholetax[i], tmp[i], sep = "__")
}

names(otu_export) <- names(tmp)


# bind meta data

# generate tree


