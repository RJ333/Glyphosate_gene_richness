# this workflow describes how to import data from mothur into phyloseq. For full
# usage, you need a shared, a constaxonomy, a metadata and a otu_rep.fasta file
# a tree is optional, but will be generated within R anyway.


# Set the working dir with mothur files in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("mothur_glyph.RData")

save.image("mothur_glyph.RData")
# if you load a workspace, you then only need to load the packages

# define required packages
# library("knitr")
# library("BiocStyle")

# .cran_packages <- c("ggplot2", "gridExtra")
# .bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn")

# .inst <- .cran_packages %in% installed.packages()
# if(any(!.inst)) {
   # install.packages(.cran_packages[!.inst])
# }

# .inst <- .bioc_packages %in% installed.packages()
# if(any(!.inst)) {
   # source("http://bioconductor.org/biocLite.R")
   # biocLite(.bioc_packages[!.inst], ask = F)
# }

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
  
# adjust taxonomy header
tax_table(mothur_ps) <- cbind(tax_table(mothur_ps), 
    rownames(tax_table(mothur_ps)))
colnames(tax_table(mothur_ps)) <- c("kingdom", 
									"phylum", 
									"class", 
									"order", 
									"family", 
									"genus", 
									"otu_id")
rank_names(mothur_ps)
											
# when you combine different objects such as otu and meta table,
# phyloseq performs an inner joint!

# read meta data, turn into phyloseq object, merge with existing ps object									
metafile <- read.delim("/data/projects/glyphosate/analysis/metadata/all_samples_with_meta_cond.tsv", 
						row.names = 1, 
						header = TRUE,
						na.strings = "")
metafile2 <- sample_data(metafile)


# read OTU representative sequences, for generation of file check gitlab #59
OTU_seqs <- readDNAStringSet(file = "OTU_reps.fasta", 
							  format = "fasta",
							  nrec = -1L, 
							  skip = 0L, 
							  seek.first.rec = FALSE, 
							  use.names = TRUE)
# add meta data and OTU representative seqs to phyloseq object
mothur_ps2 <- merge_phyloseq(mothur_ps, metafile2, refseq(OTU_seqs))
# remove OTUs with less than 3 reads
mothur_ps3 <- filter_taxa(mothur_ps2, function (x) {sum(x > 0) > 2}, prune = TRUE)
# transform into relative abundance, displayed in percentage!
mothur_ps3_ra <- transform_sample_counts(mothur_ps3, function(x){(x / sum(x)) * 100})
mothur_ra_melt <- psmelt(mothur_ps3_ra)
# add absolute abundance (rel abundance * total cell counts)
mothur_ra_melt$abs_Abundance <- (mothur_ra_melt$Abundance * mothur_ra_melt$cell_counts)/100
# turn OTUs from char to factor
mothur_ra_melt$OTU <- as.factor(mothur_ra_melt$OTU)

# depending on the plot, we need the parallels separately or averaged
# need to remove "Sample" to average the parallels
mothur_ra_melt_mean <- aggregate(Abundance ~ OTU + time + days + new_day
								+ treatment + nucleic_acid + habitat + disturbance 
								+ cell_counts + glyphosate + glyphosate_gone 
								+ kingdom + phylum + class + order + family + genus + otu_id, 
								data = mothur_ra_melt, 
								mean)

# generate a tree to add to the phyloseq object
# we use filtered sequences (OTU > 2 reads)
seqs_pruned  <- refseq(mothur_ps3)
alignment <- AlignSeqs(seqs_pruned, anchor = NA, processors = NULL)

# phangorn package: we first construct a neighbor-joining tree,
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
# this takes longer
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data = phang.align, multicore = TRUE)

### below: running at the moment, not tested with this dataset, but worked before


# then fit a GTR+G+I maximum likelihood tree 
# using the neighbor-joining tree as a starting point.
fitGTR <- update(fit, k = 4, inv = 0.2)
# the following step takes longer
fitGTR <- optim.pml(fitGTR, 
					model = "GTR", 
					optInv = TRUE, 
					optGamma = TRUE,
                    rearrangement = "stochastic", 
					multicore = TRUE, 
					control = pml.control(trace = 1))
detach("package:phangorn", unload = TRUE)

# add the generated tree to phyloseq
mothur_ps4_ra <- merge_phyloseq(mothur_ps3_ra, phy_tree(fitGTR$tree))