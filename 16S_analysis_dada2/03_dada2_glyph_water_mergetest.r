# in this script I tested how to combine objects from different workspaces
# namely the processed dada2 objects, which were treated separately 


# define working directory to story RData image
setwd("/data/projects/glyphosate/reads/dada2_processed/water_dna")

# load first workspace representing first amplicon sequencing run
load(file = "dada2_water_dna.RData")

# this will be the workspace for the combined dada2 objects
save.image(file = "../dada2_water_comb.RData")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)

# rename respective object from workspace to differentiate from next workspace
dna_water_nochim <- seqtab2.nochim

# add other work spaces to temporary environment
tmp.env <- new.env()
load("../water_cdna/dada2_water_cdna.RData", envir = tmp.env)

# get the objects you need into your globalenv and remove tmp.env
cdna_water_nochim <- get("seqtab2.nochim", pos = tmp.env)
rm(tmp.env)

# merge the sequencetable from two differen workspaces into one table
water_seqtable <- mergeSequenceTables(table1 = dna_water_nochim, 
									  table2 = cdna_water_nochim,
									  repeats = "error", 
									  orderBy = "abundance")
  
# assign taxonomy with specially prepared silva database
# NOTE: it should be possible to use the mothur-like amplicon-trimmed database
# to reduce memory amount. Not tested

ref_fasta <- "/data/db/silva_nr_v132_train_set.fa.gz"
taxtab <- assignTaxonomy(water_seqtable, 
						 refFasta = ref_fasta, 
						 multithread = TRUE)
colnames(taxtab) <- c("Kingdom", 
					  "Phylum", 
					  "Class", 
					  "Order", 
					  "Family", 
					  "Genus")
seqs <- getSequences(water_seqtable)

# the tree should probably be calculated later, with singletons removed
# as this will reduce the computation time

names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor = NA)

# phangorn package: we first construct a neighbor-joining tree,
phang.align <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data = phang.align)

# then fit a GTR+G+I maximum likelihood tree using the neighbor-joining tree as a starting point.
fitGTR <- update(fit, k = 4, inv = 0.2)
# the following step takes longer
fitGTR <- optim.pml(fitGTR, 
					model = "GTR", 
					optInv = TRUE, 
					optGamma = TRUE,
                    rearrangement = "stochastic", 
					control = pml.control(trace = 0))
detach("package:phangorn", unload = TRUE)

# get this into phyloseq as phyloseq object (or do tree later)
ps_dada <- phyloseq(tax_table(taxtab), 
					otu_table(water_seqtable, 
							  taxa_are_rows = FALSE), 
					#phy_tree(fitGTR$tree),
					refseq(seqs))


					
### dada2 uses actual sequences as OTU headers
# test to change the sequence headers					
# adjust phyloseq object name to be not overwritten
ps_new <- ps_dada

# this changes the header from the actual sequence to Seq_001, Seq_002 etc
taxa_names(ps_new)
n_seqs <- seq(ntaxa(ps_new))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps_new) <- paste("Seq", 
							formatC(n_seqs, 
									width = len_n_seqs, 
									flag = "0"), 
							sep = "_")
taxa_names(ps_new)

# generate taxonomy path to add to Seq_001 etc
wholetax <- do.call(paste, c(as.data.frame(tax_table(ps_new))
                  [c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 
				  sep = "__"))
				  
# generate otu table with taxonomy-headers for export
otu_export <- as.data.frame(otu_table(ps_new))
tmp <- names(otu_export)

for(i in 1:length(tmp)){
names(tmp)[i] = paste(wholetax[i], tmp[i], sep = "__")
}

names(otu_export) <- names(tmp)