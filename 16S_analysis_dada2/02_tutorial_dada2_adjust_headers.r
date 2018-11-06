# test the adjustment of the OTU table
# originally, the OTU table contains the actual amplicon sequence as header
# I would like to have a numbered header such "OTU001" or "Seq001" etc
# prefixed by the full taxonomic path

# define working directory to story RData image
setwd("/data/projects/dada2_tutorial")
load(file = "dada2_tutorial.RData")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)


# adjust phyloseq object name to be not overwritten
ps_new <- ps

# this changes the header from the actual sequence to Seq_001, Seq_002 etc
taxa_names(ps_new)
n_seqs <- seq(ntaxa(ps_new))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(ps_new) <- paste("Seq", formatC(n_seqs, 
											width = len_n_seqs, 
											flag = "0"), sep = "_")
taxa_names(ps_new)

# generate taxonomy from tax table as df
wholetax <- do.call(paste, c(as.data.frame(tax_table(ps_new))
                  [c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")], 
				  sep = "__"))
				  
# generate otu table
otu_export <- as.data.frame(otu_table(ps_new))
tmp <- names(otu_export)

for(i in 1:length(tmp)){
names(tmp)[i] = paste(wholetax[i], tmp[i], sep = "__")
}

names(otu_export) <- names(tmp)