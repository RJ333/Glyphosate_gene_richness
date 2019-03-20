#!/usr/bin/env Rscript

# read mothur files and prepare data subsets
# constaxonomy contains the taxonomy of the OTUs
our_shared_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.shared"
our_cons_taxonomy_file <- "stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.0.02.cons.taxonomy"

# sample data is in the metafile
metafile <- read.delim("metafile.tsv",
						row.names = 1,
						header = TRUE,
						na.strings = "")
metafile <- sample_data(metafile)

# read table for cell count plot
cell_counts_glyph <- read.csv("cell_counts_glyph.csv", sep = ";")

# read OTU representative sequences
OTU_seqs <- readDNAStringSet(file = "OTU_rephyloseq_object_fasta_002.fasta",
							  format = "fasta",
							  nrec = -1L,
							  skip = 0L,
							  seek.first.rec = FALSE,
							  use.names = TRUE)

# import mothur output into phyloseq
mothur_phyloseq_object<- import_mothur(mothur_list_file = NULL,
						   mothur_group_file = NULL,
						   mothur_tree_file = NULL,
						   cutoff = NULL,
						   mothur_shared_file = our_shared_file,
						   mothur_constaxonomy_file = our_cons_taxonomy_file,
						   parseFunction = parse_taxonomy_default)

# now all files are imported, we can adjust them to our needs
# add further taxonomy columns "OTU" and "wholetax" and adjust column names
wholetax <- do.call(paste, c(as.data.frame(tax_table(mothur_phyloseq_object))
                  [c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6")],
				  sep = "_"))
tax_table(mothur_phyloseq_object) <- cbind(tax_table(mothur_phyloseq_object),
							  rownames(tax_table(mothur_phyloseq_object)),
							  wholetax)
colnames(tax_table(mothur_phyloseq_object)) <- c("kingdom",
									"phylum",
									"class",
									"order",
									"family",
									"genus",
									"otu_id",
									"wholetax")
                                    
# add meta data and OTU representative seqs to phyloseq_objectobject
mothur_full <- merge_phyloseq(mothur_phyloseq_object, metafile, refseq(OTU_seqs))
