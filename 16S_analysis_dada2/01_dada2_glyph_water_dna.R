# adjust dada2 for own samples

# check meta data csv

# define working directory to story RData image
setwd("/data/projects/glyphosate/reads/dada2_processed/water_dna")
load(file = "dada2_water_dna.RData")

save.image(file = "dada2_water_dna.RData")
### the outcommented steps are only required the first time

# first installing these packages manually
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("knitr", "BiocStyle"))

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

set.seed(100)

miseq_path <- file.path("/data/projects/glyphosate/reads/reads_16S_cutadapt", 
						"water_dna")
filt_path <- file.path("/data/projects/glyphosate/reads/dada2_processed", 
						"water_dna")


# here starts the processing
fns <- sort(list.files(miseq_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]


ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }

if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

# Trim and Filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(280, 225),
                  maxN = 0, maxEE = c(2.5, 5), truncQ = 2, rm.phix = FALSE,
                  trimLeft = c(10, 0), minLen = c(270, 220), 
				  compress = TRUE, multithread = TRUE)

# reads after trimming
head(out)

# as percentage 
out[,2]/out[,1]*100

# now check the quality after trimming
ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(filtFs[i]) + ggtitle("Fwd filt")) }
for(i in ii) { print(plotQualityProfile(filtRs[i]) + ggtitle("Rev filt")) }


# dereplicate reads, keep abundance and quality information
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
# split character adjusted
sam.names <- sapply(strsplit(basename(filtFs), "_S"), `[`, 1) # split character adjusted
names(derepFs) <- sam.names
names(derepRs) <- sam.names

# Infer sequence variants with training subset
ddF <- dada(derepFs[1:60], err = NULL, 
			selfConsist = TRUE, multithread = TRUE) 
			# Convergence after 6 rounds
		
ddR <- dada(derepRs[1:60], err = NULL, 
			selfConsist = TRUE, multithread = TRUE) 
			# Convergence after 6 rounds
			
plotErrors(ddF, nominalQ = TRUE)
plotErrors(ddR, nominalQ = TRUE)	

dadaFs <- dada(derepFs, err = ddF[[1]]$err_out, pool = TRUE, multithread = TRUE) 
dadaRs <- dada(derepRs, err = ddR[[1]]$err_out, pool = TRUE, multithread = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# Construct sequence table and remove chimeras
seqtab.all <- makeSequenceTable(mergers)
# number of samples and uniq seqs
dim(seqtab.all)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.all)))
hist(nchar(getSequences(seqtab.all)), 
	main = "Distribution of sequence lengths", 
	breaks = (seq(250, 500, by = 10)))

# pick expected and abundant sequence lengths and check again
seqtab2 <- seqtab.all[,nchar(colnames(seqtab.all)) %in% seq(392, 448)]

table(nchar(getSequences(seqtab2)))
hist(nchar(getSequences(seqtab2)), 
	main = "Distribution of sequence lengths", 
	breaks = (seq(390, 450, by = 5)))
	
# remove chimera based on denoised sequences
# checks if both parts of a contig represents either part of more abundant seqs
seqtab2.nochim <- removeBimeraDenovo(seqtab2, method = "consensus", 
					multithread = TRUE, verbose = TRUE)

# number of samples and uniq seqs
dim(seqtab2.nochim)

# ratio of non-chimeric to total seqs
# chimeras can make up many but low abundant seqs
sum(seqtab2.nochim)/sum(seqtab2)


# not tested
table(nchar(getSequences(seqtab2.nochim)))
hist(nchar(getSequences(seqtab2.nochim)), 
	main = "Distribution of sequence lengths", 
	breaks = (seq(390, 450, by = 5)))
	
	
	
	
# include deseq2 here?
# we could probably combine the different sequencing runs here before assigning taxonomy

# Assign taxonomy (I used silva instead of rdp classifier)
# train set and species assignment can also be created from the original silva release file
# path <- "data/db/Silva.nr_v132"
# dada2:::makeTaxonomyFasta_Silva(file.path(path, "silva.nr_v132.align"), file.path(path, "silva.nr_v132.tax"), "data/db/silva_nr_v132_train_set.fa.gz")
# dada2:::makeSpeciesFasta_Silva("data/db/SILVA_132_SSURef_tax_silva.fasta.gz", "~/tax/silva_species_assignment_v132.fa.gz")

ref_fasta <- "/data/db/silva_nr_v132_train_set.fa.gz"
taxtab <- assignTaxonomy(seqtab, refFasta = ref_fasta)
colnames(taxtab) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
seqs <- getSequences(seqtab)
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
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)



# don't forget the controls in meta data? 
# check MIMARKS example

# Combine data into a phyloseq object (csv download link only displayed in V1 of the protocol)					  
# mimarks_path <- file.path(tempdir(), "MIMARKS.csv")
# download.file("https://www.dropbox.com/s/wpxnvsui4mjj8z3/MIMARKS_Data_combined.csv?raw=1",
#                 destfile = mimarks_path)
# samdf <- read.csv(mimarks_path, header = TRUE)
samdf <- read.csv("/data/projects/glyphosate/analysis/metadata/meta_dna_water.csv",
					header = TRUE, row.names = 1)
all(rownames(seqtab) %in% rownames(samdf)) # TRUE
				
ps_water_dna <- phyloseq(tax_table(taxtab), sample_data(samdf),
                 otu_table(seqtab, taxa_are_rows = FALSE),phy_tree(fitGTR$tree))

# To begin, create a table of read counts for each e.g. Genus present in the dataset.

# Show available ranks in the dataset
rank_names(ps_water_dna)

## [1] "Kingdom" "Phylum" "Class" "Order" "Family" "Genus"

# Create table, number of features for each phyla
table(tax_table(ps_water_dna)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps_water_dna),
                 MARGIN = ifelse(taxa_are_rows(ps_water_dna), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps_water_dna),
                      tax_table(ps_water_dna))
					  
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# how to get an OTU table? 
# https://rdrr.io/bioc/phyloseq/man/otu_table-methods.html
# Important informationh here! https://joey711.github.io/phyloseq/import-data.html

otu_table(ps_water_dna) # without taxonomy, but actual sequence as header
table_water_dna <- otu_table(ps_water_dna)
write.csv(table_water_dna, file = "/data/projects/glyphosate/analysis/dada2/otu_table_water_dna.csv")
