# adjust dada2 for own samples

# check meta data csv

# define working directory to story RData image
setwd("/data/projects/glyphosate/reads/dada2_processed/water_cdna")
load(file = "dada2_water_cdna.RData")

### the outcommented steps are only required the first time

# first installing these packages manually
source("http://bioconductor.org/biocLite.R")
biocLite(c("knitr", "BiocStyle"))

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

miseq_path <- file.path("/data/projects/glyphosate/reads/reads_16S_cutadapt", "water_cdna")
filt_path <- file.path("/data/projects/glyphosate/reads/dada2_processed", "water_cdna")


# here starts the processing
fns <- sort(list.files(miseq_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

# Trim and Filter
ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }

# forward reads 10 to 280
# reverse reads 10 to 240

if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))			# be careful, the original script 
filtRs <- file.path(filt_path, basename(fnRs))			# contains filtFs and filt"s"Fs!!
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
		      c(filtFs[[i]], filtRs[[i]]),
                      trimLeft = 10, truncLen = c(280, 240),
                      maxN = 0, maxEE = 2, truncQ = 2,
                      compress = TRUE)
}

# Trim and Filter
ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(filtFs[i]) + ggtitle("Fwd filt")) }
for(i in ii) { print(plotQualityProfile(filtRs[i]) + ggtitle("Rev filt")) }


# Infer sequence variants
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtFs), "_S"), `[`, 1) # split character adjusted
names(derepFs) <- sam.names
names(derepRs) <- sam.names

# adjusted to 50 samples for training
ddF <- dada(derepFs[1:50], err = NULL, 
			selfConsist = TRUE) # Convergence after  6  rounds. about 2 hours?

			
ddR <- dada(derepRs[1:50], err = NULL, 
			selfConsist = TRUE) # Convergence after 

# plotErrors(ddF)
# plotErrors(ddR)	

dadaFs <- dada(derepFs, err = ddF[[1]]$err_out, pool = TRUE, multithread = TRUE) # for multiple cores
dadaRs <- dada(derepRs, err = ddR[[1]]$err_out, pool = TRUE, multithread = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# Construct sequence table and remove chimeras
seqtab.all <- makeSequenceTable(mergers)
seqtab <- removeBimeraDenovo(seqtab.all)

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
alignment <- AlignSeqs(cdnaStringSet(seqs), anchor = NA)

# phangorn package: we first construct a neighbor-joining tree,
phang.align <- phyDat(as(alignment, "matrix"), type = "cdna")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data = phang.align)

# then fit a GTR+G+I maximum likelihood tree using the neighbor-joining tree as a starting point.
fitGTR <- update(fit, k = 4, inv = 0.2)
# the following step takes longer
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                      rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)

save.image(file = "dada2_water_cdna.RData")

# don't forget the controls in meta data? 
# check MIMARKS example

# Combine data into a phyloseq object (csv download link only displayed in V1 of the protocol)					  
# mimarks_path <- file.path(tempdir(), "MIMARKS.csv")
# download.file("https://www.dropbox.com/s/wpxnvsui4mjj8z3/MIMARKS_Data_combined.csv?raw=1",
#                 destfile = mimarks_path)
# samdf <- read.csv(mimarks_path, header = TRUE)
samdf <- read.csv("/data/projects/glyphosate/analysis/metadata/meta_cdna_water.csv",
					header = TRUE, row.names = 1)
all(rownames(seqtab) %in% rownames(samdf)) # TRUE
				
ps_water_cdna <- phyloseq(tax_table(taxtab), sample_data(samdf),
                 otu_table(seqtab, taxa_are_rows = FALSE),phy_tree(fitGTR$tree))

# To begin, create a table of read counts for each e.g. Genus present in the dataset.

# Show available ranks in the dataset
rank_names(ps_water_cdna)

## [1] "Kingdom" "Phylum" "Class" "Order" "Family" "Genus"

# Create table, number of features for each phyla
table(tax_table(ps_water_cdna)[, "Phylum"], exclude = NULL)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps_water_cdna),
                 MARGIN = ifelse(taxa_are_rows(ps_water_cdna), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps_water_cdna),
                      tax_table(ps_water_cdna))
					  
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# how to get an OTU table? 
# https://rdrr.io/bioc/phyloseq/man/otu_table-methods.html
# Important informationh here! https://joey711.github.io/phyloseq/import-data.html

otu_table(ps_water_cdna) # without taxonomy, but actual sequence as header
table_water_cdna <- otu_table(ps_water_cdna)
write.csv(table_water_cdna, file = "/data/projects/glyphosate/analysis/dada2/otu_table_water_cdna.csv")
