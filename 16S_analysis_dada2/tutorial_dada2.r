# https://benjjneb.github.io/dada2/bigdata.html
# https://benjjneb.github.io/dada2/tutorial.html

# Samples have been demultiplexed, i.e. split into individual per-sample fastq files.
# Non-biological nucleotides have been removed, e.g. primers, adapters, linkers, etc.
# If paired-end sequencing data, the forward and reverse fastq files contain reads in matched order

# following tutorial on https://f1000research.com/articles/5-1492/v2

# define working directory to story RData image
setwd("/data/projects/dada2_tutorial")
load(file = "dada2_tutorial.RData")


### the outcommented steps are only required the first time

# first installing these packages manually
# source("http://bioconductor.org/biocLite.R")
# biocLite(c("knitr", "BiocStyle"))

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

# set.seed(100)

# miseq_path <- file.path("/data/projects/dada2_tutorial", "MiSeq_SOP")
# filt_path <- file.path("/data/projects/dada2_tutorial", "filtered")

# download sample data

# if(!file_test("-d", miseq_path)) {
  # dir.create(miseq_path)
  # download.file("http://www.mothur.org/MiSeqDevelopmentData/StabilityNoMetaG.tar",
                 # destfile = file.path(miseq_path, "StabilityNoMetaG.tar"))
  # system(paste0("tar -xvf ", file.path(miseq_path, "StabilityNoMetaG.tar"),
                 # " -C ", miseq_path, "/"))
# }

# here starts the processing
fns <- sort(list.files(miseq_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]

# Trim and Filter
ii <- sample(length(fnFs), 3)
for(i in ii) { print(plotQualityProfile(fnFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(fnRs[i]) + ggtitle("Rev")) }

# be careful, the original script contains typo "filtsFs"
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, basename(fnFs))			         
filtRs <- file.path(filt_path, basename(fnRs))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[[i]], fnRs[[i]]),
		      c(filtFs[[i]], filtRs[[i]]),
                      trimLeft = 10, truncLen = c(245, 160),
                      maxN = 0, maxEE = 2, truncQ = 2,
                      compress = TRUE)
}

# Infer sequence variants
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names

ddF <- dada(derepFs[1:40], err = NULL, 
			selfConsist = TRUE)
ddR <- dada(derepRs[1:40], err = NULL, 
			selfConsist = TRUE)

plotErrors(ddF) 
plotErrors(ddR)	

dadaFs <- dada(derepFs, err = ddF[[1]]$err_out, pool = TRUE, multithread = TRUE)
dadaRs <- dada(derepRs, err = ddR[[1]]$err_out, pool = TRUE, multithread = TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)

# Construct sequence table and remove chimeras
seqtab.all <- makeSequenceTable(mergers[!grepl("Mock", names(mergers))])
seqtab <- removeBimeraDenovo(seqtab.all)

# Assign taxonomy (I used silva instead of rdp classifier)
ref_fasta <- "/data/projects/dada2_tutorial/db/silva_nr_v123_train_set.fa.gz"
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

# Combine data into a phyloseq object (csv download link only displayed in V1 of the protocol)					  
mimarks_path <- file.path(tempdir(), "MIMARKS.csv")
download.file("https://www.dropbox.com/s/wpxnvsui4mjj8z3/MIMARKS_Data_combined.csv?raw=1",
                destfile = mimarks_path)
samdf <- read.csv(mimarks_path, header = TRUE)
samdf$SampleID <- paste0(gsub("00", "", samdf$host_subject_id), "D", samdf$age-21)
samdf <- samdf[!duplicated(samdf$SampleID),] # Remove dupicate entries for reverse reads
rownames(seqtab) <- gsub("124", "125", rownames(seqtab)) # Fixing an odd discrepancy
all(rownames(seqtab) %in% samdf$SampleID) # TRUE
				
rownames(samdf) <- samdf$SampleID
keep.cols <- c("collection_date", "biome", "target_gene", "target_subfragment",
"host_common_name", "host_subject_id", "age", "sex", "body_product", "tot_mass",
"diet", "family_relationship", "genotype", "SampleID") 
samdf <- samdf[rownames(seqtab), keep.cols]
				
ps <- phyloseq(tax_table(taxtab), sample_data(samdf),
                 otu_table(seqtab, taxa_are_rows = FALSE),phy_tree(fitGTR$tree))

# change sequence header to OTU names, check gitlab
ps_new <- ps
head(otu_table(ps_new)[,1])
taxa_names(ps_new) <- paste0("Seq", seq(ntaxa(ps_new)))
head(otu_table(ps_new)[,1])				 
# if you want to start from this checkpoint, reload object ps

# library("phyloseq")
# library("gridExtra")
# ps = readRDS("data/ps.rds")

rank_names(ps)
table(tax_table(ps)[, "Phylum"], exclude = NULL)
ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
# prevalence in the dataset, which we will define here as 
# the number of samples in which a taxa appears at least once.
prevdf = apply(X = otu_table(ps0),
                 MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2),
                 FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                      TotalAbundance = taxa_sums(ps0),
                      tax_table(ps0))
# Compute the total and average prevalences of the features in each phylum.					  
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
filterPhyla = c("Fusobacteria", "Deinococcus-Thermus")
# Filter entries with unidentified Phylum.
ps1 = subset_taxa(ps0, !Phylum %in% filterPhyla)
ps1

# Prevalence Filtering
# Subset to the remaining phyla
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position = "none")
#  Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(ps0)
prevalenceThreshold

## [1] 18

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps2 = prune_taxa(keepTaxa, ps0)

# Agglomerate taxa
save.image(file = "dada2_tutorial.RData")