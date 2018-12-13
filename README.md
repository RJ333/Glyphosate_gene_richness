#### set up `dada2` V1.8

`dada2` is a R package provided by the bioconductor platform. Apart from primer removal, it can perform all steps required for the generation of an OTU table with taxonomic classification. `dada2` does not call them OTUs, however, they use machine learning approaches to infer sequencing errors and sequence variants, resulting in one specific sequence for each ASV (amplicon sequence variant). Traditional OTUs use several parameters to cluster similar sequences and then define one representative sequence. Therefore, `dada2` claims to find the exact sequence per organism.

The only step prior to using `dada2` is to remove primers from the sequences. We will use `cutadapt` for this step. 

`dada2` infers sequencing errors, and all sequencing runs and sample habitats have their own set of error rates. Therefore, I put all data which should be corrected separately (e.g. from different runs or different habitats) into separate directories. My 3 data sets are 
* water column, DNA
* water column, cDNA
* biofilm, DNA and cDNA

The following directories were set up:

```bash
# these dirs contain the gzipped raw reads
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/water_dna
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/water_cdna
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/biofilm

# for reads without primers (after cutadapt)
/data/projects/glyphosate/reads/reads_16S_cutadapt/water_dna
...

# for reads filtered and trimmed by dada2
/data/projects/glyphosate/reads/dada2_processed/water_dna
...

# we also need databases for taxonomic classification in a dada2 appropriate format
# dada2 also provides functions to generate them out of the standard formatted Silva databases
mkdir -p /data/db
wget -O /data/db/silva_nr_v132_train_set.fa.gz\
  'https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1'
wget -O /data/db/silva_species_assignment_v132.fa.gz\
  'https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1'

# print download time stamp for databases
cd /data/db
for db in silva*fa.gz; do stat $db | grep "File" -A 5 | sed -e 1b -e '$!d'; done
```
Regarding the data bases: It should be possible to reduce the sequences within the taxonomic data base to the range of your primer object. This is e.g. done in `mothur` by applying the `pcr.seqs()` command. This would reduce the alignment effort. But the taxonomic annotation did not seem to take so long, so currently I will not follow this idea. But it could be important for someone with very large and diverse data sets. 

We are then ready to install `dada2` and `cutadapt` in a `conda` environment
```bash
# install Bioconductor and dada2 in conda
conda env create -f conda_dada2.conf
conda activate dada2
R
```
Within R, you probably have to install more packages
```r
# setup bioconductor software including dada2 and phyloseq
source("http://bioconductor.org/biocLite.R")

# these two are dependencies for the later steps, let's install them first
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
```

Some steps may take long, e.g. the error rate training. You can save a workspace to continue working from an advanced step. You should specify a working dir containing the primerless reads before saving or loading
* set working directory `setwd("/data/projects/glyphosate/reads/dada2_processed/water_dna")`
* save the workspace with `save.image(file = "dada2_water_dna.RData")`
* and load it with `load(file = "dada2_water_dna.RData")`



#### Removal of primers

The primers are removed using the provided cutadapt script. See the manual: https://cutadapt.readthedocs.io/en/stable/guide.html

To run the script, you need to provide input and output dir, direction of reads (R1 or R2), number of cores and the primer sequence to remove.

#### The `dada2` tutorial

My workflows are based upon this published one https://f1000research.com/articles/5-1492/v2 and additionally https://benjjneb.github.io/dada2/tutorial.html. If you want to follow this tutorial, check out `tutorial_dada2.r`

#### Running `dada2` on own data

TODO: the parallel runs somehow need to be optimized, only very small parts of code differ between the data sets.

I performed parallel runs of `dada2`, as each sequencing run data needs interactive and individual caretaking: the plotting of the sequence quality will show you where to trim your reads. The outcome of the scripts `02_dada2...r` are `.RData`-files, which contain a sequence table (Seqtab). At the moment, all workspaces need to be loaded to merge the Seqtabs from the parallel runs into one object (`03_dada2_glyph_water_mergetest.r`). Those sequences will then be assigned taxonomically and turned into an `phyloseq`-object.

#### File export such as sequences or count tables

If you want to export files from the `dada2`-object such as the corrected ASV-sequences or the ASV-table with the counts, you can use the following code, here shown for a sequence table called `seqtab2.nochimera`. 

```r
# export the counts and the sequences, which also have their sequence as fasta header
# this allows us to identify same OTUs when we combine more data sets, as the sequence is the same 
asv_seqs <- colnames(seqtab2.nochimera)
asv_headers <- paste0(">", asv_seqs)

# combine fasta header and seq before saving ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file = "ASV.fasta")

# adjust and save according count table:
asv_tab <- t(seqtab2.nochimera)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, file = "ASVs_counts.tsv", 
					 sep = "\t", 
					 quote = FALSE)
```

#### Exporting ASV table for phyloseq analysis

TODO: explain `phyloseq`, the import of data, the repoducability and comparability of results processed by different workflows such as `mothur`, `dada2`, `qiime` or `Silva`

#### Version info:

```r
> sessionInfo()
R version 3.5.1 (2018-07-02)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS: /home/centos/miniconda3/envs/dada2/lib/R/lib/libRblas.so
LAPACK: /home/centos/miniconda3/envs/dada2/lib/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets
[8] methods   base

other attached packages:
 [1] ShortRead_1.40.0            GenomicAlignments_1.18.0
 [3] SummarizedExperiment_1.12.0 DelayedArray_0.8.0
 [5] matrixStats_0.54.0          Biobase_2.42.0
 [7] Rsamtools_1.34.0            GenomicRanges_1.34.0
 [9] GenomeInfoDb_1.18.1         Biostrings_2.50.1
[11] XVector_0.22.0              IRanges_2.16.0
[13] S4Vectors_0.20.1            BiocParallel_1.16.1
[15] BiocGenerics_0.28.0         dada2_1.10.0
[17] Rcpp_1.0.0

loaded via a namespace (and not attached):
 [1] RColorBrewer_1.1-2     pillar_1.3.0           compiler_3.5.1
 [4] plyr_1.8.4             bindr_0.1.1            bitops_1.0-6
 [7] tools_3.5.1            zlibbioc_1.28.0        lattice_0.20-38
[10] tibble_1.4.2           gtable_0.2.0           pkgconfig_2.0.2
[13] rlang_0.3.0.1          Matrix_1.2-15          bindrcpp_0.2.2
[16] GenomeInfoDbData_1.2.0 stringr_1.3.1          hwriter_1.3.2
[19] dplyr_0.7.8            grid_3.5.1             tidyselect_0.2.5
[22] data.table_1.11.8      glue_1.3.0             R6_2.3.0
[25] latticeExtra_0.6-28    reshape2_1.4.3         ggplot2_3.1.0
[28] purrr_0.2.5            magrittr_1.5           scales_1.0.0
[31] assertthat_0.2.0       colorspace_1.3-2       stringi_1.2.4
[34] RCurl_1.95-4.11        lazyeval_0.2.1         RcppParallel_4.4.1
[37] munsell_0.5.0          crayon_1.3.4
```