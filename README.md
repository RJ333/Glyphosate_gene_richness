#### set up dada2

dada2 is a R package provided by the bioconductor platform. Apart from primer removal, it can perform all steps required for the generation of an OTU table with taxonomic classification. dada2 does not call them OTUs, however, they use machine learning approaches to infer sequencing errors and sequence variants, resulting in one specific sequence for each ASV (amplicon sequence variant). Traditional OTUs use several parameters to cluster similar sequences and then define one representative sequence.

The only step prior to using dada2 is to remove primers from the sequences. We will use `cutadapt` for this step. 

dada2 infers sequencing errors, which usually differ between sequencing runs and habitats. Therefore, I put all data which should be corrected separately (e.g. from different runs or different habitats) into separate directories. My 3 data sets are 
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
for db in silva*.tgz; do stat $db | grep "File" -A 5 | sed -e 1b -e '$!d'; done
```
We are then ready to install dada2 and cutadapt in a conda environment
```bash
# install Bioconductor and dada2 in conda, 
conda create -n dada2 bioconductor-shortread=1.38.0 r-base=3.5.1 bioconductor-dada2=1.8 cutadapt=1.18
conda activate dada2
R
```
within R, you probably still have to install more packages

```r
# setup bioconductor software including dada2 and phyloseq
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
```
#### removal of primers

The primers are removed using the provided cutadapt script. See the manual: https://cutadapt.readthedocs.io/en/stable/guide.html

To run the script, you need to provide input and output dir, direction of reads (R1 or R2), number of cores and the primer sequence to remove 

#### dada2 tutorial

My workflows are based upon this published one https://f1000research.com/articles/5-1492/v2 and additionally https://benjjneb.github.io/dada2/tutorial.html. If you want to follow this tutorial, check out `tutorial_dada2.r`

#### dada2 on own data

TODO: the parallel runs somehow need to be optimized, only very small parts of code differ between the data sets.

I performed parallel runs of dada2, as each sequencing run data needs interactive and individual caretaking: the plotting of the sequence quality will show you where to trim your reads. The outcome of the scripts `02_dada2...r` are `.RData`-files, which contain a sequence table (Seqtab). At the moment, all workspaces need to be loaded to merge the Seqtabs from the parallel runs into one object (`03_dada2_glyph_water_mergetest.r`). Those sequences will then be assigned taxonomically and turned into an `phyloseq`-object.

#### phyloseq analysis

TODO: explain phyloseq, the import of data, the repoducability and comparability of results processed by different workflows such as mothur, dada2, qiime or Silva