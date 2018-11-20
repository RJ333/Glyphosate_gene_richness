# this script installs dada2 within a conda environment and downloads the
# required databases
# the sequence reads from different runs/environments 
# should be stored in separate directories (or have a grep-able pattern)
# raw reads in 
/data/projects/glyphosate/reads/raw_reads_16S/water_dna
/data/projects/glyphosate/reads/raw_reads_16S/water_cdna
/data/projects/glyphosate/reads/raw_reads_16S/biofilm
# cutadapt reads in
/data/projects/glyphosate/reads/reads_16S_cutadapt/water_dna
/data/projects/glyphosate/reads/reads_16S_cutadapt/water_cdna
/data/projects/glyphosate/reads/reads_16S_cutadapt/biofilm
# dada2 filtered and trimmed reads in
/data/projects/glyphosate/reads/dada2_processed/water_dna
/data/projects/glyphosate/reads/dada2_processed/water_cdna
/data/projects/glyphosate/reads/dada2_processed/biofilm
# download the silva training set and the silva species assignment
# with the links above, those file can also be created from original silva db
mkdir -p /data/db
wget -O /data/db/silva_nr_v132_train_set.fa.gz\
  'https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1'
wget -O /data/db/silva_species_assignment_v132.fa.gz\
  'https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1'

# we can add a meta data file to phyloseq  
mkdir -p /data/projects/glyphosate/analysis/metadata
mkdir -p /data/projects/glyphosate/analysis/dada2
scp -i ~/.ssh/denbi.key /mnt/d/denbi/meta_dna_water.csv centos@193.196.20.111:/data/projects/glyphosate/analysis/metadata/

# first set up conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# install Bioconductor and dada2 in conda, 
conda create -n dada2 bioconductor-shortread=1.38.0 r-base=3.5.1 bioconductor-dada2=1.8 cutadapt=1.18
conda activate dada2
R

