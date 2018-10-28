# download silva training set from 
# https://www.dropbox.com/sh/mfcivbudmc21cqt/AAB1l-AUM5uKvjrR33ct-cTXa?dl=0

# https://rdrr.io/bioc/dada2/man/makeSpeciesFasta_Silva.html
# https://zenodo.org/record/1172783 prepared silva 132 training and species assignment


# download the silva training set and the silva species assignment
# with the links above, those can also be created

mkdir -p /data/db

wget -O /data/db/silva_nr_v132_train_set.fa.gz\
  'https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1'
wget -O /data/db/silva_species_assignment_v132.fa.gz\
  'https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1'
 
mkdir -p /data/projects/glyphosate/analysis/metadata
mkdir .p /data/projects/glyphosate/analysis/dada2
scp -i ~/.ssh/denbi.key /mnt/d/denbi/meta_dna_water.csv centos@193.196.20.111:/data/projects/glyphosate/analysis/metadata/

scp -i ~/.ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/analysis/dada2/otu_table_water_dna.csv /mnt/d/denbi/
 
# install Bioconductor and dada2 in conda dada2
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n dada2 bioconductor-shortread=1.38.0 r-base=3.5.1 bioconductor-dada2=1.8
conda activate dada2
conda install cutadapt=1.18
sudo yum install pigz  # allows parallel gzip


cd /data/projects/glyphosate/reads/raw_reads_16S/water_dna

vim cutadapt_dna_water.sh
#!/bin/bash
# this script assumes that f and r reads are still separated
# and primer set 341f-805r (V3-V4)

input=$1 # argument passing has not been tested yet
output=$2

for xy in $1/*R*.gz
do
  cutadapt -j 20 -g CCTACGGGNGGCWGCAG -g GACTACHVGGGTATCTAATCC -o $2/cut_${xy} ${xy} 
done

R