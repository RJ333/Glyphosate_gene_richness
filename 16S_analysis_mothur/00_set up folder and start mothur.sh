# getting data and folders set up

mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/biofilm
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/water_dna
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/water_cdna
mkdir /data/projects/glyphosate/reads/mothur_processed

wget -O /data/db/Silva.seed_v132.tgz 'https://mothur.org/w/images/7/71/Silva.seed_v132.tgz'
wget -O /data/db/Silva.nr_v132.tgz 'https://mothur.org/w/images/3/32/Silva.nr_v132.tgz'
cd /data/db 
tar -xvzf /data/db/Silva.seed_v132.tgz
tar -xvzf /data/db/Silva.nr_v132.tgz
 
# x for extract
# v for verbose
# z for gnuzip
# f for file, should come at last just before file name
ln -s /data/db/silva.seed_v132.align /data/projects/glyphosate/reads/mothur_processed

# create soft links to gzipped read files in one target folder
# the other folders are kept as they are the source for dada2 separate processing
ln -s /data/projects/glyphosate/reads/raw_reads_16S/*/*.gz /data/projects/glyphosate/reads/raw_reads_16S/

# check for current silva database
ls /data/db/silva

# for quick working mothur on chandler-1
conda create -n mothur_1395 mothur=1.39.5

# if vseach is missing in conda environment and you have another mothur version installed
cp /data/programs/mothur/vsearch /home/centos/miniconda3/envs/mothur_1395/bin/