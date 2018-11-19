# create directories with 16S reads for mothur analysis
# one directory for each MiSeq run, so we can use the same structure for dada2
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/biofilm
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/water_dna
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/water_cdna
mkdir /data/projects/glyphosate/reads/mothur_processed
# download, extract and link the desired taxonomy database (here: Silva 132)
wget -O /data/db/Silva.seed_v132.tgz 'https://mothur.org/w/images/7/71/Silva.seed_v132.tgz'
wget -O /data/db/Silva.nr_v132.tgz 'https://mothur.org/w/images/3/32/Silva.nr_v132.tgz'
cd /data/db 
tar -xvzf /data/db/Silva.seed_v132.tgz
tar -xvzf /data/db/Silva.nr_v132.tgz
# link 
ln -s /data/db/silva.seed_v132.align /data/projects/glyphosate/reads/mothur_processed

# create soft links to gzipped read files in one mothur target folder
# the physical reads are kept in seperate folders
# as they are the source for dada2 run-separated processing
ln -s /data/projects/glyphosate/reads/raw_reads_16S/*/*.gz /data/projects/glyphosate/reads/raw_reads_16S/

# create conda environment for mothur on chandler-1
conda create -n mothur_1395 mothur=1.39.5

# if vsearch is missing in conda environment and you have another mothur version installed
cp /data/programs/mothur/vsearch /home/centos/miniconda3/envs/mothur_1395/bin/