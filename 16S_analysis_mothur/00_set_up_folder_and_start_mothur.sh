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
# link (make sure databases are present, as link will be created anyway)
ln -s /data/db/silva.seed_v132.align /data/projects/glyphosate/reads/mothur_processed

# create soft links to gzipped read files in one mothur target folder
# the physical reads are kept in seperate folders
# as they are the source for dada2 run-separated processing (only water samples, no controls)
ln -s /data/projects/glyphosate/reads/raw_reads_16S/water*/*.gz /data/projects/glyphosate/reads/raw_reads_16S/
cd /data/projects/glyphosate/reads/raw_reads_16S/
rm *pos*.gz
rm *neg*.gz
# create conda environment for mothur, make sure channels are configured
conda create -n mothur_1395 mothur=1.39.5

# if vsearch is missing in conda environment and you have another mothur version installed
cp /data/programs/mothur/vsearch /home/centos/miniconda3/envs/mothur_1395/bin/
# or from bio-48
scp -i /data/Rene/ssh/denbi.key /data/bin/mothur/vsearch centos@193.196.20.103:/home/centos/miniconda3/envs/mothur_1395/bin/
# start mothur
conda activate mothur_1395
mothur

# meta data for phyloseq import from local pc
scp -i /drives/d/ssh/denbi.key /mnt/d/all_samples_with_meta.tsv centos@193.196.20.103:/data/projects/glyphosate/analysis/metadata/all_samples_with_meta_cond2.tsv

# run script to adjust otu reps to proper fasta sequences for import to phyloseq
cd /data/projects/glyphosate/reads/mothur_processed

python /home/centos/scripts/otu_rep_to_fasta.py \
  -i stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.fasta \
  -o OTU_reps_fasta_water_003.fasta