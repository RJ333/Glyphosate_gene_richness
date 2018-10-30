
# getting data and folders set up

mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/all_runs
mkdir -p /data/projects/glyphosate/reads/processed

# create soft links to gzipped read files
cp /data/Rene/miseq/biofilm_dna_cdna_600/fastq/*.fastq /data/projects/glyphosate/reads/raw_reads_16S/all_runs
cp /data/Rene/miseq/water2_dna_600/fastq/raw/*.fastq /data/projects/glyphosate/reads/raw_reads_16S/all_runs
cp /data/Rene/miseq/water3_cdna_600/fastq/*.fastq /data/projects/glyphosate/reads/raw_reads_16S/all_runs
cp /data/Rene/miseq/lars_janine/fastq/R*/*.fastq /data/projects/glyphosate/reads/raw_reads_16S/all_runs


ln -s /data/Rene/miseq/biofilm_dna_cdna_600/fastq/*.fastq /data/projects/glyphosate/reads/raw_reads_16S

# test link for each origin folder

zcat Dpositiv_S47_L001_R1_001.fastq | head

# remove sample Mo7 (Christin B.)
rm Mo7*.gz

# check for current silva database
ls /data/db/silva

# unzip silva database ??

##### this is for conda for the old version 1.40.05, which throws an 
##### error using count.seqs() if raw reads were not unzipped

# for quick working mothur on chandler-1
conda create -n mothur_1395 mothur=1.39.5

##### this is for the prerelease 1.41.0, which works so far
PATH=/data/Rene/mothur_prerelease_1_41_0/mothur/mothur:$PATH 
# did not work, maybe $PATH added to wrong start file?
# bashrc is not read at startup
# so we call manually
/data/Rene/mothur_prerelease_1_41_0/mothur/mothur