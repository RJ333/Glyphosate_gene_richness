# getting data and folders set up

mkdir -p /data/projects/glyphosate/reads/raw_reads_16S
mkdir -p /data/projects/glyphosate/reads/processed
# create soft links to gzipped read files

ln -s /data/Rene/miseq/biofilm_dna_cdna_600/fastq/*.gz /data/projects/glyphosate/reads/raw_reads

# test link for each origin folder

zcat Dpositiv_S47_L001_R1_001.fastq.gz | head

ln -s /data/Rene/miseq/water2_dna_600/fastq/raw/*.gz /data/projects/glyphosate/reads/raw_reads

zcat 65_S47_L001_R1_001.fastq.gz | head

ln -s /data/Rene/miseq/water3_cdna_600/fastq/*.gz /data/projects/glyphosate/reads/raw_reads

zcat c65_S46_L001_R1_001.fastq.gz | head

ln -s /data/Rene/miseq/lars_janine/fastq/R*/*.gz /data/projects/glyphosate/reads/raw_reads

zcat R12_S96_L001_R1_001.fastq.gz | head

# remove undetermined reads and sample Mo7

rm *determined*.gz
rm Mo7*.gz

for file in c-*; do mv "$file" `echo $file | tr '-' '_'` ; done

# current silva database
mv /data/Rene/SILVA_132_SSURef_NR99_13_12_17_opt.arb.gz /data/db/silva

# unzip silva database ??


# to start mothur in conda, add the following lines to bashrc on bio-49 and save
vi ~/.bashrc

# . /data/tools/miniconda3/etc/profile.d/conda.sh
# conda activate

source ~/.bashrc

# start mothur environment
conda activate mothur

mothur

# set output dir
set.dir(output=/data/projects/glyphosate/reads/processed)

# make stability file
# Group names should not include :, -, or / characters
make.file(inputdir=/data/projects/glyphosate/reads/raw_reads_16S, type = gz)

make.contigs(file=stability.files, processors=45)	

