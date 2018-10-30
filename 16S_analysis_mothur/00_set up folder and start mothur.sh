# getting data and folders set up

mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/biofilm
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/water_dna
mkdir -p /data/projects/glyphosate/reads/raw_reads_16S/water_cdna
mkdir /data/projects/glyphosate/reads/mothur_processed


# create soft links to gzipped read files in one target folder
# the other folders are kept as they are the source for dada2 separate processing
ln -s /data/projects/glyphosate/reads/raw_reads_16S/*/*.gz /data/projects/glyphosate/reads/raw_reads_16S/

# check for current silva database
ls /data/db/silva

# for quick working mothur on chandler-1
conda create -n mothur_1395 mothur=1.39.5