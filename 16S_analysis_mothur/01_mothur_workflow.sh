# getting data and folders set up

mkdir -p /data/projects/glyphosate/reads/raw_reads

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

# remove undetermined reads

rm *determined*.gz