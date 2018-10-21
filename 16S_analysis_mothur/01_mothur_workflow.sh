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

cat Dpositiv_S47_L001_R1_001.fastq | head

ln -s /data/Rene/miseq/water2_dna_600/fastq/raw/*.fastq /data/projects/glyphosate/reads/raw_reads_16S

cat 65_S47_L001_R1_001.fastq | head

ln -s /data/Rene/miseq/water3_cdna_600/fastq/*.fastq /data/projects/glyphosate/reads/raw_reads_16S

cat c65_S46_L001_R1_001.fastq | head

ln -s /data/Rene/miseq/lars_janine/fastq/R*/*.fastq /data/projects/glyphosate/reads/raw_reads_16S

cat R12_S96_L001_R1_001.fastq | head

# remove sample Mo7 (Christin B.)

rm Mo7*.gz

# current silva database
mv /data/Rene/SILVA_132_SSURef_NR99_13_12_17_opt.arb.gz /data/db/silva

# unzip silva database ??

##### this is for the old version 1.40.05, which throws an error using count.seqs() if raw reads
##### were not unzipped
# to start mothur in conda, add the following lines to bashrc on bio-49 and save
# vim ~/.bashrc
# . /data/tools/miniconda3/etc/profile.d/conda.sh
# conda activate
# source ~/.bashrc
# start mothur environment
# conda activate mothur
#mothur

# we need to make sure to call the newest version, prerelease 1.41.0
PATH=/data/Rene/mothur_prerelease_1_41_0/mothur/mothur:$PATH 
# did not work, bashrc is also not used instartup, maybe profile?
# so we call manually
/data/Rene/mothur_prerelease_1_41_0/mothur/mothur


# set output dir
set.dir(input=/data/projects/glyphosate/reads/processed)
set.dir(output=/data/projects/glyphosate/reads/processed)

# make stability file
# Group names should not include :, -, or / characters
# c-negativ und c-positiv wurden nicht bearbeitet
make.file(inputdir=/data/projects/glyphosate/reads/raw_reads_16S/all_runs, type = gz)

# stability file cleaves names after first "_", therefore c_negative and c_positive are both "c" now
make.contigs(file = stability.files, processors = 22)
	
# 1.40.5, gz-files: It took 56538 secs to process 40476806 sequences.
# 1.41.0, fastq-files
summary.seqs(fasta = stability.trim.contigs.fasta, processors = 24)
# saved output to overview_contigs.txt and stability.trim.contigs.summary

count.groups(group = stability.contigs.groups)

# file oligo.txt copied from Lars Möller, 341f-805r described in paper-style (not reversed-complementary)
trim.seqs(fasta = stability.trim.contigs.fasta, oligos = oligo.txt, processors = 24)	

# It took 1128 secs to trim 41418035 sequences.

# check quality of primer-trimmed reads
summary.seqs(fasta = stability.trim.contigs.trim.fasta, processors = 24)
summary.seqs(fasta = stability.trim.contigs.scrap.fasta, processors = 24)

# we have to create a new group file now
list.seqs(fasta = stability.trim.contigs.trim.fasta)
get.seqs(accnos = stability.trim.contigs.trim.accnos, group = stability.contigs.groups)

# mothur does not like line breaks, even with \
# removes(!) sequences which do not apply to the defined parameters from dataset 
# (can be found in the accnos file)
screen.seqs(fasta = stability.trim.contigs.trim.fasta, group = stability.contigs.pick.groups, maxambig = 0, maxhomop = 8, maxlength = 450, minlength = 390, processors = 24) 

# optional, summarises number of sequences per group (eg sample) 
count.groups(group = stability.contigs.pick.good.groups)

# looks for unique sequences and merges identical sequences under one reference unique sequence
unique.seqs(fasta = stability.trim.contigs.trim.good.fasta)

# creates a table with the unique sequences in each group 
# and the number of times each unique sequence 
# is found in total and per group (sample)-> important for abundance data!
count.seqs(name = stability.trim.contigs.trim.good.names, group = stability.contigs.pick.good.groups)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.fasta, count = stability.trim.contigs.trim.good.count_table)

Using 24 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       390     390     0       3       1
2.5%-tile:      1       402     402     0       4       614332
25%-tile:       1       402     402     0       4       6143312
Median:         1       402     402     0       5       12286623
75%-tile:       1       421     421     0       5       18429934
97.5%-tile:     1       427     427     0       6       23958914
Maximum:        1       450     450     0       8       24573245
Mean:   1       408     408     0       4
# of unique seqs:       1552629
total # of seqs:        24573245

It took 186 secs to summarize 24573245 sequences.

Output File Names:
 /data/projects/glyphosate/reads/processed/stability.trim.contigs.trim.good.unique.summary
 
# I copied the pcr.seqs and silva files from lmoeller, as he used the same primers.
# should still be checked again

#pcr.seqs(fasta=silva.nr_v132.align, start=?????, end=????, keepdots=F, processors=20)
#-> Output File Names: 
#silva.nr_v132.pcr.align	
#trims the curated sequences of the reference database to a certain section (i.e. start/end) --> this needs to be done only once, you can use the silva.nr_v132.pcr.align file for all your next analyses, as long as you are working with the same primer set (-->same fragment!) and there is no newer version of the reference database (check https://mothur.org/wiki/Silva_reference_files)

# align.seqs(fasta = stability.trim.contigs.trim.good.unique.fasta, reference = silva.nr_v132.pcr.align)
align.seqs(fasta = stability.trim.contigs.trim.good.unique.fasta, reference = silva.seed_v132.pcr.align, processors = 24)



