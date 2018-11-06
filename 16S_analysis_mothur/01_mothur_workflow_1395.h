# this script uses mothur interactively! It is not designed to run on its own

# the first steps are a little complicated: 1) the mothur default path are not 
# correctly pasted for the output folder, therefore you can't have the 
# stability.files in the output folder

# 2) generating the stability.files in the first place sets a new input directory
# the stability.files later needs to be adjusted manually, as the reads 
# contain the whole path

# when the stability.files is adjusted and back in the folder, 
# where the original reads are, you
# can continue 

# this still requires a lot of playing around, not final yet

# mothur does not like line breaks, even with \

# You need to set the input and output in the same command. 
# Mothur clears the the values that are not set when set.dir is run.
# mothur > set.dir(input=../data/MiSeq_SOP, output=../analyses/)
# set input dir for raw reads and output dir for all further files
set.dir(input = /data/projects/glyphosate/reads/mothur_processed, output = /data/projects/glyphosate/reads/mothur_processed)
set.logfile(name = glyphosate_start_1395.log)

# make stability file
# in this version, the stability file needs to be adjusted manually
# from local shell
scp -i ~/.ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/stability.files /mnt/d/mothur1395/
# adjusting e.g. in excel
scp -i ~/.ssh/denbi.key /mnt/d/mothur1395/stability2.files.txt centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/

# back on VM
mv /data/projects/glyphosate/reads/mothur_processed/stability2.files.txt /data/projects/glyphosate/reads/mothur_processed/stability.files
# Group names should not include :, -, or / characters
make.file(inputdir = /data/projects/glyphosate/reads/raw_reads_16S/, type = gz)

# make.file resets the directories, so reset them too!
set.dir(input = /data/projects/glyphosate/reads/raw_reads_16S/, output = /data/projects/glyphosate/reads/mothur_processed)
make.contigs(file = stability.files, processors = 28)

# now we don't use the folder of the raw reads anymore
set.dir(input = /data/projects/glyphosate/reads/mothur_processed, output = /data/projects/glyphosate/reads/mothur_processed)
summary.seqs(fasta = stability.trim.contigs.fasta, processors = 28)
count.groups(group = stability.contigs.groups)

# file oligo.txt copied from Lars Möller, 
# sequence of 341f-805r as would be described in paper-style 
# --> (805r not reversed-complementary)
trim.seqs(fasta = stability.trim.contigs.fasta, oligos = oligo.txt, processors = 28)	

# check quality of primer-trimmed reads
summary.seqs(fasta = stability.trim.contigs.trim.fasta, processors = 28)
# and the removed reads
summary.seqs(fasta = stability.trim.contigs.scrap.fasta, processors = 28)

# we have to create a new group file now
list.seqs(fasta = stability.trim.contigs.trim.fasta)
get.seqs(accnos = stability.trim.contigs.trim.accnos, group = stability.contigs.groups)

# removes(!) sequences which do not apply to the defined parameters from dataset 
# (can be found in the accnos file)
screen.seqs(fasta = stability.trim.contigs.trim.fasta, group = stability.contigs.pick.groups, maxambig = 0, maxhomop = 8, maxlength = 450, minlength = 390, processors = 28) 

# optional, summarises number of sequences per group (eg sample) 
count.groups(group = stability.contigs.pick.good.groups)

# looks for unique sequences and merges identical sequences 
# under one reference unique sequence
unique.seqs(fasta = stability.trim.contigs.trim.good.fasta)

# creates a table with the unique sequences in each group 
# and the number of times each unique sequence 
# is found in total and per group (sample)-> important for abundance data!
count.seqs(name = stability.trim.contigs.trim.good.names, group = stability.contigs.pick.good.groups)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.fasta, count = stability.trim.contigs.trim.good.count_table, processors = 28)

# I copied the pcr.seqs and silva files from Lars Möller, 
# as he used the same primers.
pcr.seqs(fasta = /data/db/silva.seed_v132.align, start = 6387, end = 23442, keepdots = F, processors = 28)
	
# trims the curated sequences of the reference database to a certain section 
# (i.e. start/end) -->  this needs to be done only once
# you can use the silva.nr_v132.pcr.align file for all your next analyses, 
# as long as you are working with the same primer set (-->same fragment!) 
# and there is no newer version of the reference database 
# (check https://mothur.org/wiki/Silva_reference_files)

# we use seed_v132 to align, this is the reduced database
align.seqs(fasta = stability.trim.contigs.trim.good.unique.fasta, reference = silva.seed_v132.pcr.align, processors = 28)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.align, count = stability.trim.contigs.trim.good.count_table, processors = 28)

# again removing sequences that dont match 
# -> if Start or End do not fit or NBases is off
screen.seqs(fasta = stability.trim.contigs.trim.good.unique.align, count = stability.trim.contigs.trim.good.count_table, summary = stability.trim.contigs.trim.good.unique.summary, start = 41, end = 17053, processors = 28)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.align, count = stability.trim.contigs.trim.good.good.count_table)

# removes gap characters within the sequences (e.g. '-') 
# to condense file and save disk space
filter.seqs(fasta = stability.trim.contigs.trim.good.unique.good.align, vertical = T)

# because after all the trimming and aligning, 
# perhaps more redundant sequences can be found, 
# thats why we use unique.seqs again
unique.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.fasta, count = stability.trim.contigs.trim.good.good.count_table)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.count_table)

# pre-clustering merges sequences that are 2 nt different (diffs=) 
# from each other (1 nt per 100 bp) -> 
# this allows to take mutations into account (Katis V4 reads were 250 bp reads)

# I will stick to 2 nt diff although I have 300 bp reads
pre.cluster(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.count_table, diffs = 2)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.count_table)

# now we want to detect chimeric sequences in our data set; 
# dereplicate=true (t) means that if one sequence gets flagged as chimeric 
# in one group, it is NOT automatically flagged as chimeric in other groups #
# (because it might be an abundant sequence in another group) 

# I copied the vsearch bin from mothur executable into the conda environment
chimera.vsearch(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.count_table, dereplicate = t, processors = 28)

# optional: count seqs per sample 
count.groups(count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)

#removing the chimeric sequences from our dataset (fasta file)
remove.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, accnos = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)

# classifying the sequences (before OTU picking!) based on the 
# complete SILVA reference set (non redundant?) 
# cutoff = 85 (85%) means that everything gets classified 
# as 'non-classified' when there not at least 85% probability (similarity?) 
# to the next phylo level; try probs=F to remove bootstrap values from 
# taxonomy
classify.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference = /data/db/silva.nr_v132.align, taxonomy = /data/db/silva.nr_v132.tax, cutoff = 85)

# removing phyla/lineages that we dont want in our dataset 
# (e.g. Eukaryota, mitochondria etc because the primers are not 
# supposed to match those)
remove.lineage(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon = Chloroplast-Mitochondria-unknown-Eukaryota)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, processors = 28)

##### OTU picking based on 98% similarity 
# (memory intensive, performed on phy-2 with 1.39.5)

# store private key at bio48/49, adjust to chmod 600, 
# use to copy data from denbi-cloud to bio-48 

# on bio-48:
cd /data/projects/glyphosate/analysis_16S/mothur_1_39_5
scp -i /data/Rene/ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta .
scp -i /data/Rene/ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table .
scp -i /data/Rene/ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy .

# then run mothur, make sure to get the right version and use screen/tmux
ssh phy-2
screen
/dss6/bio49/bin/mothur/mothur
set.dir(input = /dss6/bio49/projects/glyphosate/analysis_16S/mothur_1_39_5, output = /dss6/bio49/projects/glyphosate/analysis_16S/mothur_1_39_5/)
cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.02, processors = 100)
/dss6/bio49/projects/glyphosate/reads/processed/

#################### this is where I stopped

make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label = 0.02)
# (label not important here?)
# a shared file that contains all the information on which OTU was found how many times in which group (sample)