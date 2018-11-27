# this script uses mothur interactively!

# the first steps are a little complicated:
# 1) the group names in the stability files are wrong (contain full path), 
# you have to adjust them using awk. "group" means sample or library 

# 2) mothur looks in different folders, but there is a pasting error with "/" so it won't work
# until it is in the inputdir. but for some reason, the stability.files will be put into the output dir

# so after creating stability.files, we move and rename it to the input folder (where the reads are)
# we modify it using awk (outside of mothur in second shell) and then can use make.contigs


# mothur does not like line breaks, even with \

# You need to set the input and output in the same command. 
# Mothur clears the values that are not set when calling set.dir


# use screen or tmux or byobu!

# set input dir for raw reads and output dir for to make stability file
set.dir(input = /data/projects/glyphosate/reads/raw_reads_16S/, output = /data/projects/glyphosate/reads/mothur_processed)
set.logfile(name = glyphosate_water_1395.log, append = T)

# Group names should not include :, -, or / characters
# gives error ( inputdir not exist or not writable) and writes to output dir
make.file(inputdir = /data/projects/glyphosate/reads/raw_reads_16S, type = gz)

# in second shell!
# use awk to remove the path and exchange the first column 
# with the second column until the first underscore
cd /data/projects/glyphosate/reads/mothur_processed
mv stability.files /data/projects/glyphosate/reads/raw_reads_16S/stabil.temp
cd  /data/projects/glyphosate/reads/raw_reads_16S/
awk 'BEGIN{OFS="\t"}; {sub(".*/", "", $2); sub(".*/", "", $3); print $1, $2, $3}' stabil.temp |\
awk 'BEGIN{FS = "_| "; OFS = "\t"}; {print $1, $0}' | awk 'BEGIN{OFS = "\t"}; {print $2, $4, $5}' > stability.files

# make contigs from forward and reverse, using the stability.files
make.contigs(file = stability.files, processors = 27)

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

# optional, summarises number of sequences per group 
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

# as it is only water samples and the files are not so big, I'll give it a first try on my VM
# (memory intensive, performed on phy-2 with 1.39.5)
cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03, processors = 10)

# store private key at bio48/49, adjust to chmod 600, 
# use to copy data from denbi-cloud to bio-48 

# on bio-48:
#cd /data/projects/glyphosate/analysis_16S/mothur_1_39_5
#scp -i /data/Rene/ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta .
#scp -i /data/Rene/ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table .
#scp -i /data/Rene/ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy .

# then run mothur, make sure to get the right version and use screen/tmux
#ssh phy-2
#screen
#/dss6/bio49/bin/mothur/mothur
#set.dir(input = /dss6/bio49/projects/glyphosate/analysis_16S/mothur_1_39_5, output = /dss6/bio49/projects/glyphosate/analysis_16S/mothur_1_39_5/)
# killed again with 100 cores due to memory exceedment
# but 60 cores and cutoff 0.02 worked
#cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.02, processors = 60)

# I'm also trying to repeat cluster.split() with cutoff 0.03. This crashed several times, now I'm using only 10 cores on phy-2
# mv output to subfolder cut_off_002, the rerun with cutoff 003
#ssh phy-2
#screen
#/dss6/bio49/bin/mothur/mothur
#set.dir(input = /dss6/bio49/projects/glyphosate/analysis_16S/mothur_1_39_5, output = /dss6/bio49/projects/glyphosate/analysis_16S/mothur_1_39_5/cut_off_003)
#cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03, processors = 10)

# you can take up the cluster.split()-steps from the sens.spec()-step. But I assume that the memory demands 
# are already determined before ("splitting file" in the beginning of cluster.split()) and it will crash again, regardless of the number of cores set with sens.spec()

# copy files back to the cloud, except for dist-file, which is huge and not needed
# on bio-48:
cd /data/projects/glyphosate/analysis_16S/mothur_1_39_5
# scp -i /data/Rene/ssh/denbi.key /data/projects/glyphosate/analysis_16S/mothur_1_39_5/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.dist centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/
scp -i /data/Rene/ssh/denbi.key /data/projects/glyphosate/analysis_16S/mothur_1_39_5/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/
scp -i /data/Rene/ssh/denbi.key /data/projects/glyphosate/analysis_16S/mothur_1_39_5/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.sensspec centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/

# this is the count table with absolute values
make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

# rename output file outside mothur, as this is with singletons
mv stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list_all.shared

# remove singletons
split.abund(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff = 1)

# check how much we loose
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.02.abund.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.02.abund.count_table, processors = 24)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.02.rare.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.02.rare.count_table, processors = 24)

# now we create a new shared file based on the new count table, where the 'singletons' have been removed
make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.02.abund.count_table, label = 0.02)

# now we can classify our OTUs
classify.otu(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.02.abund.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, label = 0.02)

count.groups(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.shared)

# OTU table with relative abundances -> total group
get.relabund(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.shared)
# OTU table with relative abundances -> total OTU
get.relabund(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.shared, scale = totalotu)
# OTU table with relative abundances -> average OTU
get.relabund(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.shared, scale = averageotu)

# generate a fasta file containing a representative sequence (here based on most abundant) per OTU, this will be imported to phyloseq for tree calculation
get.oturep(method = abundance, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.02.abund.count_table, fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.02.abund.fasta)

# to copy files
scp -i /drives/d/ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/reads/mothur_processed/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.abund.relabund /mnt/d/data/mothur_sample/