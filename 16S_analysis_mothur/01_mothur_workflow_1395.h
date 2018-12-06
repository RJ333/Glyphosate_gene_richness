# this script uses mothur interactively!

# use screen, tmux or byobu!

# set input dir for raw reads and output dir for to make stability file
set.dir(input = /data/projects/glyphosate/reads/raw_reads_16S/, output = /data/projects/glyphosate/reads/mothur_processed)
set.logfile(name = glyphosate_1395.log, append = T)

# Group names should not include :, -, or / characters
# gives error ( inputdir not exist or not writable) and writes to output dir
make.file(inputdir = /data/projects/glyphosate/reads/raw_reads_16S, type = gz)

# the generated stability.files is not correct, use the awk command from README.md (in another terminal) to adjust the file

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
	
# pcr.seqs() trims the curated sequences of the reference database to a 
# certain section (i.e. start/end) -->  this needs to be done only once
# you can use the silva.nr_v132.pcr.align file for all your next analyses, 
# as long as you are working with the same primer set (--> same fragment!) 
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

# performend on phy-2
# check the README.md if you want to transfer the files to phy-2, which has more RAM
cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03, processors = 10)

# this is the count table with absolute values (check whether the file names include "0.02" or similar)
make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

# now we can classify our OTUs
classify.otu(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, label = 0.02)

count.groups(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.shared)

# generate a fasta file containing a representative sequence (here based on most abundant) per OTU, this will be imported to phyloseq for tree calculation
get.oturep(method = abundance, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta)