# this script uses mothur interactively!

# use screen, tmux or byobu!
# Group names (sample names) must not include :, -, or / characters

# set input dir for raw reads and output dir for to make stability file
set.dir(input = /data/projects/glyphosate/reads/raw_reads_16S/, output = /data/projects/glyphosate/reads/mothur_processed)
set.logfile(name = glyphosate_1395.log, append = T)


# this command returns error and then writes to output dir
make.file(inputdir = /data/projects/glyphosate/reads/raw_reads_16S, type = gz)

# the generated stability.files is not correct, use the awk command (see below) from README.md (in another terminal) to adjust the file
// # move file into required dir
// cd /data/projects/glyphosate/reads/mothur_processed
// mv stability.files /data/projects/glyphosate/reads/raw_reads_16S/stabil.temp
// cd /data/projects/glyphosate/reads/raw_reads_16S/
// # adjust group name and path
// awk 'BEGIN{OFS="\t"}; {sub(".*/", "", $2); sub(".*/", "", $3); print $1, $2, $3}' stabil.temp |\
// awk 'BEGIN{FS = "_| "; OFS = "\t"}; {print $1, $0}' | awk 'BEGIN{OFS = "\t"}; {print $2, $4, $5}' > stability.files

# make contigs from forward and reverse, using the stability.files
make.contigs(file = stability.files, processors = 27)

# now we don't use the folder of the raw reads anymore
set.dir(input = /data/projects/glyphosate/reads/mothur_processed, output = /data/projects/glyphosate/reads/mothur_processed)
summary.seqs(fasta = stability.trim.contigs.fasta, processors = 28)
count.groups(group = stability.contigs.groups)

# oligo.txt contains the primer sequences to be removed 
# see outcommented example below for formatting
// cat oligo.txt
// forward CCTACGGGNGGCWGCAG
// reverse GACTACHVGGGTATCTAATCC 
trim.seqs(fasta = stability.trim.contigs.fasta, oligos = oligo.txt, processors = 28)	

# check stats of remaining reads after primer removal
summary.seqs(fasta = stability.trim.contigs.trim.fasta, processors = 28)
# and the stats of the removed reads
summary.seqs(fasta = stability.trim.contigs.scrap.fasta, processors = 28)

# update the files after primer removal
list.seqs(fasta = stability.trim.contigs.trim.fasta)
get.seqs(accnos = stability.trim.contigs.trim.accnos, group = stability.contigs.groups)

# perform quality trimming on the sequences 
screen.seqs(fasta = stability.trim.contigs.trim.fasta, group = stability.contigs.pick.groups, maxambig = 0, maxhomop = 8, maxlength = 450, minlength = 390, processors = 28) 

# optional, summarises number of sequences per group 
count.groups(group = stability.contigs.pick.good.groups)

# looks for unique sequences and dereplicates identical sequences 
# under one reference unique sequence
unique.seqs(fasta = stability.trim.contigs.trim.good.fasta)

# creates a table with the unique sequences in each group and how often each unique sequence 
# is found in total and per group (sample) -> important for abundance data!
count.seqs(name = stability.trim.contigs.trim.good.names, group = stability.contigs.pick.good.groups)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.fasta, count = stability.trim.contigs.trim.good.count_table, processors = 28)

# pcr.seqs() trims the curated sequences of the reference database to a certain section (i.e. start/end)
# I copied the pcr.seqs and silva files from Lars Möller (we both used 341f-805r primer set).
# If you don't know the start/end parameters of your primers, use SINA to identify the regions
pcr.seqs(fasta = /data/db/silva.seed_v132.align, start = 6387, end = 23442, keepdots = F, processors = 28)
	
# we use seed_v132 (seed is reduced SILVA database) to align our sequences 
align.seqs(fasta = stability.trim.contigs.trim.good.unique.fasta, reference = silva.seed_v132.pcr.align, processors = 28)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.align, count = stability.trim.contigs.trim.good.count_table, processors = 28)

# again removing sequences that don't match start or end values (based on previous summary.seqs()-output)
screen.seqs(fasta = stability.trim.contigs.trim.good.unique.align, count = stability.trim.contigs.trim.good.count_table, summary = stability.trim.contigs.trim.good.unique.summary, start = 41, end = 17053, processors = 28)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.align, count = stability.trim.contigs.trim.good.good.count_table)

# the aligned seqs contain gap characters (e.g. '-'), remove them to save disk space
filter.seqs(fasta = stability.trim.contigs.trim.good.unique.good.align, vertical = T)

# after trimming and aligning more redundant sequences might be found, 
# so we dereplicate again
unique.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.fasta, count = stability.trim.contigs.trim.good.good.count_table)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.count_table)

# pre-clustering merges sequences that are 2 nt different (diffs=) 
# from each other (1 nt per 100 bp) -> 
# this allows to take mutations into account (Katis V4 reads were 250 bp reads)

# I will stick to diffs=2 (which is conservative for 300 bp reads) to keep more possible OTUs
pre.cluster(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.count_table, diffs = 2)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.count_table)

# detect chimeric sequences in our data set; 
# dereplicate=true (t) means that if one sequence gets flagged as chimeric 
# in one group, it is NOT automatically flagged as chimeric in other groups
# (because it might be an abundant sequence in another group) 
chimera.vsearch(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.count_table, dereplicate = t, processors = 28)

# optional: count seqs per sample 
count.groups(count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)

# removing the chimeric sequences from our dataset (fasta file)
remove.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, accnos = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)

# classifying the sequences (before OTU picking!) based on the 
# complete SILVA reference set (non redundant?) 
# cutoff = 85 (85%) means that everything gets classified 
# as 'non-classified' below 85% probability (similarity?) 
# to the next phylo level; try probs=F to remove bootstrap values from 
# taxonomy
classify.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference = /data/db/silva.nr_v132.align, taxonomy = /data/db/silva.nr_v132.tax, cutoff = 85)

# removing phyla/lineages that we dont want in our dataset (e.g. Eukaryota, mitochondria) 
# or should not appear because the primers are not supposed to match those
remove.lineage(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon = Chloroplast-Mitochondria-unknown-Eukaryota)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, processors = 28)

##### OTU picking based on 98% similarity 

# performend on phy-2
# check the README.md if you want to transfer the files to phy-2, which has more RAM
cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03, processors = 10)

# this is the count table with absolute values (check whether the file names include the label "0.02" or similar)
make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)

# now we can classify our OTUs
classify.otu(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, label = 0.02)

count.groups(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.shared)

# generate a fasta file containing a representative sequence (here based on most abundant) per OTU, this will be imported to phyloseq for tree calculation
get.oturep(method = abundance, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta)