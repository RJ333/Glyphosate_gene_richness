### bullet list with each step of mothur pipeline used

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

# aligns sequences with the reference database
# below is the large database
# align.seqs(fasta = stability.trim.contigs.trim.good.unique.fasta, reference = silva.nr_v132.pcr.align)
# the seed_v132 is the reduced database
align.seqs(fasta = stability.trim.contigs.trim.good.unique.fasta, reference = silva.seed_v132.pcr.align, processors = 24)

It took 7816 secs to align 1552629 sequences.

summary.seqs(fasta = stability.trim.contigs.trim.good.unique.align, count = stability.trim.contigs.trim.good.count_table)

Using 24 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       37      5       0       2       1
2.5%-tile:      41      17053   402     0       4       614332
25%-tile:       41      17053   402     0       4       6143312
Median:         41      17053   402     0       5       12286623
75%-tile:       41      17053   415     0       5       18429934
97.5%-tile:     41      17053   427     0       6       23958914
Maximum:        16198   17053   447     0       8       24573245
Mean:   45      17052   408     0       4
# of unique seqs:       1552629
total # of seqs:        24573245

It took 296 secs to summarize 24573245 sequences.

Output File Names:
 /data/projects/glyphosate/reads/processed/stability.trim.contigs.trim.good.unique.summary

# again removing sequences that dont match 
# -> starting and ending at the wrong position (look at minimun/maximum in summary) 
screen.seqs(fasta = stability.trim.contigs.trim.good.unique.align, count = stability.trim.contigs.trim.good.count_table, summary = stability.trim.contigs.trim.good.unique.summary, start = 41, end = 17053)

summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.align, count = stability.trim.contigs.trim.good.good.count_table)

Using 24 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       17053   390     0       3       1
2.5%-tile:      41      17053   402     0       4       612231
25%-tile:       41      17053   402     0       4       6122302
Median:         41      17053   402     0       5       12244603
75%-tile:       41      17053   421     0       5       18366904
97.5%-tile:     41      17053   427     0       6       23876974
Maximum:        41      17053   447     0       8       24489204
Mean:   		40      17053   408     0       4
# of unique seqs:       1526850
total # of seqs:        24489204

It took 232 secs to summarize 24489204 sequences.

Output File Names:
 /data/projects/glyphosate/reads/processed/stability.trim.contigs.trim.good.unique.good.summary

# removes gap characters within the sequences (e.g. '-') 
# to condense file and save disk space
filter.seqs(fasta = stability.trim.contigs.trim.good.unique.good.align, vertical = T)

# because after all the trimming and aligning, 
# perhaps more redundant sequences can be found, 
# thats why we use unique.seqs again
unique.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.fasta, count = stability.trim.contigs.trim.good.good.count_table)
summary.seqs(fasta=stability.trim.contigs.trim.good.unique.good.filter.unique.fasta, count=stability.trim.contigs.trim.good.unique.good.filter.count_table)

Using 24 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       1027    390     0       3       1
2.5%-tile:      22      1027    402     0       4       612231
25%-tile:       22      1027    402     0       4       6122302
Median:         22      1027    402     0       5       12244603
75%-tile:       22      1027    421     0       5       18366904
97.5%-tile:     22      1027    427     0       6       23876974
Maximum:        22      1027    447     0       8       24489204
Mean:   		21      1027    408     0       4
# of unique seqs:       1526774
total # of seqs:        24489204

It took 179 secs to summarize 24489204 sequences.

Output File Names:
 /data/projects/glyphosate/reads/processed/stability.trim.contigs.trim.good.unique.good.filter.unique.summary
 
# pre-clustering merges sequences that are 2 nt different (diffs=) 
# from each other (1 nt per 100 bp) -> 
# this allows to take mutations into account (Katis 250 bp reads)

# I will stick to 2 nt diff although I have 300 bp reads
pre.cluster(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.count_table, diffs = 2)

summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.count_table)

Using 24 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       1027    390     0       3       1
2.5%-tile:      22      1027    402     0       4       612231
25%-tile:       22      1027    402     0       4       6122302
Median:         22      1027    402     0       5       12244603
75%-tile:       22      1027    422     0       5       18366904
97.5%-tile:     22      1027    427     0       6       23876974
Maximum:        22      1027    447     0       8       24489204
Mean:   21      1027    408     0       4
# of unique seqs:       638888
total # of seqs:        24489204

It took 76 secs to summarize 24489204 sequences.

Output File Names:
 /data/projects/glyphosate/reads/processed/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.summary

#now we want to detect chimeric sequences in our data set; 
# dereplicate=true (t) means that if one sequence gets flagged as chimeric 
# in one group, it is NOT automatically flagged as chimeric in other groups #
# (because it might be an abundant sequence in another group) 
chimera.vsearch(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.count_table, dereplicate = t)

# optional count seqs per sample 
count.groups(count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)

#removing the chimeric sequences from our dataset (fasta file)
remove.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, accnos = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)

summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)

Using 50 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       1027    390     0       3       1
2.5%-tile:      22      1027    402     0       4       592324
25%-tile:       22      1027    402     0       4       5923237
Median:         22      1027    402     0       5       11846474
75%-tile:       22      1027    407     0       5       17769710
97.5%-tile:     22      1027    427     0       6       23100623
Maximum:        22      1027    440     0       8       23692946
Mean:   21      1027    408     0       4
# of unique seqs:       523074
total # of seqs:        23692946

It took 61 secs to summarize 23692946 sequences.

Output File Names:
 /data/projects/glyphosate/reads/processed/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.su
mmary

# classifying the sequences (before OTU picking!) based on the 
# complete SILVA reference set (non redundant?) 
# cutoff = 85 (85%) means that everything gets classified 
# as 'non-classified' when there not at least 85% probability (similarity?) 
# to the next phylo level; try probs=F to remove bootstrap values from 
# taxonomy
classify.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference = silva.nr_v132.align, taxonomy = silva.nr_v132.tax, cutoff = 85)
It took 4548 secs to classify 523074 sequences.

# removing phyla/lineages that we dont want in our dataset 
# (e.g. Eukaryota, mitochondria etc because the primers are not 
# supposed to match those)

# remove.lineage(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon = Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
# keeping Archea
remove.lineage(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon = Chloroplast-Mitochondria-unknown-Eukaryota)

summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, processors = 50)

# summary with Archea removed

Using 30 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       1027    390     0       3       1
2.5%-tile:      22      1027    402     0       4       592324
25%-tile:       22      1027    402     0       4       5923234
Median:         22      1027    402     0       5       11846467
75%-tile:       22      1027    407     0       5       17769700
97.5%-tile:     22      1027    427     0       6       23100610
Maximum:        22      1027    440     0       8       23692933
Mean:   21      1027    408     0       4
# of unique seqs:       523068
total # of seqs:        23692933

It took 40 secs to summarize 23692933 sequences.

# summary including Archea

Using 50 processors.

                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       1027    390     0       3       1
2.5%-tile:      22      1027    402     0       4       592324
25%-tile:       22      1027    402     0       4       5923234
Median:         22      1027    402     0       5       11846467
75%-tile:       22      1027    407     0       5       17769700
97.5%-tile:     22      1027    427     0       6       23100610
Maximum:        22      1027    440     0       8       23692933
Mean:   21      1027    408     0       4
# of unique seqs:       523068
total # of seqs:        23692933

It took 54 secs to summarize 23692933 sequences.

# OTU picking based on 98% similarity: (very long step)
# crashed on bio-49 with all cores, no message except for "killed"
# using bio-48 with 22 cores
cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.02, processors = 22)

# Using 22 processors.
# Splitting the file...
# Running command: 
# dist.seqs(fasta=/data/projects/glyphosate/reads/processed/stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta.0.temp, 
# processors=22, cutoff=0.02, outputdir=/data/projects/glyphosate/reads/processed/)

Sequence  Time    Num_Dists_Below_Cutoff
44400     122     247322

#################### this is where I left off


make.shared(list=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label=0.03)
(label hier vll nicht wichtig)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared		# a shared file that contains all the information on which OTU was found how many times in which group (sample)
-> rename to: stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list_all.shared

split.abund(fasta=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, list=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff=3)
Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rare.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.rare.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.rare.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table


summary.seqs(fasta=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.fasta, count=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table)
-> output
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        19      741     244     0       3       1
2.5%-tile:      23      741     253     0       3       349349
25%-tile:       23      741     253     0       4       3493489
Median:         23      741     253     0       4       6986978
75%-tile:       23      741     253     0       5       10480466
97.5%-tile:     23      741     253     0       6       13624606
Maximum:        23      741     259     0       8       13973954
Mean:   23      741     253.012 0       4.41639
# of unique seqs:       147018
total # of seqs:        13973954

Output File Names: 
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.summary

#now we create a new shared file based on the new count table, where the 'singletons' have been removed
make.shared(list=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list, count=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table, label=0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared			#this is the actual OTU table with reads per OTU per sample

#now we can classify our OTUs
classify.otu(list=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list, count=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table, taxonomy=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, label=0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.taxonomy		#taxonomy file with taxonomical info for each OTU!
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.tax.summary	#different version of taxonomy file with taxonomical info for each OTU!

count.groups(shared=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.count.summary		#table with reads per sample 

get.relabund(shared=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.relabund -> totalgroup		#OTU table with relative abundances

get.relabund(shared=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared, scale=totalotu)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.relabund -> totalotu

get.relabund(shared=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared, scale=averageotu)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.relabund -> averageotu


get.groups(count=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table, list=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list, taxonomy=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.taxonomy, groups=B1-B2-B3-B4-DEPC-leer1-leer2-Negativ-NegativPlatte4-PB1-X247-X248-X249-X274-X275-X276-X301-X302-X303)
-> Output File names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.pick.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.pick.taxonomy

split.abund(list=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.list, count=stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.pick.count_table, cutoff=10)
