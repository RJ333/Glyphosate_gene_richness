
# make sure to the right version of mothur
/data/bin/mothur/mothur  # 1.39.5
set.dir(input = , output = )
	
# creates stability file of unzipped fastq files --> 
# change group names in your stability file manually with excel
make.file(inputdir = POS488_all, type = fastq)

# needs to be reset after make.file
set.dir(input = , output = )

# assembly of forward and reverse
make.contigs(file = stability.files, processors = 50)	
summary.seqs(fasta = stability.trim.contigs.fasta, processors = 20)

# how many bad sequences:
mothur > summary.seqs(fasta = stability.trim.contigs.scrap.fasta, processors = 20) 

# optional, summarises number of sequences per group (eg sample)
count.groups(group = stability.contigs.groups)

# removes primers (or other oligos defined by oligo file)
trim.seqs(fasta = stability.trim.contigs.fasta, oligos = oligo.txt)	
summary.seqs(fasta = stability.trim.contigs.trim.fasta)

# we have to create a new group file now
list.seqs(fasta = stability.trim.contigs.trim.fasta)
get.seqs(accnos = stability.trim.contigs.trim.accnos, group = stability.contigs.groups)

# removes(!) sequences which do not apply to the defined parameters from dataset 
# (can be found in the accnos file)
screen.seqs(fasta = stability.trim.contigs.trim.fasta, group = stability.contigs.pick.groups, maxambig = 0, maxhomop = 8, maxlength = 275, minlength = 250, processors = 40) 
summary.seqs(fasta = stability.trim.contigs.trim.good.fasta)

# optional, summarises number of sequences per group (eg sample) 
count.groups(group = stability.contigs.pick.good.groups)

# looks for unique sequences and merges identical sequences 
# under one reference unique sequence
unique.seqs(fasta = stability.trim.contigs.trim.good.fasta)	

# creates a table with the unique sequences in each group 
# and the number of times each unique sequence 
# is found in total and per group (sample)-> important for abundance data!
count.seqs(name = stability.trim.contigs.trim.good.names, group = stability.contigs.pick.good.groups)	
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.fasta, count = stability.trim.contigs.trim.good.count_table)

# trims the curated sequences of the reference database to 
# a certain section (i.e. start/end) --> this needs to be done only once, 
# you can use the silva.nr_v132.pcr.align file for all your next analyses, 
# as long as you are working with the same primer set (-->same fragment!) and 
# there is no newer version of the reference database 
# (check https://mothur.org/wiki/Silva_reference_files)
pcr.seqs(fasta = silva.nr_v132.align, start = 11894, end = 23445, keepdots = F, processors = 20)

# aligns sequences with the reference database
align.seqs(fasta = stability.trim.contigs.trim.good.unique.fasta, reference = silva.nr_v132.pcr.align)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.align, count = stability.trim.contigs.trim.good.count_table)

# again removing sequences that dont match -> starting and ending at the 
# wrong position (look at minimun/maximum in summary)
screen.seqs(fasta = stability.trim.contigs.trim.good.unique.align, count = stability.trim.contigs.trim.good.count_table, summary = stability.trim.contigs.trim.good.unique.summary, start = 1968, end = 11550)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.align, count = stability.trim.contigs.trim.good.good.count_table)

# removes gap characters within the sequences (e.g. '-') 
# to condense file and save disk space
filter.seqs(fasta = stability.trim.contigs.trim.good.unique.good.align, vertical = T)

# because after all the trimming and aligning, perhaps more redundant sequences 
# can be found, thats why we use unique.seqs again
unique.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.fasta, count = stability.trim.contigs.trim.good.good.count_table)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.count_table)

# pre-clustering merges sequences that are 2 nt different (diffs = ) 
# from each other (1 nt per 100 bp) 
# -> this allows to take mutations into account 
pre.cluster(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.count_table, diffs = 2)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.count_table)

# now we want to detect chimeric sequences in our data set; 
# dereplicate = true (t) means that if one sequence gets flagged as chimeric 
# in one group, it is NOT automatically flagged as chimeric in other groups 
# (because it might be an abundant sequence in another group) 
chimera.vsearch(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.count_table, dereplicate = t)

# optional
count.groups(count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)

# removing the chimeric sequences from our dataset (fasta file)
remove.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.fasta, accnos = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.accnos)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)

# classifying the sequences (before OTU picking!) based on the complete SILVA 
# reference set (non redundant?), cutoff  =  85 (85%) means that everything gets 
# classified as 'non-classified' when there not at least 85% probability 
# (similarity?) to the next phylo level; try probs = F to remove bootstrap 
# values from taxonomy
classify.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference = silva.nr_v132.align, taxonomy = silva.nr_v132.tax, cutoff = 85)

# removing phyla/lineages that we dont want in our dataset (e.g. Eukaryota, 
# mitochondria etc because the primers are not supposed to match those)
remove.lineage(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon = Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, processors = 30)

# OTU picking based on 97% similarity:
cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03)


make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label = 0.03)
(label hier vll nicht wichtig)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared		# a shared file that contains all the information on which OTU was found how many times in which group (sample)
-> rename to: stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list_all.shared

split.abund(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff = 3)
Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rare.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.rare.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.rare.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table


summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table)
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
make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table, label = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared			#this is the actual OTU table with reads per OTU per sample

#now we can classify our OTUs
classify.otu(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, label = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.taxonomy		#taxonomy file with taxonomical info for each OTU!
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.tax.summary	#different version of taxonomy file with taxonomical info for each OTU!

count.groups(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.count.summary		#table with reads per sample 

get.relabund(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.relabund -> totalgroup		#OTU table with relative abundances

get.relabund(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared, scale = totalotu)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.relabund -> totalotu

get.relabund(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared, scale = averageotu)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.relabund -> averageotu


get.groups(count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.taxonomy, groups = B1-B2-B3-B4-DEPC-leer1-leer2-Negativ-NegativPlatte4-PB1-X247-X248-X249-X274-X275-X276-X301-X302-X303)
-> Output File names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.pick.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.pick.taxonomy

split.abund(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.pick.count_table, cutoff = 10)

#####################################################################################################

############# classifying only sequences from OTU00001 again ############################


bin.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.fasta, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table, label = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.fasta -> with count
#fasta-formatted file where sequences are ordered according to the OTU that they belong to (but only the unique sequences!) -> get sequences for OTU00001 manually

bin.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.fasta, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list, label = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.fasta -> without counts

get.oturep(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, column = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.dist, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.fasta
#fasta-formatted sequence file containing only a representative sequence for each OTU

classify.seqs(fasta = OTU_1_uniques.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table, reference = silva.nr_v132.align, taxonomy = silva.nr_v132.tax, cutoff = 85, processors = 15)
-> Output File Names:
OTU_1_uniques.nr_v132.wang.taxonomy
OTU_1_uniques.nr_v132.wang.tax.summary

#using method  =  knn
classify.seqs(fasta = OTU_1_uniques_without_OTU.txt, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference = silva.nr_v132.align, taxonomy = silva.nr_v132.tax, cutoff = 85, method = knn)

((((list.seqs(fasta = OTU_1_uniques.fasta)
-> Output File Names:
OTU_1_uniques.accnos

get.seqs(accnos = OTU_1_uniques.accnos, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.pick.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.pick.count_table)))))

############ getting only Planktomycete sequences #####################################

#getting all sequences of certian lineage after classify.seqs
get.lineage(taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.taxonomy, taxon = Planctomycetes, fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy -> planctomycetes
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta -> planctomycetes
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table -> planctomycetes

list.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctomycetes.fasta)


#getting only represantative sequences for each planctomycetean OTU
cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.wang.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03, processors = 14)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.sensspec

get.oturep(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, column = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.dist, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.fasta


########## classifying with method = knn ###########

#classifying the sequences (before OTU picking!) based on the complete SILVA reference set (non redundant?), cutoff  =  85 (85%) means that everything gets classified as 'non-classified' when there not at least 85% probability (similarity?) to the next phylo level; try probs = F to remove bootstrap values from taxonomy
classify.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, reference = silva.nr_v132.align, taxonomy = silva.nr_v132.tax, cutoff = 85, method = knn)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.knn.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.knn.tax.summary		#taxonomy file with the taxonomical info for each sequence


#removing phyla/lineages that we dont want in our dataset (e.g. Eukaryota, mitochondria etc because the primers are not supposed to match those)
remove.lineage(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.knn.taxonomy, taxon = Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.knn.pick.taxonomy		#new taxonomy file with the lineages removed
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta						#fasta file with specified lineages removed
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table	#new count file with the lineages removed
	

summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
-> output
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        1       741     233     0       3       1
2.5%-tile:      23      741     253     0       3       355993
25%-tile:       23      741     253     0       4       3559925
Median:         23      741     253     0       4       7119850
75%-tile:       23      741     253     0       5       10679775
97.5%-tile:     23      741     254     0       6       13883707
Maximum:        23      742     272     0       8       14239699
Mean:   23      741     253.036 0       4.45044
# of unique seqs:       177463
total # of seqs:        14239699
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.summary

#OTU picking based on 97% similarity:
cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.knn.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.sensspec

make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, label = 0.03)
(label hier vll nicht wichtig)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared		#a shared file that contains all the information on which OTU was found how many times in which group (sample)
-> rename to: stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list_all.shared

split.abund(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table, cutoff = 3)
Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rare.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.rare.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.rare.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table

summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table)
-> output
                Start   End     NBases  Ambigs  Polymer NumSeqs
Minimum:        3       741     244     0       3       1
2.5%-tile:      23      741     253     0       3       354995
25%-tile:       23      741     253     0       4       3549946
Median:         23      741     253     0       4       7099892
75%-tile:       23      741     253     0       5       10649838
97.5%-tile:     23      741     254     0       6       13844789
Maximum:        23      741     272     0       8       14199783
Mean:   23      741     253.035 0       4.45026
# of unique seqs:       150466
total # of seqs:        14199783

Output File Names: 
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.summary
#now we create a new shared file based on the new count table, where the 'singletons' have been removed
make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table, label = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared			#this is the actual OTU table with reads per OTU per sample

#now we can classify our OTUs
classify.otu(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.0.03.abund.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.knn.pick.taxonomy, label = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.taxonomy		#taxonomy file with taxonomical info for each OTU!
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.tax.summary	#different version of taxonomy file with taxonomical info for each OTU!

count.groups(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.count.summary		#table with reads per sample 

get.relabund(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.relabund -> totalgroup		#OTU table with relative abundances

# try to find a core microbiome
get.coremicrobiome(relabund = OTU_rel_komplett_bereinigt_means.relabund, abundance = 1)
-> Output File Names:
OTU_rel_komplett_bereinigt_means.0.03.core.microbiome
OTU_rel_komplett_bereinigt_means.0.03.core.microbiomelist

# trying to find indicator species (all)
indicator(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.shared, design = meta_design_all.design, processors = 15)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.indicator.summary

# trying to find indicator species (subsampled)
indicator(shared = OTU_subsampled_roundend_non_zero.shared, design = metadata_komplett_bereinigt_for_mothur_sample_type_only.design, processors = 15)
-> Output File Names:
OTU_subsampled_roundend_non_zero.indicator.summary

# indicator fot t7 samples
indicator(shared = t7_OTU_subsample_rounded_non_zero.shared, design = design_t7.design)
-> Output File Names:
t7_OTU_subsample_rounded_non_zero.indicator.summary

# indicator fot t7 samples
indicator(shared = t7_OTU_subsample_rounded.shared, design = design_t7.design)
-> Output File Names:
t7_OTU_subsample_rounded.indicator.summary

# indicator fot treatment samples
indicator(shared = treat_OTU_subsample_rounded_non_zero.shared, design = design_treat.design)
-> Output File Name:
treat_OTU_subsample_rounded_non_zero.indicator.summary

#create a file for LefSe input
make.lefse(relabund = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.relabund, constaxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.taxonomy, label = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.lefse

make.lefse(relabund = subsampled_OTUs_for_mothur.txt, constaxonomy = taxonomy_komplett_bereinigt_for_mothur.txt, label = 0.03)
-> Output File Names:
subsampled_OTUs_for_mothur.0.03.lefse

#make a tree based on the subsampled OTU table
list.seqs(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list, taxonomy = taxonomy_komplett_bereinigt_for_mothur.txt)
-> Output File Names:
taxonomy_komplett_bereinigt_for_mothur.accnos

get.seqs(accnos = taxonomy_komplett_bereinigt_for_mothur.accnos, fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.fasta)
# did not work

get.otus(accnos = taxonomy_komplett_bereinigt_for_mothur.accnos, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.list

list.seqs(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.list)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.accnos

get.seqs(accnos = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.accnos, fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.fasta)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.0.03.abund.pick.fasta


#make tree only with rep sequences of subsampled otus
get.oturep(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, column = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.dist, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.fasta

get.seqs(accnos = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.accnos, fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.count_table

clearcut(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.fasta, DNA = T)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre

#using the unifrac.weighted command with above tree and all samples in the count file... making a tree and using the count table from the subsampled otus in R did not work - seems to be a bug of some sort
unifrac.weighted(tree = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.count_table, groups = P009-P055-P069-P082-P109-P171-P218-P238-P251-P261-P262-X004-X005-X006-X007-X008-X010-X011-X012-X013-X014-X015-X016-X017-X018-X019-X020-X021-X022-X023-X024-X025-X026-X027-X028-X029-X030-X031-X032-X039-X040-X041-X042-X043-X044-X045-X046-X047-X048-X049-X050-X051-X052-X053-X054-X056-X057-X058-X059-X060-X061-X062-X063-X064-X065-X066-X067-X068-X070-X071-X072-X073-X074-X075-X076-X077-X078-X079-X080-X081-X083-X084-X085-X086-X087-X088-X089-X090-X091-X092-X093-X094-X095-X096-X097-X098-X099-X100-X101-X102_1-X102_2-X102_3-X103-X104-X105-X106-X107-X108-X110-X111-X115-X116-X117-X118-X119-X120-X121-X122-X123-X124-X125-X126-X127-X128-X129-X130-X131-X132-X133-X134-X135-X136-X137-X138-X139-X140-X141-X142-X143-X150-X151-X152-X153-X154-X155-X156-X157-X158-X159-X160-X161-X162-X164-X165-X166-X167-X168-X169-X170-X172-X173-X174-X175-X176-X177-X178-X179-X180-X181-X182-X183-X184-X185-X186-X187-X188-X189-X190-X191-X192-X193-X194-X195-X196-X197-X198-X199-X200-X201-X202-X203-X204-X205-X206-X207-X208-X209-X210-X211-X212-X213_1-X213_2-X213_3-X214-X215-X216-X217-X219-X220-X221-X222-X223-X224-X225-X226-X227-X228-X229-X230-X231-X232-X233-X234-X235-X236-X237-X239-X240-X241-X242-X243-X244-X245-X246-X247-X248-X249-X250-X252-X253-X254-X255-X256-X257-X258-X259-X260-X263-X264-X265-X266-X267-X268-X269-X270-X271-X272-X273-X274-X275-X276-X277-X278-X279-X280-X281-X282-X283-X284-X285-X286-X287-X288-X289-X290-X291-X292-X293-X294-X295-X296-X297-X298-X299-X300-X301-X302-X303, subsample = 13926, iters = 100, distance = lt, processors = 20)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.full.trewsummary
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre1.weighted.ave.full.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre1.weighted.std.full.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre1.weighted.phylip.full.dist

unifrac.weighted(tree = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.count_table, groups = P009-P055-P082-P109-P171-P218-P238-P251-P261-P262-X007-X008-X010-X011-X012-X023-X024-X025-X026-X027-X028-X051-X052-X053-X054-X056-X060-X061-X062-X063-X064-X065-X081-X083-X084-X085-X086-X096-X097-X098-X099-X100-X101-X106-X107-X108-X110-X111-X118-X119-X120-X121-X122-X123-X134-X135-X136-X137-X138-X139-X162-X164-X165-X166-X167-X172-X173-X174-X175-X176-X192-X193-X194-X195-X196-X197-X207-X208-X209-X210-X211-X212-X217-X219-X220-X221-X222-X223-X224-X225-X226-X227-X228-X229-X230-X231-X232-X233-X234-X235-X236-X237-X239-X240-X241-X242-X243-X244-X245-X246-X250-X252-X253-X254-X255-X256-X257-X258-X259-X260-X263-X264-X265-X266-X267-X268-X269-X270-X271-X272-X273-X277-X278-X279-X280-X281-X282-X283-X284-X285-X286-X287-X288-X289-X290-X291-X292-X293-X294-X295-X296-X297-X298-X299-X300, subsample = 13926, iters = 100, distance = lt, processors = 20)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.trewsummary
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre1.weighted.ave.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre1.weighted.std.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre1.weighted.phylip.dist

unifrac.weighted(tree = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.count_table, groups = P238-P251-P261-P262-X223-X224-X225-X226-X227-X228-X229-X230-X231-X232-X233-X234-X235-X236-X237-X239-X240-X241-X242-X243-X244-X245-X246-X250-X252-X253-X254-X255-X256-X257-X258-X259-X260-X263-X264-X265-X266-X267-X268-X269-X270-X271-X272-X273-X277-X278-X279-X280-X281-X282-X283-X284-X285-X286-X287-X288-X289-X290-X291-X292-X293-X294-X295-X296-X297-X298-X299-X300, subsample = 13926, iters = 100, distance = lt, processors = 15)
-< Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.trewsummary
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre1.weighted.ave.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre1.weighted.std.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre1.weighted.phylip.dist

unifrac.weighted(tree = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre, count = subsample_OTU_sequences_names_total_rounded_for_mothur.count_table, distance = lt)

list.seqs(count = subsample_OTU_sequences_names_total_non_zero_rounded_for_mothur.count_table)
-> Output File Names:
subsample_OTU_sequences_names_total_non_zero_rounded_for_mothur.accnos

#make a tree where all sequences that have zero counts in subsampled OTU table have been removed and total counts have been rounded
get.seqs(accnos = subsample_OTU_sequences_names_total_non_zero_rounded_for_mothur.accnos, fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.fasta, count = subsample_OTU_sequences_names_total_non_zero_rounded_for_mothur.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.pick.fasta
subsample_OTU_sequences_names_total_non_zero_rounded_for_mothur.pick.count_table



summary.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.pick.fasta, count = subsample_OTU_sequences_names_total_non_zero_rounded_for_mothur.pick.count_table, processors = 20)

clearcut(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.pick.fasta, DNA = T)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.pick.tre

unifrac.weighted(tree = reduced_tree_edit.TRE, count = subsample_OTU_sequences_names_total_non_zero_rounded_for_mothur.pick_edit.count_table, distance = lt)


# make a tree with only incubation dataset and OTUs with rel abundance > 0.1 % per sample (folder tree_0.1)
list.seqs(taxonomy = incubation_0.1_clean_tax_table.txt)
-> Output File Names:
incubation_0.1_clean_tax_table.accnos

get.otus(accnos = incubation_0.1_clean_tax_table.accnos, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.list

list.seqs(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.list)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.accnos

get.seqs(accnos = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.accnos, fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.fasta)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.fasta

clearcut(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.fasta, DNA = T)

# make a tree with only t7 samples and OTUs with rel abundance > 0.1 % per sample type! (folder tree_0.1 -> t7)
list.seqs(taxonomy = t7_0.1_clean_mean_sample_type_ordered_tax_table_for_tree.txt)
-> Output File Names:
t7_0.1_clean_mean_sample_type_ordered_tax_table_for_tree.accnos

get.otus(accnos = t7_0.1_clean_mean_sample_type_ordered_tax_table_for_tree.accnos, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.list)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.list

list.seqs(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.list)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.accnos

get.seqs(accnos = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.abund.0.03.pick.accnos, fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.fasta)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.fasta

clearcut(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.fasta, DNA = T)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.pick.tre

####getting reprsentative sequence for OTU00276###########################################
get.oturep(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.fasta, column = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.dist, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.rep.fasta

############ getting only Planktomycete sequences #####################################

#getting all sequences of certain lineage after classify.seqs
get.lineage(taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.knn.taxonomy, taxon = Planctomycetes, fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.knn.pick_planctos.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick_planctos.count_table


list.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.fasta)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.accnos

#getting only represantative sequences for each plynctomycetean OTU
cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick_planctos.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.nr_v132.knn.pick_planctos.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03, processors = 10)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.opti_mcc.unique_list.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.opti_mcc.unique_list.sensspec

get.oturep(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.fasta, column = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.dist, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.pick_planctos.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.opti_mcc.unique_list.0.03.rep.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.opti_mcc.unique_list.0.03.rep.fasta

make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick.pick_planctos.opti_mcc.unique_list.0.03.rep.count_table)


##############classifying with knn and phycomix sequences included in database########################
classify.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.count_table, reference = new.align, taxonomy = new.tax, cutoff = 85, method = knn, processors = 40)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.knn.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.knn.tax.summary

#removing phyla/lineages that we dont want in our dataset (e.g. Eukaryota, mitochondria etc because the primers are not supposed to match those)
remove.lineage(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.knn.taxonomy, taxon = Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.knn.pick.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.count_table

#OTU picking based on 97% similarity:
cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.knn.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03, processors = 60)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.sensspec

split.abund(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.fasta, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.count_table, cutoff = 3)
Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.rare.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.0.03.rare.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.0.03.abund.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.0.03.rare.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.0.03.abund.count_table

#now we create a new shared file based on the new count table, where the 'singletons' have been removed
make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.0.03.abund.count_table, label = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.shared			#this is the actual OTU table with reads per OTU per sample

#now we can classify our OTUs
classify.otu(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.0.03.abund.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.knn.pick.taxonomy, label = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.0.03.cons.tax.summary	

count.groups(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.shared)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.count.summary		#table with reads per sample 

get.relabund(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.shared)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.relabund          #OTU table with relative abundances

get.lineage(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.knn.pick.taxonomy, taxon = Bacteria;Planctomycetes;Planctomycetacia;Planctomycetales;Planctomycetaceae;Phycimox_1;-Bacteria;Planctomycetes;Planctomycetacia;Planctomycetales;Planctomycetaceae;Phycimox_2;-Bacteria;Planctomycetes;Planctomycetacia;Planctomycetales;Planctomycetaceae;Phycimox_3;)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.knn.pick.pick.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.pick.count_table

##############classifying with wang and phycomix sequences included in database########################
classify.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.count_table, reference = new.align, taxonomy = new.tax, cutoff = 85, method = wang, processors = 40)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.tax.summary

remove.lineage(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.taxonomy, taxon = Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.count_table

cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.01)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.sensspec

split.abund(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.fasta, list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.count_table, cutoff = 3)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.01.rare.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.01.abund.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.0.01.rare.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.0.01.abund.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.0.01.rare.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.0.01.abund.count_table

make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.01.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.0.01.abund.count_table, label = 0.01)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.01.abund.shared

classify.otu(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.01.abund.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.0.01.abund.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.taxonomy, label = 0.01)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.0.01.cons.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.03.abund.0.01.cons.tax.summary

count.groups(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.01.abund.shared)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.01.abund.count.summary

get.relabund(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.01.abund.shared)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.opti_mcc.unique_list.0.01.abund.relabund

####subsampling with mothur before cluster.split#### SOMETHING WENT WRONG HERE WITH THE SIZE PARAMETER!!!
split.abund(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.count_table, cutoff = 3)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.rare.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.abund.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.rare.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.abund.fasta

count.groups(count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.abund.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.abund.count.summary

sub.sample(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.taxonomy, processors = 10)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.subsample.count_table
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.subsample.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.fasta

count.groups(count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.subsample.count_table)

##getting only Phycomix sequences from subsampled dataset##
get.lineage(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.subsample.taxonomy, taxon = Bacteria;Planctomycetes;Planctomycetacia;Planctomycetales;Planctomycetaceae;Phycimox_1;-Bacteria;Planctomycetes;Planctomycetacia;Planctomycetales;Planctomycetaceae;Phycimox_2;-Bacteria;Planctomycetes;Planctomycetacia;Planctomycetales;Planctomycetaceae;Phycimox_3;)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.subsample.pick.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.pick.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.pick.count_table

cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.subsample.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03)
-> Output File Names:
Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.pick.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.pick.opti_mcc.unique_list.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.pick.opti_mcc.unique_list.sensspec

make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.pick.count_table, label = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.pick.opti_mcc.unique_list.shared

classify.otu(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.subsample.pick.taxonomy, label = 0.03)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.pick.opti_mcc.unique_list.0.03.cons.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.subsample.pick.opti_mcc.unique_list.0.03.cons.tax.summary

##getting only Phycomix sequences## THIS IS THE FINAL APPROACH USED! (but recheck count data, if from subsampled or full dataset)
get.lineage(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.subsample.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.taxonomy, taxon = Bacteria;Planctomycetes;Planctomycetacia;Planctomycetales;Planctomycetaceae;Phycimox_1;-Bacteria;Planctomycetes;Planctomycetacia;Planctomycetales;Planctomycetaceae;Phycimox_2;-Bacteria;Planctomycetes;Planctomycetacia;Planctomycetales;Planctomycetaceae;Phycimox_3;)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.fasta
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.pick.count_table

cluster.split(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.fasta, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.taxonomy, splitmethod = classify, taxlevel = 4, cutoff = 0.03, processors = 40)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.dist
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.opti_mcc.unique_list.list
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.opti_mcc.unique_list.sensspec

make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.pick.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.opti_mcc.unique_list.shared

get.relabund(shared = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.opti_mcc.unique_list.shared)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.opti_mcc.unique_list.relabund

classify.otu(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.opti_mcc.unique_list.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.pick.count_table, taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.taxonomy)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.opti_mcc.unique_list.0.03.cons.tax.summary

###phylotype approach###
list.seqs(fasta = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.fasta)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.accnos

get.seqs(accnos = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.pick.pick.accnos, name = stability.trim.contigs.trim.good.names)
-> Output File Names:
stability.trim.contigs.trim.good.pick.names

phylotype(taxonomy = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.taxonomy, name = stability.trim.contigs.trim.good.pick.names)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.tx.sabund
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.tx.rabund
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.tx.list

list.seqs(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.tx.list)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.tx.accnos

get.seqs(accnos = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.tx.accnos, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.pick.count_table

##issued a warning -> recheck compitability of list and count file!
make.shared(list = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.tx.list, count = stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.denovo.vsearch.pick_planctos.pick.pick.count_table)
-> Output File Names:
stability.trim.contigs.trim.good.unique.good.filter.unique.precluster.pick_planctos.new.wang.pick.pick.tx.shared