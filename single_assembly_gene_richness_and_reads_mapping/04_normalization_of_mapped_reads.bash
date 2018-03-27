# issue #23 script for gene abundance normalization 

# contig_gene_sample_reads.tsv sums multiple genes on the same contig

# one script, that counts the individual genes per contig per sample
# but the number behind could be used to distinguish, so everything until the ; needs to be included
contig_1409	Prodigal:2.6	CDS	3029	3844	.	+	0	ID=PROKKA_03938;gene=phnE_1;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:A0QQ68;locus_tag=PROKKA_03938;product=Phosphate-import permease protein PhnE
contig_1409	Prodigal:2.6	CDS	3844	4662	.	+	0	ID=PROKKA_03939;gene=phnE_2;inference=ab initio prediction:Prodigal:2.6,similar to AA sequence:UniProtKB:A0QQ68;locus_tag=PROKKA_03939;product=Phosphate-import permease protein PhnE


# output of mapping script

contig	gene	sample	reads_per_gene
contig_15283	araD	A1	237
contig_36332	araD	A1	330
contig_4095	araD	A1	1419
contig_42768	araD	A1	36
contig_60642	araD	A1	53
contig_65521	araD	A1	991

# get gene length per individual genes, summarize gene length per contig, then use above file for number of mapped reads

# benötigte infos


-> sample A1 contig_1409 phnE_1 gene_length ~800 reads per gene per contig in this sample 10000 sequencing depth 40000000  
-> sample A1 contig_1409 phnE_2 gene_length ~800 reads per gene per contig in this sample 10000 sequencing depth 40000000 