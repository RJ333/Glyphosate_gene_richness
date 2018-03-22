# this script uses the reads per gene data, which has been processed before
# it generates richness and abundance data from it and,also relation to day 0
# and plots it. 
# gene richness should be at least above 5 once to be taken into account
# mapped reads are not normalized by length of gene nor sequencing depth yet 

library(ggplot2)
library(data.table)

# read in mapping file
gene_reads <- read.delim(file.choose(), header = TRUE)  # contigs_genes_sample_reads.tsv

# generate meta data and merge with mapping file
sample <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10")
days <- c(0, 3, 7, 14, 22, 43, 71, 0, 22, 71)
treatment <- c("glyph", "glyph", "glyph", "glyph", "glyph", "glyph", "glyph", "control", "control", "control")
sample_meta_data <- data.frame(sample, days, treatment)

gene_reads2 <- merge(gene_reads, sample_meta_data, by = "sample")

# calculating richness data per gene
gene_richness <- aggregate( . ~  gene + sample + days + treatment, data = gene_reads2, length)
gene_richness <- gene_richness[with(gene_richness, order(gene)), ]
gene_richness <- gene_richness[,c(1:5)]
colnames(gene_richness)[colnames(gene_richness)=="contig_id"] <- "gene_richness"


# add column with relative richness values with data.table
gene_richness <- setDT(gene_richness)
gene_richness[,gene_richness_relative := gene_richness/gene_richness[days == 0]*100, by = .(gene,treatment)]

# create subsets with and without glyphosate related genes
phn_operon_richness_subset <- subset(gene_richness, grepl("phn[C-P]|sox[A-B]", gene))
ref_richness_subset <- subset(gene_richness, !grepl("phn[C-P]|sox[A-B]", gene))

# plot subsets
ggplot(phn_operon_richness_subset, aes(x = days, y = gene_richness_relative, group = gene, colour = gene)) +
  geom_line(size = 2) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)
  
ggplot(ref_richness_subset, aes(x = days, y = gene_richness_relative, group = gene, colour = gene)) +
  geom_line(size = 2) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)


# generating read abundance data per gene
gene_abundance <- aggregate( . ~  gene + sample + days + treatment, data = gene_reads2, sum)
gene_abundance <- gene_abundance[,c(1:4, 6)]
gene_abundance <- gene_abundance[with(gene_abundance, order(gene)), ]

# add column with relative richness values
gene_abundance <- setDT(gene_abundance)
gene_abundance[,reads_per_gene_relative := reads_per_gene/reads_per_gene[days == 0]*100, by = .(gene,treatment)]

# create subsets
phn_operon_abundance_subset <- subset(gene_abundance, grepl("phn[C-P]|sox[A-B]", gene))
ref_abundance_subset <- subset(gene_abundance, !grepl("phn[C-P]|sox[A-B]", gene))

# plot subsets
ggplot(phn_operon_abundance_subset, aes(x = days, y = reads_per_gene_relative, group = gene, colour = gene)) +
  geom_line(size = 1.5) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)

ggplot(ref_abundance_subset, aes(x = days, y = reads_per_gene_relative, group = gene, colour = gene)) +
  coord_cartesian(ylim = c(0, 400)) +
  geom_line(size = 1.5) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)





