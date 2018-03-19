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
colnames(gene_richness)[colnames(gene_richness)=="contig"] <- "gene_richness"


# add column with relative richness values with data.table
gene_richness <- setDT(gene_richness)
gene_richness[,gene_richness_relative := gene_richness/gene_richness[days == 0]*100, by = .(gene,treatment)]

# create subsets with and without glyphosate related genes
phn_operon_richness_subset <- subset(gene_richness, grepl("phn[C-P]", gene))
phn_operon_richness_subset_glyph <- subset(gene_richness, grepl("phn[C-P]", gene) & treatment == "glyph")
sox_richness_subset <- subset(gene_richness, grepl("sox[A-Z]", gene))
ref_richness_subset <- subset(gene_richness, !grepl("phn[C-P]|sox[A-B]", gene))
ref_richness_subset_glyph <- subset(gene_richness, !grepl("phn[C-P]|sox[A-B]", gene) & treatment == "glyph")

nitrogen_rich_subset <- subset(gene_richness, grepl("nirS|norB|nosZ|nuoF", gene))

# plot subsets
ggplot(nitrogen_rich_subset, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 2) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)

ggplot(phn_operon_richness_subset, aes(x = days, y = gene_richness_relative, group = gene, colour = gene)) +
  geom_line(size = 2) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)

ggplot(sox_richness_subset, aes(x = days, y = gene_richness_relative, group = gene, colour = gene)) +
  geom_line(size = 2) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)
  
ggplot(ref_richness_subset_glyph, aes(x = days, y = gene_richness_relative, group = gene, colour = gene)) +
  geom_line(size = 1) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3.5) +
  geom_line(data = phn_operon_richness_subset_glyph, aes(x = days, y = gene_richness_relative, group = gene),size = 0.8)
  geom_text(data = phn_operon_richness_subset_glyph, label = gene, show.legend = FALSE, size = 3.5)

# generating read abundance data per gene

# reads are not normalized against sequencing depth yet!!!

gene_abundance <- aggregate( . ~  gene + sample + days + treatment, data = gene_reads2, sum)
gene_abundance <- gene_abundance[,c(1:4,6)]
gene_abundance <- gene_abundance[with(gene_abundance, order(gene)), ]

# add column with relative richness values
gene_abundance <- setDT(gene_abundance)
gene_abundance[,reads_per_gene_relative := reads_per_gene/reads_per_gene[days == 0]*100, by = .(gene,treatment)]

# create subsets
phn_operon_abundance_subset <- subset(gene_abundance, grepl("phn[C-P]", gene))
sox_abundance_subset <- subset(gene_abundance, grepl("sox[A-Z]", gene))

# ref_abundance_subset <- subset(gene_abundance, !grepl("phn[C-P]|sox[A-B]", gene))

# plot subsets
ggplot(sox_abundance_subset, aes(x = days, y = reads_per_gene, group = gene, colour = gene)) +
  geom_line(size = 1.5) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)

  
  
ggplot(ref_abundance_subset, aes(x = days, y = reads_per_gene_relative, group = gene, colour = gene)) +
  coord_cartesian(ylim = c(0,400)) +
  geom_line(size = 1.5) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)





