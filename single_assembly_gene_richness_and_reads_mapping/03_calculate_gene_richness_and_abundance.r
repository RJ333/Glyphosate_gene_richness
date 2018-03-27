# this script uses the reads per gene data, which has been processed before
# it generates richness and abundance data from it and,also relation to day 0
# and plots it. 
# gene richness should be at least above 5 once to be taken into account
# mapped reads are not normalized by length of gene nor sequencing depth yet 

library(ggplot2)
library(data.table)

# read in mapping file
gene_reads <- read.delim(file.choose(), header = TRUE)  # contigs_genes_sample_reads.tsv

# read in taxonomic richness file
tax_rich <- read.delim(file.choose(), header = TRUE)  # tax_richness.txt

# plot taxonomic richness

ggplot(tax_rich, aes(x = days, colour = treatment))+
  geom_line(aes(y = richness), linetype = 1) +
  geom_line(aes(y = richness_rel), linetype = 2) +
  
  facet_wrap(~ treatment, nrow = 2)


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
# check if absolute richness values are above > 5 at least once in treatment
# ydiF, malF, ktrB from "phnM_like" too small
# norB from "nitrogen" too small
# pphA, phnR, phnA, phnT from "more_phosphonate"
# chtA not found at all, chiA always 2
# all glycosyl hydrolases too snall

phnM_like <- subset(gene_richness, grepl("phnM|nuoF|mntB|araD|tfdB|artI|arfA|qedA|mlhB|purB", gene))
phn_operon_richness_subset <- subset(gene_richness, grepl("phn[C-N]", gene))
phn_operon_richness_subset_glyph <- subset(gene_richness, grepl("phn[C-N]", gene) & treatment == "glyph")
sox_richness_subset <- subset(gene_richness, grepl("sox[A-Z]", gene))
not_phn_richness_subset <- subset(gene_richness, !grepl("phn[C-P]", gene))
ref_richness_subset_glyph <- subset(gene_richness, !grepl("phn[C-P]|sox[A-B]", gene) & treatment == "glyph")
housekeeping_richness <- subset(gene_richness, grepl("gyrA|gyrB|purB|rpoC|recG", gene))
nitrogen_rich_subset <- subset(gene_richness, grepl("nirS|nosZ", gene))
p_starve_rich_subset <- subset(gene_richness, grepl("pho|pst|psp", gene))
katalase_rich <- subset(gene_richness, grepl("kat", gene))
more_phosphonate <- subset(gene_richness, grepl("phn[USWX]|thiO", gene))
#chitinase <- subset(gene_richness, grepl("chiA|chtA", gene))
#glycosyl_hydrolases <- subset(gene_richness, grepl("xynB|bcsZ|celB|cenC", gene))
metal <- subset(gene_richness, grepl("merA|czcD", gene))  # differing in control, very similar in treatment, also to phnE
monooxy <- subset(gene_richness, grepl("tmoS", gene))  # could be useful
# plot subsets

ggplot(chitinase, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
ggplot(glycosyl_hydrolases, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)


###### plot containing tax richness
  
ggplot(phn_operon_richness_subset, aes(x = days, y = gene_richness_relative)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  #geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  geom_line(data = tax_rich, aes(x = days, y = richness_rel))+
  #geom_text(data = tax_rich, aes(x = days, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

# produce subset for average and SE calculation
write.csv(phn_operon_richness_subset, file = "phn_op_richness.csv") 
write.csv(not_phn_richness_subset, file = "not_phn_op_richness.csv") 
  
########

ggplot(monooxy, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(katalase_rich, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(more_phosphonate, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
ggplot(p_starve_rich_subset, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(phnM_like, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)  
  
ggplot(housekeeping_richness, aes(x = days, y = gene_richness_relative, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(nitrogen_rich_subset, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 2) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)

ggplot(phn_operon_richness_subset_glyph, aes(x = days, y = gene_richness_relative, group = gene, colour = gene)) +
  geom_line(size = 1.0) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(sox_richness_subset, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 2) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)
  
ggplot(ref_richness_subset_glyph, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3.5) +
  geom_line(data = phn_operon_richness_subset_glyph, aes(x = days, y = gene_richness_relative, group = gene),size = 0.8)
  geom_text(data = phn_operon_richness_subset_glyph, label = gene, show.legend = FALSE, size = 3.5)

# generating read abundance data per gene

# reads are not normalized against sequencing depth yet!!!

gene_abundance <- aggregate( . ~  gene + sample + days + treatment, data = gene_reads2, sum)
gene_abundance <- gene_abundance[,c(1:4, 6)]
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
  coord_cartesian(ylim = c(0, 400)) +
  geom_line(size = 1.5) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)





