# plots based on prokka_all (contains reads per product per taxonomy) 
# and product_reads2 (contains richness per product)

# data is prepared in 03_calculate_and_plot_gene_richness.R
#and 01_funct_tax_data_merging.bash

library(ggplot2)
library(dplyr)
# plot taxonomic richness

ggplot(tax_rich_omics, aes(x = days, colour = treatment))+
  geom_line(aes(y = richness), linetype = 1) +
  geom_line(aes(y = richness_rel), linetype = 2) +
  facet_wrap(~ treatment, nrow = 2)
  
  
  
# plot gene richness
sarc <- subset(gene_richness2, grepl("arcosin", product2))
sarc <- subset(sarc, grepl("oxidase", product2))
sarc <- subset(sarc, grepl("delta", product2))

# visualization per product
ggplot(sarc, aes (x = new_day, colour = product2, fill = contig_id))+
  geom_bar(stat = "identity", position = position_dodge(), aes(y = product_rpm))+
  #geom_line(aes(y = gene_richness_relative), size = 1.5)+
  #geom_line(data = tax_rich_omics, aes(x = days, y = richness_rel, colour = "black"), linetype = 2, size = 2, alpha = 0.7) +
  theme(legend.position="none")+
  facet_wrap(~treatment, nrow = 2)

# cumulated product's visualization  
ggplot(sarc, aes(colour = product2))+
  geom_bar(stat = "identity", aes(x = new_day - 0.75, y = product_rpm), fill = "black", width = 1)+
  geom_bar(stat = "identity", aes(x = new_day + 0.75, y = gene_richness_relative), fill = "blue", width = 1)+
  geom_line(data = tax_rich_omics, aes(x = days, y = richness_rel, colour = "black"), linetype = 2) + 
  facet_wrap(~treatment, nrow = 2) 


phn_op <- subset(prokka_all, treatment == "glyph" & grepl("phn[C-N]", gene))
sox_all <- subset(prokka_all, treatment == "glyph" & grepl("sarcosine@oxidase", product2))


test <- sox_all[with(sox_all, order(sox_all$product2, sox_all$sample, -sox_all$average_coverage)),]

test <- table(prokka_all$average_coverage) 
# with this we can check what average coverage is neccessary to
# generate a contig and check, how close our interesting genes are to that
# lowest value  for phn_op is 2.296, for sox_all is 2.55, 
# but we have at least ~7500 contigs with lower average coverage
sum(head(test, 4600))

phn_op_rich <- subset(gene_richness2, treatment == "glyph" & grepl("phn[C-N]", gene))
sox_all_rich <- subset(gene_richness2, treatment == "glyph" & grepl("sarcosine@oxidase", product2))



#Anteil von phn Genen und sox Genen ohne taxonomische Zuordnung?
table(is.na(phn_op$genus))
FALSE  TRUE 
  126   958 
  
table(is.na(sox_all$genus))
FALSE  TRUE 
   71   415 


#------------------------------------------------------------------ assessment of single copy and housekeeping genes

   
# what genes do we have from Wu et al marker genes? 
# additionally Mads Albertsen et al, "recovery of genomes", 2013, supplementary Table S4
levels(droplevels(subset(gene_richness2, grepl("rpl", gene))$gene))
levels(droplevels(subset(gene_richness2, grepl("rps", gene))$gene))   

table(droplevels(subset(gene_richness2, grepl("rps", gene))$gene))
table(droplevels(subset(gene_richness2, grepl("rpl", gene))$gene))

table(droplevels(subset(prokka_all, grepl("rpl", gene))$gene))
table(droplevels(subset(prokka_all, grepl("rpm", gene))$gene))
table(droplevels(subset(prokka_all, grepl("rps", gene))$gene))


table(droplevels(subset(prokka_all, grepl("era", gene))$product2))
table(droplevels(subset(prokka_all, grepl("trna", product2) & grepl("class", product2))$product2))

table(droplevels(subset(prokka_all, grepl("trna", product2) & grepl("ligase", product2))$product2))

#translation table product2 gene
gene_product2 <- unique(gene_richness2[, c(1, 2)])

rpl_A1 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@l", product2) & sample == "A1")$product2)))
rpl_A2 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@l", product2) & sample == "A2")$product2)))
rpl_A3 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@l", product2) & sample == "A3")$product2)))
rpl_A4 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@l", product2) & sample == "A4")$product2)))
rpl_A5 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@l", product2) & sample == "A5")$product2)))
rpl_A6 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@l", product2) & sample == "A6")$product2)))
rpl_A7 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@l", product2) & sample == "A7")$product2)))
rpl_B8 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@l", product2) & sample == "B8")$product2)))
rpl_B9 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@l", product2) & sample == "B9")$product2)))
rpl_B10 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@l", product2) & sample == "B10")$product2)))

rpl <- Reduce(function(x,y) merge(x = x, y = y, by = "Var1"), 
  list(rpl_A1, rpl_A2, rpl_A3, rpl_A4, rpl_A5, rpl_A6, rpl_A7, rpl_B8, rpl_B9, rpl_B10))
names(rpl) <- c("product2", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10")
rpl <- merge(rpl, gene_product2, by = "product2")
rpl[with(rpl, order(rpl$product2, as.character(rpl$gene), decreasing = TRUE)),]
rpl %>% arrange(product2, desc(gene)) %>% distinct(product2, A1, A2, A3, A4, A5, .keep_all = TRUE)


rps_A1 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@s", product2) & sample == "A1")$product2)))
rps_A2 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@s", product2) & sample == "A2")$product2)))
rps_A3 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@s", product2) & sample == "A3")$product2)))
rps_A4 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@s", product2) & sample == "A4")$product2)))
rps_A5 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@s", product2) & sample == "A5")$product2)))
rps_A6 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@s", product2) & sample == "A6")$product2)))
rps_A7 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@s", product2) & sample == "A7")$product2)))
rps_B8 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@s", product2) & sample == "B8")$product2)))
rps_B9 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@s", product2) & sample == "B9")$product2)))
rps_B10 <- as.data.frame(table(droplevels(subset(prokka_all, grepl("ribosomal@protein@s", product2) & sample == "B10")$product2)))

rps <- Reduce(function(x,y) merge(x = x, y = y, by = "Var1"), 
  list(rps_A1, rps_A2, rps_A3, rps_A4, rps_A5, rps_A6, rps_A7, rps_B8, rps_B9, rps_B10))
names(rps) <- c("product2", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10")
rps <- merge(rps, gene_product2, by = "product2")
rps[with(rps, order(rps$product2, as.character(rps$gene), decreasing = TRUE)),]
rps %>% arrange(product2, desc(gene)) %>% distinct(product2, A1, A2, A3, A4, A5, .keep_all = TRUE)

#-----------------------------------------------------------------------



#fragen?

#ist richness empfindlicher als abundanz? wie teste ich das im vergleich
#zu "normalen" Genen?

# relative abundance and relative richness are not correlated?? (phnC, phnD, phnH
# more correlation for phnE,
# there is an excel table on that findings "richness_abundance comparison"
# differences in reaction time? because of detection limit per contig?
# back to starting niveau for abundance? or richness? elevated richness over time?

# To do
#plots aller relevanten OTUS

# taxonomy + gene:

# gallaeci besitzt soxBCD gene auf vier contigs und phnCDE 1 contig H noch ein contig
# methylo shikimate dehydrogenase? aroE, 1.1.1.25
#Pseudomonas hat phn[C-N] und alle sox gene
#Hoeflea soxBCD, alle phn außer JK
#taxonomie sollte nicht überbewerttet werden





as.data.frame(table(droplevels(phn_op$genus)))

as.data.frame(table(droplevels(sox_op$genus)))
 as.data.frame(table(droplevels(subset(phn_op, genus == "Pseudomonas")$gene)))


phn <- subset(gene_richness2, grepl("phnN", gene))
#phn <- subset(gene_richness2, grepl("phn[C-N]", gene))

# visualization per product



# cumulated product's visualization  
ggplot(sarc, aes (fill = product2, colour = product2))+
  geom_bar(stat = "identity", aes(x = new_day - 0.75, y = product_rpm), fill = "black", width = 1)+
  geom_bar(stat = "identity", aes(x = new_day + 0.75, y = gene_richness_relative), fill = "blue", width = 1)+
  facet_wrap(~treatment, nrow = 2) 
  
  
###### plots for quick investigation: relative richness containing tax richness
  
ggplot(phn_operon_sox_richness_subset, aes(x = new_day, y = gene_richness2)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  #geom_line(data = tax_rich_omics, aes(x = new_day, y = richness_rel))+
  #geom_text(data = tax_rich_omics, aes(x = new_day, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
ggplot(sox_richness_subset, aes(x = new_day, y = gene_richness2_relative)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  #geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  #geom_line(data = tax_rich_omics, aes(x = new_day, y = richness_rel))+
  #geom_text(data = tax_rich_omics, aes(x = new_day, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
ggplot(housekeeping_richness, aes(x = new_day, y = gene_richness2_relative)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  #geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  geom_line(data = tax_rich_omics, aes(x = new_day, y = richness_rel), size = 1.3)+
  #geom_text(data = tax_rich_omics, aes(x = new_day, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(new_genes, aes(x = new_day, y = gene_richness2)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  #geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  geom_line(data = tax_rich_omics, aes(x = new_day, y = richness_rel))+
  #geom_text(data = tax_rich_omics, aes(x = new_day, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
  
# plot species abundance
Paracoccus <- subset(prokka_all, genus == "Paracoccus")
ggplot(Paracoccus, aes(x = new_day, y = kaiju_rpm, fill = treatment)) +
geom_bar(stat = "identity") +
facet_wrap(~ treatment)


# plot gene abundance with taxonomic annotation
phnJ <- subset(prokka_all, gene == "phnJ")

ggplot(phnJ, aes(x = new_day, y = product_rpm, fill = genus)) +
geom_bar(stat = "identity") +
facet_wrap(~ treatment)


# plot multiple genes or products
phn <- subset(prokka_all, grepl("phn", gene))

ggplot(phn, aes(x = new_day, y = product_rpm, fill = order, colour = gene)) +
#geom_bar(stat = "identity") +
geom_bar(stat = "identity", position = position_dodge(), size = 2) +
facet_wrap(~ treatment)

# use multiple subsets to refine results
sarcos <- subset(prokka_all, grepl("sarcosine", product2))
sarcos <- subset(sarcos, grepl("oxidase", product2))

ggplot(sarcos, aes(x = new_day, y = product_rpm, fill = product2, colour = genus)) +
geom_bar(stat = "identity") +
#geom_bar(stat = "identity", position = position_dodge(), size = 1.5) +
facet_wrap(~ treatment, nrow = 2)

# search for specific gene and specific organism
sarcos_methylo <- subset(sarcos, grepl("Methylo", genus))

ggplot(sarcos_methylo, aes(x = new_day, y = product_rpm, fill = product2, colour = genus)) +
geom_bar(stat = "identity") +
#geom_bar(stat = "identity", position = position_dodge(), size = 1.5) +
facet_wrap(~ treatment, nrow = 2)