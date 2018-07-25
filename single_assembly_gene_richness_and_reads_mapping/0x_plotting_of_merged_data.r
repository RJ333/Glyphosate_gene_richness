# plots based on prokka_all (contains reads per product per taxonomy) 
# and product_reads2 (contains richness per product)

# data is prepared in 03_calculate_and_plot_gene_richness.R
#and 01_funct_tax_data_merging.bash


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

# relative abundance and relative richness are not correlated?? (phnC, phnD, phnH
# more correlation for phnE,
# differences in reaction time? because of detection limit per contig?
# back to starting niveau for abundance? or richness? elevated richness over time?
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

# gallaeci besitzt soxBCD gene auf vier contigs und phnCDE 1 contig H noch ein contig

# methylo shikimate dehydrogenase? aroE, 1.1.1.25

#Pseudomonas hat phn[C-N] und alle sox gene
#Hoeflea soxBCD, alle phn außer JK

#plots aller relevanten OTUS
#Anteil von phn Genen und sox Genen ohne taxonomische Zuordnung?
table(is.na(phn_op$genus))
FALSE  TRUE 
  126   958 
  
table(is.na(sox_all$genus))
FALSE  TRUE 
   71   415 

#taxonomie sollte nicht überbewerttet werden

#ist richness empfindlicher als abundanz? wie teste ich das im vergleich
#zu "normalen" Genen?
   
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