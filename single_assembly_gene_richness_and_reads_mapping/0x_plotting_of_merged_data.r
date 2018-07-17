# plots based on prokka_all (contains reads per product per taxonomy) 
# and product_reads2 (contains richness per product)

# data is prepared in 03_calculate_and_plot_gene_richness.R
#and 01_funct_tax_data_merging.bash

# plot taxonomic richness

ggplot(tax_rich_omics, aes(x = new_day, colour = treatment))+
  geom_line(aes(y = richness), linetype = 1) +
  geom_line(aes(y = richness_rel), linetype = 2) +
  facet_wrap(~ treatment, nrow = 2)
  
  
  
# plot gene richness

sarc <- subset(gene_richness2, grepl("arcosin", product2))
sarc <- subset(sarc, grepl("oxidase", product2))

# visualization per product
ggplot(sarc, aes (x = new_day, fill = product2, colour = product2))+
  geom_bar(stat = "identity", position = position_dodge(), aes(y = product_rpm))+
  geom_line(aes(y = gene_richness_relative), size = 1.5)+
  facet_wrap(~treatment, nrow = 2)

# cumulated product's visualization  
ggplot(sarc, aes (fill = product2, colour = product2))+
  geom_bar(stat = "identity", aes(x = new_day - 0.75, y = product_rpm), fill = "black", width = 1)+
  geom_bar(stat = "identity", aes(x = new_day + 0.75, y = gene_richness_relative), fill = "blue", width = 1)+
  facet_wrap(~treatment, nrow = 2) 


phn <- subset(gene_richness2, grepl("phn", gene))
phn <- subset(gene_richness2, grepl("phn[C-N]", gene))

# visualization per product
ggplot(phn, aes (x = new_day, fill = product2, colour = product2))+
  geom_bar(stat = "identity", position = position_dodge(), aes(y = product_rpm))+
  geom_line(aes(y = gene_richness_relative), size = 1.5)+
  facet_wrap(~treatment, nrow = 2)

# cumulated product's visualization  
ggplot(phn, aes (fill = product2, colour = product2))+
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