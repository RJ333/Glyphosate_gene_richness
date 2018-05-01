# data is collected in prokka_all in workspace data_merging_single_assemblies.RData
library(ggplot2)
library(data.table)


# check out phn genes normalized abundance
degrade.set <- droplevels(subset(prokka_all, grepl("phn", gene) | grepl("sarcosine", adjusted_products)))
levels(degrade.set$gene)
# it includes some phn genes and products which we don't want, as they are not degrading glyphosate
degrade.set <- droplevels(subset(degrade.set, !(gene %in% c("phnX","phnW","phnS","phnT","phnU","phnV","phnR","phnA"))))
degrade.set <- droplevels(subset(degrade.set, !(grepl("carbam", adjusted_products))))
str(degrade.set)
head(degrade.set)
# combine gene and product name for addressing in plot
degrade.set$ident <- do.call(paste, c(degrade.set[c("gene","adjusted_products")], sep = "_"))
degrade.set$ident <- as.factor(degrade.set$ident)
degrade.set$ident <- factor(degrade.set$ident, levels(degrade.set$ident)[order(degrade.set$ident)])

ggplot(degrade.set, aes(x = new_day, y = rpm, group = ident, fill = ident)) +
  geom_bar(width = 1.7, stat = "identity")+
  scale_fill_discrete(name  ="Gene products") +
  xlab("Days")+
 # theme(legend.position="none")+
  ylab("normalized reads")+
  facet_wrap(~treatment)
  
 

degrade.set2 <- degrade.set[,c(4, 11, 14, 15, 19)]
 
# calculating richness data per gene
degrade.set2_richness <- aggregate( . ~  adjusted_products + new_day + treatment, data = degrade.set2, length)
degrade.set2_richness <- degrade.set2_richness[with(degrade.set2_richness, order(adjusted_products)), ]
degrade.set2_richness <- degrade.set2_richness[,c(1:4)]
colnames(degrade.set2_richness)[colnames(degrade.set2_richness)=="contig_id"] <- "gene_richness"
degrade.set2_richness <- setDT(degrade.set2_richness)
degrade.set2_richness[,gene_richness_relative := gene_richness/gene_richness[new_day == 0]*100, by = .(adjusted_products, treatment)]

ggplot(degrade.set2_richness, aes(x = new_day, y = gene_richness)) +
  geom_line(aes(group = adjusted_products, colour = adjusted_products), size = 1.2) +
  #geom_text(aes(label = adjusted_products), show.legend = FALSE, size = 3) +
  #geom_line(data = tax_rich_omics, aes(x = days, y = richness_rel))+
  #geom_text(data = tax_rich_omics, aes(x = days, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
ggplot(degrade.set2_richness, aes(x = new_day, y = gene_richness_relative)) +
  geom_line(aes(group = adjusted_products, colour = adjusted_products), size = 1.2) +
  #geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  #geom_line(data = tax_rich_omics, aes(x = days, y = richness_rel))+
  #geom_text(data = tax_rich_omics, aes(x = days, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)