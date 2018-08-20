# display abundance of phn gene per contig
library(ggplot2)


# this is for gene abbreviations

phn_genes <- c("phnC", "phnD", "phnE", "phnF", "phnG", "phnH", "phnI", "phnJ", "phnK", "phnL", "phnM",
				"phnN", "phnP")
				
phn_plot_preparation <- function(selection){
  gene_subset <- subset(prokka_all, gene == selection)
  gene_agg <- aggregate(gene ~ sample + new_day + treatment, gene_subset, length)
  list(gene_subset, gene_agg)
}

for (phn_gene in phn_genes) {
  result_list <- phn_plot_preparation(phn_gene)

print(ggplot(result_list[[1]]) +
  geom_col(aes(x = new_day, y = product_rpm/5, group = sample), width = 1.5, alpha = 0.5) +
  geom_violin(aes(x = new_day, y = product_rpm, group = sample), size = 0.5, alpha = 0.5) +
  geom_text(data = result_list[[2]], size = 5, colour = "grey50", 
    aes(x = new_day, y = -2, label = paste0("n=", gene))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Cumulated rpm"))+
  ggtitle(phn_gene) +
  labs(y = "Rpm per contig", 
	   x = "Days after glyphosate addition") +
  theme_bw()+
	theme(axis.text = element_text(size = 22))+
	theme(axis.text.x = element_text(angle = 0, vjust = 0.5))+
	theme(axis.title = element_text(size = 25)) +
	theme(axis.title.x = element_text(angle = 0, vjust = 0, margin = margin(t = 15, r = 0, b = 0, l = 0)))+
	theme(panel.grid.major = element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor = element_line(colour = NA, size = 0.5))+
	theme(legend.position = "none") +
  facet_wrap(~treatment, nrow = 2))
  ggsave(file = paste(phn_gene,"per_contig_plot.png"), width = 14, height = 12)
}

# this is for gene product names

sox_genes <- c("monomeric@sarcosine@oxidase", 
			   "sarcosine@oxidase%2c@delta@subunit@family", 
			   "sarcosine@oxidase%2c@gamma@subunit@family", 
			   "sarcosine@oxidase@subunit@beta")

sox_plot_preparation <- function(selection){
  gene_subset <- subset(prokka_all, product2 == selection)
  gene_agg <- aggregate(product2 ~ sample + new_day + treatment, gene_subset, length)
  list(gene_subset, gene_agg)
}
# the product names are too long for the ggplot savings...
# so we use a counter j to not overwrite the plots
j <- 0

for (sox_gene in sox_genes) {
  result_list <- sox_plot_preparation(sox_gene)
  
j <- j + 1

print(ggplot(result_list[[1]]) +
  geom_col(aes(x = new_day, y = product_rpm/5, group = sample), width = 1.5, alpha = 0.5) +
  geom_violin(aes(x = new_day, y = product_rpm, group = sample), size = 0.5, alpha = 0.5) +
  geom_text(data = result_list[[2]], size = 5, colour = "grey50", 
    aes(x = new_day, y = -2, label = paste0("n=", product2))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_y_continuous(sec.axis = sec_axis(~.*5, name = "Cumulated rpm"))+
  ggtitle(sox_gene) +
  labs(y = "Rpm per contig", 
	   x = "Days after glyphosate addition") +
  theme_bw()+
	theme(axis.text = element_text(size = 22))+
	theme(axis.text.x = element_text(angle = 0, vjust = 0.5))+
	theme(axis.title = element_text(size = 25)) +
	theme(axis.title.x = element_text(angle = 0, vjust = 0, margin = margin(t = 15, r = 0, b = 0, l = 0)))+
	theme(panel.grid.major = element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor = element_line(colour = NA, size = 0.5))+
	theme(legend.position = "none") +
  facet_wrap(~treatment, nrow = 2))
  ggsave(file = paste0("sox",j,".png"), device = "png", width = 14, height = 12)
}