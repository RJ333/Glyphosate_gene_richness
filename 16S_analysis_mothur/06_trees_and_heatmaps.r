# install ggtree
source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")

# fine, but gives me a lot of sample points which I don't need
plot_tree(ps, label.tips = "genus")
# this is also fine, but the sequence header alone does not help me
plot_tree(ps, "treeonly", label.tips = "taxa_names")
# this would be what I want:
plot_tree(ps, "treeonly", label.tips = "genus")
# but gives

ps <- ps_water_glyph2 


require(ggplot2)
require(ggtree)
require(phyloseq)
require(scales)
# plot phylogenetic trees for data subset

ps_water_glyph <- subset_samples(mothur_ps4_ra,	habitat == "water" &
													treatment == "glyph",
													prune = TRUE)

badTaxa = c("Otu000641", "Otu000674")
goodTaxa <- setdiff(taxa_names(ps_water_glyph), badTaxa)
ps_water_glyph2 <- prune_taxa(goodTaxa, ps_water_glyph)
													
ps_water_glyph2 <- filter_taxa(ps_water_glyph2, function (x) {sum(x > 0) > 15}, prune = TRUE)											
			
tree_plot <- plot_tree(ps_water_glyph2, 
								label.tips = "genus", 
							  color = "disturbance",
							   #size = "Abundance",
							   text.size = 2)

plot_folder <- "/data/projects/glyphosate/plots/R/trees_16S/"	

ggsave(tree_plot, file = paste(plot_folder, "water_glyph_16.svg", 
								  sep = ""),
								  height = 20,
								  width = 13)


## pull the tree out of the phyloseq object
GP.tree <- phy_tree(ps_water_glyph2)

## make a trait annotation for those 50 sequences
x <- rep(c("T","F"), 10000 * c(0.5, 0.5))
x <- sample(x, 251) 
trait <- c(x)
trait <- as.factor(trait) # doesn't work as logical

## put the trait in a dataframe
require(tibble)
x <- data_frame(label = GP.tree$tip.label, trait = trait)

## convert the phylo object to a treeio::treedata object
require(geiger)
GP.tree <- treedata(phy = GP.tree, data = GP.tree$tip.label)

## add the annotation
GP.tree <- full_join(GP.tree, x, by="label")

## take a look at the tree
ggtree(GP.tree) + 
	geom_text2(aes(label=label), size = 2) + 
	geom_tiplab(aes(label = Genus), hjust = -0.3)

ps_genus <- tax_glom(ps_water_glyph2, "Genus", NArm = TRUE)
	
	

								  
								  
								  
								  
#scp -i /drives/d/ssh/denbi.key centos@193.196.20.111:/data/projects/glyphosate/plots/R/trees_16S/* /mnt/d/denbi/chandler/trees_16S/

# tree with otus marked only during disturbance = high?
						   
# there might be a conflict between phyloseq and ggtree? but function needs "psmelt()" from phyloseq...
# detach("package:phyloseq", unload = TRUE)
											  
draw_phyloseq_ggtree <- function(phyloseq) {	
    tree <- phyloseq@phy_tree
    p <- ggtree(tree, ladderize = F)
    p <- p + geom_text(subset=.(!isTip), aes(label=label), hjust=-.2, size=4)
    
    dd <- psmelt(phyloseq)
    dd <- dd[dd$Abundance > 0, ]
    data <- merge(p$data, dd, by.x="label", by.y="OTU")
    
    spacing <- 0.02
    idx <- with(data, sapply(table(node)[unique(node)], function(i) 1:i)) %>% unlist
    hjust <- spacing * idx * max(data$x)
    data$xdodge <- data$x + hjust
    
    p <- p + geom_point(data=data, aes(x=xdodge, color=days,
                            shape=family, size=Abundance), na.rm=T) + 
        theme(legend.position="right") + scale_size_continuous(trans=log_trans(5))
    
    d2 <- unique(data[, c("x", "y", "genus")])
    p + geom_text(data=d2, aes(label=genus), hjust=-.3, na.rm=T, size=4)
}



draw_phyloseq_ggtree(ps_dna_water_glyph2)
										  
ps_dna_water_glyph5 <- filter_taxa(ps_dna_water_glyph, function (x) {sum(x > 0) > 5}, prune = TRUE)


# one edge is super long, needs to be removed OTU000641
badTaxa = c("Otu000641", "Otu000674")
goodTaxa <- setdiff(taxa_names(ps_dna_water_glyph5), badTaxa)
ps_dna_water_glyph5 <- prune_taxa(goodTaxa, ps_dna_water_glyph5)

plot_tree(ex2, label.tips = "genus", 
			   color = "disturbance",
			   shape = "treatment",
			   size = "Abundance",
			   text.size = 2)

# can't combine "treeonly" and meta data (as the samples are to be colored)
plot_tree(ps_dna_water_glyph5, "treeonly", label.tips = "genus", text.size = 3, color = "new_day")			   

# agglomerate by difference (phylogeny)
h1 = 0.2
ps_dna_water_glyph3_0.2 <- tip_glom(ps_dna_water_glyph3, h = h1)

# agglomerate by rank (taxonomy)
ps_dna_water_glyph3_genus <- tax_glom(ps_dna_water_glyph3, "genus", NArm = FALSE)
taxa_names(ps_dna_water_glyph3_genus) <- tax_table(ps_dna_water_glyph3_genus)[,6]

# after glomming, only taxa_names -> Otu000... can be used as label, not genus or so
a <- plot_tree(ps_dna_water_glyph3_0.2, "treeonly", label.tips = "taxa_names", text.size = 3)
b <- plot_tree(ps_dna_water_glyph3_genus, "treeonly", label.tips = "taxa_names", text.size = 3)

# group plots together
grid.arrange(nrow = 1, a, b)

# heatmaps from https://joey711.github.io/phyloseq/plot_heatmap-examples.html

# plot sample similarity (disturbance high e.g.) and phylogeny
plot_heatmap(ex2, sample.label = "new_day", "genus") + facet_wrap(~treatment)
heatmap(otu_table(ex2))