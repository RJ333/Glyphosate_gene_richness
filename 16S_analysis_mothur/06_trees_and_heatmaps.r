## ggtree?

# plot phylogenetic trees for data subset
ps_dna_water_glyph <- subset_samples(mothur_ps4_ra, nucleic_acid == "dna" & 
													habitat == "water" &
													treatment == "glyph",
													prune = TRUE)
											  

										  
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