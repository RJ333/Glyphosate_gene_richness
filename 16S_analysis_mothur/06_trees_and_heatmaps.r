## ggtree?

# plot phylogenetic trees for data subset
ps_dna_water <- subset_samples(mothur_ps4_ra, nucleic_acid == "dna" & 
											  habitat == "water", 
											  prune = TRUE)
											  
											  
ps_dna_water50 <- filter_taxa(ps_dna_water, function (x) {sum(x > 0) > 50}, prune = TRUE)

# can't combine "treeonly" and meta data (as the samples are to be colored)
#plot_tree(ps_dna_water50, "treeonly", label.tips = "taxa_names", text.size = 3)

# one edge is super long, needs to be removed OTU000641
badTaxa = ("Otu000641")
goodTaxa <- setdiff(taxa_names(ps_dna_water50), badTaxa)
ex2 <- prune_taxa(goodTaxa, ps_dna_water50)

plot_tree(ex2, label.tips = "genus", 
			   color = "disturbance",
			   shape = "treatment",
			   size = "Abundance",
			   text.size = 2)

			   
			   
# heatmaps from https://joey711.github.io/phyloseq/plot_heatmap-examples.html
plot_heatmap(ex2, sample.label = "new_day", "genus") + facet_wrap(~treatment)
heatmap(otu_table(ex2))