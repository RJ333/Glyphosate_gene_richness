# this is a first test of how to create ordination plots 
# for glyphosate data in R phyloseq

# use the relative abundance phyloseq object, no transformation needed

# Set the working dir with mothur files in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("mothur_glyph.RData")

# Load packages into session, and print package version
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)



# this is copied from a tutorial from:
# https://joey711.github.io/phyloseq/plot_ordination-examples.html

# subset_samples(GlobalPatterns, SampleType=="Ocean")




# Remove OTUs that do not show appear more than 5 times in more than half the samples
# wh0 = genefilter_sample(GP, filterfun_sample(function(x) x > 5), A=0.5*nsamples(GP))
# GP1 = prune_taxa(wh0, GP)



ps_dna_water <- subset_samples(mothur_ps4_ra, nucleic_acid == "dna" & 
											  habitat == "water")
# turn days into factors for plotting
sample_data(ps_dna_water)$days <- as.factor(sample_data(ps_dna_water)$days)
GP1 = ps_dna_water

# (1) Just OTUs

# Let’s start by plotting just the OTUs, and shading the points by Phylum. 
# Note that even in our “trimmed” dataset there are ntaxa(GP1)= 204 OTUs.

GP.ord <- ordinate(GP1, "NMDS", "bray")
p1 <- plot_ordination(GP1, GP.ord, type="taxa", color="phylum", title="taxa") +
			facet_wrap(~phylum, 3)
		
print(p1)

# (2) Just samples

# Next, let’s plot only the samples, and shade the points by “SampleType” 
# while also modifying the shape according to whether they are human-associated. 
# There are a few additional ggplot2 layers added to make the plot even nicer…

p2 = plot_ordination(GP1, GP.ord, type="sample", color="days", shape="treatment") 
p2 + geom_polygon(aes(fill=nucleic_acid)) + geom_point(size=5) + ggtitle("samples")


##### below not possible without phytree object!!

# In this section I loop through different method parameter options to the 
# plot_ordination function, store the plot results in a list, and then plot 
# these results in a combined graphic using ggplot2.

require(plyr)
dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
        ordi = ordinate(physeq, method = i, distance = dist)
        plot_ordination(physeq, ordi, "samples", color = "days", shape = "treatment")
}, GP1, dist)
names(plist) <- ord_meths

# The following chunk will extract the data from each of those individual plots, 
# and put it back together in one big data.frame suitable 
# for including all plots in one graphic.

pdataframe = ldply(plist, function(x){
    df = x$data[, 1:2]
    colnames(df) = c("Axis_1", "Axis_2")
    return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

# Now that all the ordination results are combined in one data.frame, called 
# pdataframe, we can use this to make a standard faceted ggplot scatterplot.

p = ggplot(pdataframe, aes(Axis_1, Axis_2, color=SampleType, shape=human, fill=SampleType))
p = p + geom_point(size=4) + geom_polygon()
p = p + facet_wrap(~method, scales="free")
p = p + scale_fill_brewer(type="qual", palette="Set1")
p = p + scale_colour_brewer(type="qual", palette="Set1")
p