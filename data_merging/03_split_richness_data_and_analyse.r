# in the previous scripts we generated the R object (R workspace "data_merging_glyph_single_assemblies.RData") prokka_all
# and all_products_richness. Now each gene/product should be assigned to a group, depending on its richness change
# to determine which group a gene belongs to, certain conditions must be fulfilled

# cell peak value higher than starting value?
# end value
# average value?

# calculate average richness per product

# check the range from highest to lowest value for all products

######################## describe relative position of richness for relevant degrading genes and housekeeping genes

# remove products that occur less than 10 times in 10 samples
all_products_richness10 <- subset(all_products_richness, gene_richness >= 10)

# make glyphosate-treated subsets for peak days
# richness_day0glyph <- subset(all_products_richness10, new_day == 0 & treatment == "glyph")  # all relative richness is 100 % here by definition
richness_day3glyph <- subset(all_products_richness10, new_day == 3 & treatment == "glyph")
richness_day7glyph <- subset(all_products_richness10, new_day == 7 & treatment == "glyph")
richness_day14glyph <- subset(all_products_richness10, new_day == 14 & treatment == "glyph")
richness_day22glyph <- subset(all_products_richness10, new_day == 22 & treatment == "glyph")
richness_day43glyph <- subset(all_products_richness10, new_day == 43 & treatment == "glyph")
richness_day71glyph <- subset(all_products_richness10, new_day == 71 & treatment == "glyph")

# richness_day0control <- subset(all_products_richness10, new_day == 0 & treatment == "control")  # all relative richness is 100 % here by definition
richness_day22control <- subset(all_products_richness10, new_day == 22 & treatment == "control")
richness_day71control <- subset(all_products_richness10, new_day == 71 & treatment == "control")


# add percentiles
# using .bincode instead of cut, as cut require unique values for breaks
richness_day14glyph <- within(richness_day14glyph, percentile <- as.integer(.bincode(gene_richness_relative, quantile(gene_richness_relative, probs=0:100/100), include.lowest=TRUE)))

# soxA can also be cytochrome SoxAX related stuff!

subset(richness_day14glyph, grepl("phosphon|sarc|monomeric", adjusted_products))
subset(richness_day14glyph, grepl("phn|sox", gene))
subset(richness_day14glyph, gene_richness_relative < 100)
subset(richness_day14glyph, percentile > 90)

subset(prokka_all, grepl("soxA", gene))

# list of housekeeping genes

# list of degrading genes (including phn controls phnXW etc and thiO etc)

### create a list to show where the genes are in relation to all other genes, i.e. phnM on position 200, top 25 %
### compare with housekeeping genes, this is the only reference where we know its purpose in this case

