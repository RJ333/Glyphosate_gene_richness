# in the previous scripts we generated the R object (R workspace "data_merging_glyph_single_assemblies.RData") prokka_all
# and all_products_richness. Now each gene/product should be assigned to a group, depending on its richness change
# to determine which group a gene belongs to, certain conditions must be fulfilled

# cell peak value higher than starting value?
# end value
# average value?

# calculate average richness per product

# check the range from highest to lowest value for all products

######################## describe relative position of richness for relevant degrading genes and housekeeping genes

# remove products that occur less than 2 times per sample
all_products_richness2 <- subset(all_products_richness, gene_richness >= 2)

# make glyphosate-treated subsets for peak days
richness_day0glyph <- subset(all_products_richness2, new_day == 0 & treatment == "glyph")  # all relative richness is 100 % here by definition
richness_day3glyph <- subset(all_products_richness2, new_day == 3 & treatment == "glyph")
richness_day7glyph <- subset(all_products_richness2, new_day == 7 & treatment == "glyph")
richness_day14glyph <- subset(all_products_richness2, new_day == 14 & treatment == "glyph")
richness_day22glyph <- subset(all_products_richness2, new_day == 22 & treatment == "glyph")
richness_day43glyph <- subset(all_products_richness2, new_day == 43 & treatment == "glyph")
richness_day71glyph <- subset(all_products_richness2, new_day == 71 & treatment == "glyph")

richness_day0control <- subset(all_products_richness2, new_day == 0 & treatment == "control")  # all relative richness is 100 % here by definition
richness_day22control <- subset(all_products_richness2, new_day == 22 & treatment == "control")
richness_day71control <- subset(all_products_richness2, new_day == 71 & treatment == "control")


#for (i in c("3glyph", "7glyph", "14glyph", "22glyph", "43glyph", "71glyph", "22control", "71control")) {
#print(eval(parse(text=paste0("richness_day",i))))
#}

# add percentiles for relative richness (per sample!)

# using .bincode instead of cut, as cut require unique values for breaks
richness_day0glyph <- within(richness_day0glyph, percentile <- as.integer(.bincode(gene_richness_relative, quantile(gene_richness_relative, probs=0:100/100), include.lowest=TRUE)))
richness_day3glyph <- within(richness_day3glyph, percentile <- as.integer(.bincode(gene_richness_relative, quantile(gene_richness_relative, probs=0:100/100), include.lowest=TRUE)))
richness_day7glyph <- within(richness_day7glyph, percentile <- as.integer(.bincode(gene_richness_relative, quantile(gene_richness_relative, probs=0:100/100), include.lowest=TRUE)))
richness_day14glyph <- within(richness_day14glyph, percentile <- as.integer(.bincode(gene_richness_relative, quantile(gene_richness_relative, probs=0:100/100), include.lowest=TRUE)))
richness_day22glyph <- within(richness_day22glyph, percentile <- as.integer(.bincode(gene_richness_relative, quantile(gene_richness_relative, probs=0:100/100), include.lowest=TRUE)))
richness_day43glyph <- within(richness_day43glyph, percentile <- as.integer(.bincode(gene_richness_relative, quantile(gene_richness_relative, probs=0:100/100), include.lowest=TRUE)))
richness_day71glyph <- within(richness_day71glyph, percentile <- as.integer(.bincode(gene_richness_relative, quantile(gene_richness_relative, probs=0:100/100), include.lowest=TRUE)))
richness_day0control <- within(richness_day0control, percentile <- as.integer(.bincode(gene_richness_relative, quantile(gene_richness_relative, probs=0:100/100), include.lowest=TRUE)))
richness_day22control <- within(richness_day22control, percentile <- as.integer(.bincode(gene_richness_relative, quantile(gene_richness_relative, probs=0:100/100), include.lowest=TRUE)))
richness_day71control <- within(richness_day71control, percentile <- as.integer(.bincode(gene_richness_relative, quantile(gene_richness_relative, probs=0:100/100), include.lowest=TRUE)))

# add percentiles for absolute richness

richness_day0glyph <- within(richness_day0glyph, percentile_abs <- as.integer(.bincode(gene_richness, quantile(gene_richness, probs=0:100/100), include.lowest=TRUE)))
richness_day3glyph <- within(richness_day3glyph, percentile_abs <- as.integer(.bincode(gene_richness, quantile(gene_richness, probs=0:100/100), include.lowest=TRUE)))
richness_day7glyph <- within(richness_day7glyph, percentile_abs <- as.integer(.bincode(gene_richness, quantile(gene_richness, probs=0:100/100), include.lowest=TRUE)))
richness_day14glyph <- within(richness_day14glyph, percentile_abs <- as.integer(.bincode(gene_richness, quantile(gene_richness, probs=0:100/100), include.lowest=TRUE)))
richness_day22glyph <- within(richness_day22glyph, percentile_abs <- as.integer(.bincode(gene_richness, quantile(gene_richness, probs=0:100/100), include.lowest=TRUE)))
richness_day43glyph <- within(richness_day43glyph, percentile_abs <- as.integer(.bincode(gene_richness, quantile(gene_richness, probs=0:100/100), include.lowest=TRUE)))
richness_day71glyph <- within(richness_day71glyph, percentile_abs <- as.integer(.bincode(gene_richness, quantile(gene_richness, probs=0:100/100), include.lowest=TRUE)))
richness_day0control <- within(richness_day0control, percentile_abs <- as.integer(.bincode(gene_richness, quantile(gene_richness, probs=0:100/100), include.lowest=TRUE)))
richness_day22control <- within(richness_day22control, percentile_abs <- as.integer(.bincode(gene_richness, quantile(gene_richness, probs=0:100/100), include.lowest=TRUE)))
richness_day71control <- within(richness_day71control, percentile_abs <- as.integer(.bincode(gene_richness, quantile(gene_richness, probs=0:100/100), include.lowest=TRUE)))

# rowbind the objects to one dataframe
richness_percentiles <- rbind(richness_day0glyph, richness_day3glyph, richness_day7glyph, richness_day14glyph, richness_day22glyph, richness_day43glyph, richness_day71glyph, richness_day0control, richness_day22control, richness_day71control)

# join gene and product names
richness_percentiles$gene_prod <- do.call(paste, c(richness_percentiles[c("adjusted_products", "gene")], sep = "_")) 

# subset genes in which we are interested
compare_percentiles <- subset(richness_percentiles, grepl("phn|gyrA|gyrB|purB|rpoC|recG|gap|dnaJ|dnaK|dmlR|ftsY|ftsZ|rplB|polA|sdhA|thiO|arcosine", gene_prod) & !grepl("carbamoylsarcosine", adjusted_products))

# not included yet:
subset(richness_percentiles, grepl("pho|pst|psp", gene))

# remove genes/products with low richness (why aren't they already excluded?)
compare_percentiles <- subset(compare_percentiles, !grepl("phnA|phnT|dnaK2|glyceraldehyde@3@phosphate@dehydrogenase@2", gene_prod))
compare_percentiles <- subset(compare_percentiles, !grepl("glyceraldehyde@3@phosphate@dehydrogenase_gapA", gene_prod))

# put data into wide format
library(reshape2)
relative_percentiles_cast <- dcast(compare_percentiles, gene_prod ~ new_day + treatment, value.var = "percentile")
abs_percentiles_cast <- dcast(compare_percentiles, gene_prod ~ new_day + treatment, value.var = "percentile_abs")
write.csv(relative_percentiles_cast, file = "relative_percentiles_cast.csv")
write.csv(abs_percentiles_cast, file = "abs_percentiles_cast.csv")

# to remove
"glyceraldehyde@3@phosphate@dehydrogenase" & gapA
"glyceraldehyde@3@phosphate@dehydrogenase@2"
"dnaK2"
"phnT"
"phnA"

# housekeeping genes
"gyrA|gyrB|purB|rpoC|recG|gap|dnaJ|dnaK|dmlR|ftsY|ftsZ|rplB|polA|sdhA"


subset(all_products_richness, grepl("monomeric@sarc", adjusted_products))
subset(all_products_richness, grepl("gapA", gene))
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

show genes that correlate or anti correlate e.g. phn control genes in treatment
which genes of a set a more variable and which are more stable?

include pho/pst/psp genes?