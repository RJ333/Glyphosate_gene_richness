# Set the working dir with mothur files in it
setwd("/data/projects/glyphosate/reads/mothur_processed/")
load("mothur_glyph_002.RData")

# load libraries
library(phyloseq)
library(data.table)
library(VennDiagram)
library(ggplot2)
library(reshape2)

# define function to plot Venn diagram with 4 categories
fourway.Venn <- function(A,B,C,D,cat.names = c("Water\nDNA",
											   "Biofilm\nDNA",
											   "Water\nRNA",
											   "Biofilm\nRNA")){
  grid.newpage()
  area1 <- length(A)
  area2 <- length(B)
  area3 <- length(C)
  area4 <- length(D)
  n12<-length(Reduce(intersect, list(A,B)))
  n13<-length(Reduce(intersect, list(A,C)))
  n14<-length(Reduce(intersect, list(A,D)))
  n23<-length(Reduce(intersect, list(B,C)))
  n24<-length(Reduce(intersect, list(B,D)))
  n34<-length(Reduce(intersect, list(C,D)))
  n123<-length(Reduce(intersect, list(A,B,C)))
  n124<-length(Reduce(intersect, list(A,B,D)))
  n134<-length(Reduce(intersect, list(A,C,D)))
  n234<-length(Reduce(intersect, list(B,C,D)))
  n1234<-length(Reduce(intersect, list(A,B,C,D)))
  
venn.plot <- draw.quad.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  area4 = area4,
  n12 = n12,
  n13 = n13,
  n14 = n14,
  n23 = n23,
  n24 = n24,
  n34 = n34,
  n123 = n123,
  n124 = n124,
  n134 = n134,
  n234 = n234,
  n1234 = n1234,
  category = cat.names,
  cat.pos = c(0,180,0,200),
  fill = c("blue", "red", "green", "yellow"),
  alpha = .3,
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "green", "black")
)
grid.draw(venn.plot)
}

plot_folder <- "/data/projects/glyphosate/plots/R/community_composition/"

# prepare phyloseq object

# remove OTUs with less than 2 reads in at least 1 sample
mothur_ps3 <- filter_taxa(mothur_ps2, function (x) {sum(x > 1) >= 1}, prune = TRUE)
# transform into relative abundance, displayed in percentage!
mothur_melt <- psmelt(mothur_ps3)
mothur_melt$OTU <- as.factor(mothur_melt$OTU)

# generate subsets for each category (based on treatment or nucleic acid)

# treatment
water_glyph_otus <- subset(mothur_melt, habitat == "water" & treatment == "glyph" & Abundance > 1)
water_glyph_unique_otus<-water_glyph_otus[which(!duplicated(water_glyph_otus[,"OTU"])),]
nrow(water_glyph_unique_otus)
# 4476

water_control_otus <- subset(mothur_melt, habitat == "water" & treatment == "control" & Abundance > 1)
water_control_unique_otus<-water_control_otus[which(!duplicated(water_control_otus[,"OTU"])),]
nrow(water_control_unique_otus)
# 6785

biofilm_glyph_otus <- subset(mothur_melt, habitat == "biofilm" & treatment == "glyph" & Abundance > 1)
biofilm_glyph_unique_otus<-biofilm_glyph_otus[which(!duplicated(biofilm_glyph_otus[,"OTU"])),]
nrow(biofilm_glyph_unique_otus)
# 1684

biofilm_control_otus <- subset(mothur_melt, habitat == "biofilm" & treatment == "control" & Abundance > 1)
biofilm_control_unique_otus<-biofilm_control_otus[which(!duplicated(biofilm_control_otus[,"OTU"])),]
nrow(biofilm_control_unique_otus)
# 1765

# nucleic_acid
water_dna_otus <- subset(mothur_melt, habitat == "water" & nucleic_acid == "dna" & Abundance > 1)
water_dna_unique_otus<-water_dna_otus[which(!duplicated(water_dna_otus[,"OTU"])),]
nrow(water_dna_unique_otus)
# 3912

water_cdna_otus <- subset(mothur_melt, habitat == "water" & nucleic_acid == "cdna" & Abundance > 1)
water_cdna_unique_otus<-water_cdna_otus[which(!duplicated(water_cdna_otus[,"OTU"])),]
nrow(water_cdna_unique_otus)
# 7447

biofilm_dna_otus <- subset(mothur_melt, habitat == "biofilm" & nucleic_acid == "dna" & Abundance > 1)
biofilm_dna_unique_otus<-biofilm_dna_otus[which(!duplicated(biofilm_dna_otus[,"OTU"])),]
nrow(biofilm_dna_unique_otus)
# 1646

biofilm_cdna_otus <- subset(mothur_melt, habitat == "biofilm" & nucleic_acid == "cdna" & Abundance > 1)
biofilm_cdna_unique_otus<-biofilm_cdna_otus[which(!duplicated(biofilm_cdna_otus[,"OTU"])),]
nrow(biofilm_cdna_unique_otus)
# 1842

# plot Venn diagram (adjust labels in function)
fourway.Venn(water_glyph_unique_otus$OTU,
			 biofilm_glyph_unique_otus$OTU,
			 water_control_unique_otus$OTU,
			 biofilm_control_unique_otus$OTU)
			 
dev.copy(png, paste(plot_folder, "4wayVenn_treatments.png"))
dev.off()

fourway.Venn(water_dna_unique_otus$OTU,
			 biofilm_dna_unique_otus$OTU,
			 water_cdna_unique_otus$OTU,
			 biofilm_cdna_unique_otus$OTU)
			 
dev.copy(png, paste(plot_folder, "4wayVenn_nucleic_acids.png"))
dev.off()


# only biofilm vs water
water_otus <- subset(mothur_melt, habitat == "water" & Abundance > 1)
water_unique_otus<-water_otus[which(!duplicated(water_otus[,"OTU"])),]
nrow(water_unique_otus)
# 10692

biofilm_otus <- subset(mothur_melt, habitat == "biofilm" & Abundance > 1)
biofilm_unique_otus<-biofilm_otus[which(!duplicated(biofilm_otus[,"OTU"])),]
nrow(biofilm_unique_otus)
# 2903

length(intersect(water_unique_otus$OTU, biofilm_unique_otus$OTU))
# 743

# plot pairwise venn diagram
grid.newpage()
venn.plot <- draw.pairwise.venn(area1        = 10692,
                                area2        = 2903,
                                cross.area   = 743,
                                scaled       = T,
                                category     = c("Water", "Biofilm"),
                                fill         = c("blue", "red"),
                                alpha        = 0.3,
                                lty          = "blank",
                                cex          = 2,
                                cat.cex      = 2,
                                cat.pos      = c(285, 105),
                                cat.dist     = 0.09,
                                cat.just     = list(c(-1, -1), c(1, 1)),
                                ext.pos      = 30,
                                ext.dist     = -0.05,
                                ext.length   = 0.85,
                                ext.line.lwd = 2,
                                ext.line.lty = "dashed")
grid.draw(venn.plot)
dev.copy(png, paste(plot_folder, "Venn_water_biofilm.png"))
dev.off()

