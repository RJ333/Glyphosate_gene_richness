# this script uses the reads per gene data, which has been processed before
# it generates richness and abundance data from it and,also relation to day 0
# and plots it. 
# gene richness should be at least above 5 once to be taken into account
# mapped reads are not normalized by length of gene nor sequencing depth yet 

library(ggplot2)
library(data.table)
library(scales)

# define function summarySE

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean   (xx[[col]], na.rm=na.rm),
          sd   = sd     (xx[[col]], na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

# read in mapping file
gene_reads <- read.delim(file.choose(), header = TRUE)  # contigs_genes_sample_reads.tsv

# vectors for later grouping between functional richness and taxonomic richness
fuck <- c("group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func")
fuck2 <- c(rep("group_tax", 32))

# read in taxonomic richness file, add variable for grouping
tax_rich <- read.delim(file.choose(), header = TRUE)
tax_rich["group"] <- fuck2 
tax_rich_omics <- subset(tax_rich, days >= 0)

# plot taxonomic richness

ggplot(tax_rich_omics, aes(x = days, colour = treatment))+
  geom_line(aes(y = richness), linetype = 1) +
  geom_line(aes(y = richness_rel), linetype = 2) +
  facet_wrap(~ treatment, nrow = 2)


# generate meta data and merge with mapping file
sample <- c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10")
days <- c(0, 3, 7, 14, 22, 43, 71, 0, 22, 71)
treatment <- c("glyph", "glyph", "glyph", "glyph", "glyph", "glyph", "glyph", "control", "control", "control")
sample_meta_data <- data.frame(sample, days, treatment)

gene_reads2 <- merge(gene_reads, sample_meta_data, by = "sample")

# calculating richness data per gene
gene_richness <- aggregate( . ~  gene + sample + days + treatment, data = gene_reads2, length)
gene_richness <- gene_richness[with(gene_richness, order(gene)), ]
gene_richness <- gene_richness[,c(1:5)]
colnames(gene_richness)[colnames(gene_richness)=="contig"] <- "gene_richness"

# add column with relative richness values with data.table
gene_richness <- setDT(gene_richness)
gene_richness[,gene_richness_relative := gene_richness/gene_richness[days == 0]*100, by = .(gene,treatment)]

# removing genes with absolute richness not being above > 5 at least once in treatment
gene_richness <- droplevels(subset(gene_richness, !grepl("fdm|ydiF|norB|malF|ktrB|pphA|phnR|phnA|phnT|chtA|chiA|xynB|bcsZ|celB|cenC", gene)))

# ydiF, malF, ktrB from "phnM_like" too small
# norB from "nitrogen" too small
# pphA, phnR, phnA, phnT from "more_phosphonate"
# chtA not found at all, chiA always 2
# all glycosyl hydrolases too snall

new_genes <- subset(gene_richness, grepl("gap|Sarc|dnaJ|dnaK|dmlR|ftsY|ftsZ|rplB|polA|sdhA", gene))

not_phn_rich <- subset(gene_richness, !grepl("phn[C-N]|soxA|soxB|Sarc|gyrA|gyrB|purB|rpoC|recG|gap|dnaJ|dnaK|dmlR|ftsY|ftsZ|rplB|polA|sdhA", gene))

sox_richness_subset <- subset(gene_richness, grepl("sox[AB]|Sarc", gene))

deg_richness <- subset(gene_richness, grepl("phn[C-N]|soxA|soxB|Sarc", gene))

housekeeping_richness <- subset(gene_richness, grepl("gyrA|gyrB|purB|rpoC|recG|gap|dnaJ|dnaK|dmlR|ftsY|ftsZ|rplB|polA|sdhA", gene))

p_starve_rich <- subset(gene_richness, grepl("pho|pst|psp", gene))


phnM_like <- subset(gene_richness, grepl("phnM|nuoF|mntB|araD|tfdB|artI|arfA|qedA|mlhB|purB", gene))
phn_operon_richness_subset_glyph <- subset(gene_richness, grepl("phn[C-N]", gene) & treatment == "glyph")
ref_richness_subset_glyph <- subset(gene_richness, !grepl("phn[C-P]|sox[A-B]", gene) & treatment == "glyph")
nitrogen_rich_subset <- subset(gene_richness, grepl("nirS|nosZ", gene))
katalase_rich <- subset(gene_richness, grepl("kat", gene))
more_phosphonate <- subset(gene_richness, grepl("phn[USWX]|thiO", gene))
metal <- subset(gene_richness, grepl("merA|czcD", gene))  # differing in control, very similar in treatment, also to phnE
monooxy <- subset(gene_richness, grepl("tmoS", gene))  # could be useful

###### plots for quick investigation: relative richness containing tax richness
  
ggplot(phn_operon_sox_richness_subset, aes(x = days, y = gene_richness)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  #geom_line(data = tax_rich_omics, aes(x = days, y = richness_rel))+
  #geom_text(data = tax_rich_omics, aes(x = days, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
ggplot(sox_richness_subset, aes(x = days, y = gene_richness_relative)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  #geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  #geom_line(data = tax_rich_omics, aes(x = days, y = richness_rel))+
  #geom_text(data = tax_rich_omics, aes(x = days, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
ggplot(housekeeping_richness, aes(x = days, y = gene_richness_relative)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  #geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  geom_line(data = tax_rich_omics, aes(x = days, y = richness_rel), size = 1.3)+
  #geom_text(data = tax_rich_omics, aes(x = days, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(new_genes, aes(x = days, y = gene_richness)) +
  geom_line(aes(group = gene, colour = gene), size = 1.2) +
  #geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  geom_line(data = tax_rich_omics, aes(x = days, y = richness_rel))+
  #geom_text(data = tax_rich_omics, aes(x = days, y = richness_rel, label = "richness"), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

  
# remove to small genes, replace NA with 0
housekeeping_richness[is.na(housekeeping_richness)] <- 0
deg_richness[is.na(deg_richness)] <- 0
not_phn_rich[is.na(not_phn_rich)] <- 0
p_starve_rich[is.na(p_starve_rich)] <- 0

# creating statistical data based on script summarySE, it does not understand "NA"! dataframe[is.na(dataframe)] <- 0
housekeeping_richness_summary <- summarySE(housekeeping_richness, measurevar = "gene_richness_relative", groupvars = c("days" ,"treatment")) 
housekeeping_richness_summary["group"] <- fuck 

deg_richness_summary <- summarySE(deg_richness, measurevar = "gene_richness_relative", groupvars = c("days" ,"treatment")) 
deg_richness_summary["group"] <- fuck 
  
not_phn_rich_summary <- summarySE(not_phn_rich, measurevar = "gene_richness_relative", groupvars = c("days" ,"treatment")) 
not_phn_rich_summary["group"] <- fuck 

p_starve_rich_summary <- summarySE(p_starve_rich, measurevar = "gene_richness_relative", groupvars = c("days" ,"treatment")) 
p_starve_rich_summary["group"] <- fuck 
# plot the averaged subset with confidence interval

#my_y_title <- expression(paste("not ",italic("phn"), " operon and species richness in relation to day 0 [%]"))
#my_legend <- expression(paste("not ",italic("phn"), " operon richness"))

my_y_title <- expression(paste("P-sensing genes and species richness in relation to day 0 [%]"))
my_legend <- expression(paste("P-sensing genes richness"))

p_starve_rich_mean_ci <- ggplot(p_starve_rich_summary, aes(x = days, y = gene_richness_relative, group = treatment, colour = treatment))+
geom_line(aes(linetype = group), size = 2)+
geom_ribbon(aes(ymax = gene_richness_relative + ci, ymin = gene_richness_relative - ci), alpha = 0.1, colour = NA)+
geom_line(data = tax_rich_omics, aes(y = richness_rel, linetype = group), size = 1.5)+
#geom_point(data = more_cell_counts_0, aes(x = new_days, y = cells_rel, group = treatment), size = 3, alpha = 0.8)+
scale_colour_manual(values = c("glyph" = "black", "control" = "grey50"),
						name = "  ",
						breaks = c("glyph", "control"),
						labels = c("Treatment", "Control"))+
scale_linetype_manual(values = c("group_func" = 1, "group_tax" = 3),
						name = "",
						breaks = c("group_func", "group_tax"),
						labels = c(my_legend, "16S rRNA genes richness"))+								
theme_bw()+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	theme(axis.title = element_text(size=20))+
	theme(axis.text.x = element_text(angle=0,vjust=0.5))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=17))+
	theme(legend.position="right")+
	theme(strip.background = element_blank(),
		strip.text.x = element_blank())+
	xlab("Days after glyphosate addition")+
	ylab(my_y_title)+
	guides(lty = guide_legend(keywidth = 1.5, keyheight = 1))+
	coord_cartesian(ylim = c(50, 260)) 
p_starve_rich_mean_ci
ggsave(file = "p_starve_rich_and_tax.png", width = 14, height = 12)
  
  
########

ggplot(monooxy, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(katalase_rich, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(more_phosphonate, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
ggplot(p_starve_rich_subset, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(phnM_like, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)  
  
ggplot(housekeeping_richness, aes(x = days, y = gene_richness_relative, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(nitrogen_rich_subset, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 2) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)

ggplot(phn_operon_richness_subset_glyph, aes(x = days, y = gene_richness_relative, group = gene, colour = gene)) +
  geom_line(size = 1.0) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(sox_richness_subset, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 2) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)
  
ggplot(ref_richness_subset_glyph, aes(x = days, y = gene_richness, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3.5) +
  geom_line(data = phn_operon_richness_subset_glyph, aes(x = days, y = gene_richness_relative, group = gene),size = 0.8)
  geom_text(data = phn_operon_richness_subset_glyph, label = gene, show.legend = FALSE, size = 3.5)







