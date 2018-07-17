# this script uses the reads per gene data, which has been processed before
# it generates richness and abundance data from it and,also relation to day 0
# and plots it. 
# gene richness should be at least above 5 once to be taken into account
# mapped reads are not normalized by length of gene nor sequencing depth yet 

# workspace is gene_richness_single.RData

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
product_reads2 <- prokka_all[,c(1:3,6,25,26,31)]  # contigs_genes_sample_reads.tsv

# vectors for later grouping between functional richness and taxonomic richness
fuck <- c("group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func")
fuck2 <- c(rep("group_tax", 32))

# read in taxonomic richness file, add variable for grouping
tax_rich <- read.delim(file.choose(), header = TRUE, row.names = 1, sep = ",")
# tax_rich["group"] <- fuck2 
tax_rich_omics <- subset(tax_rich, new_day >= 0)



# calculating richness data per gene
# missing genes are a problem for aggregrate
product_reads2$gene <- as.character(product_reads2$gene)
product_reads2[is.na(product_reads2)] <- "hello"
product_reads2$gene <- as.factor(product_reads2$gene)

product_reads2$product_rpm <- as.numeric(product_reads2$product_rpm)

gene_richness <- aggregate( . ~  product2 + gene + sample + new_day + treatment, data = product_reads2, length)

gene_richness <- gene_richness[with(gene_richness, order(product2)), ]
gene_richness <- gene_richness[,c(1:6)]
colnames(gene_richness)[colnames(gene_richness)=="contig_id"] <- "gene_richness"

abu_sum <- aggregate( . ~  product2 + gene + sample + new_day + treatment, data = product_reads2, sum)
abu_sum <- abu_sum[with(abu_sum, order(product2)), ]
abu_sum <- abu_sum[, c(1:5,7)]

gene_richness2 <- merge(gene_richness, abu_sum, by = c("product2", "gene", "sample", "new_day", "treatment"))
# add column with relative richness values with data.table
gene_richness2 <- setDT(gene_richness2)
gene_richness2[,gene_richness_relative := gene_richness/gene_richness[new_day == 0]*100, by = .(product2, gene, treatment)]


  
  
  
# removing genes with absolute richness not being above > 5 at least once in treatment
gene_richness2 <- droplevels(subset(gene_richness2, !grepl("fdm|ydiF|norB|malF|ktrB|pphA|phnR|phnA|phnT|chtA|chiA|xynB|bcsZ|celB|cenC", gene)))

# ydiF, malF, ktrB from "phnM_like" too small
# norB from "nitrogen" too small
# pphA, phnR, phnA, phnT from "more_phosphonate"
# chtA not found at all, chiA always 2
# all glycosyl hydrolases too snall

new_genes <- subset(gene_richness2, grepl("gap|Sarc|dnaJ|dnaK|dmlR|ftsY|ftsZ|rplB|polA|sdhA", gene))

not_phn_rich <- subset(gene_richness2, !grepl("phn[C-N]|soxA|soxB|Sarc|gyrA|gyrB|purB|rpoC|recG|gap|dnaJ|dnaK|dmlR|ftsY|ftsZ|rplB|polA|sdhA", gene))

sox_richness_subset <- subset(gene_richness2, grepl("sox[AB]|Sarc", gene))

deg_richness <- subset(gene_richness2, grepl("phn[C-N]|soxA|soxB|Sarc", gene))

phn_operon_richness <- subset(gene_richness2, grepl("phn[C-N]", gene))

housekeeping_richness <- subset(gene_richness2, grepl("gyrA|gyrB|purB|rpoC|recG|gap|dnaJ|dnaK|dmlR|ftsY|ftsZ|rplB|polA|sdhA", gene))

p_starve_rich <- subset(gene_richness2, grepl("pho|pst|psp", gene))





phnM_like <- subset(gene_richness2, grepl("phnM|nuoF|mntB|araD|tfdB|artI|arfA|qedA|mlhB|purB", gene))
phn_operon_richness_subset_glyph <- subset(gene_richness2, grepl("phn[C-N]", gene) & treatment == "glyph")
ref_richness_subset_glyph <- subset(gene_richness2, !grepl("phn[C-P]|sox[A-B]", gene) & treatment == "glyph")
nitrogen_rich_subset <- subset(gene_richness2, grepl("nirS|nosZ", gene))
katalase_rich <- subset(gene_richness2, grepl("kat", gene))
more_phosphonate <- subset(gene_richness2, grepl("phn[USWX]|thiO", gene))
metal <- subset(gene_richness2, grepl("merA|czcD", gene))  # differing in control, very similar in treatment, also to phnE
monooxy <- subset(gene_richness2, grepl("tmoS", gene))  # could be useful



  
# remove to small genes, replace NA with 0
housekeeping_richness[is.na(housekeeping_richness)] <- 0
deg_richness[is.na(deg_richness)] <- 0
not_phn_rich[is.na(not_phn_rich)] <- 0
p_starve_rich[is.na(p_starve_rich)] <- 0
phn_operon_richness[is.na(phn_operon_richness)] <- 0

# creating statistical data based on script summarySE, it does not understand "NA"! dataframe[is.na(dataframe)] <- 0
housekeeping_richness_summary <- summarySE(housekeeping_richness, measurevar = "gene_richness2_relative", groupvars = c("new_day" ,"treatment")) 
housekeeping_richness_summary["group"] <- fuck 

deg_richness_summary <- summarySE(deg_richness, measurevar = "gene_richness2_relative", groupvars = c("new_day" ,"treatment")) 
deg_richness_summary["group"] <- fuck 
  
not_phn_rich_summary <- summarySE(not_phn_rich, measurevar = "gene_richness2_relative", groupvars = c("new_day" ,"treatment")) 
not_phn_rich_summary["group"] <- fuck 

p_starve_rich_summary <- summarySE(p_starve_rich, measurevar = "gene_richness2_relative", groupvars = c("new_day" ,"treatment")) 
p_starve_rich_summary["group"] <- fuck 

phn_operon_richness_summary <- summarySE(phn_operon_richness, measurevar = "gene_richness2_relative", groupvars = c("new_day" ,"treatment")) 
phn_operon_richness_summary["group"] <- fuck 

# plot the averaged subset with confidence interval

#my_y_title <- expression(paste("not ",italic("phn"), " operon and species richness in relation to day 0 [%]"))
#my_legend <- expression(paste("not ",italic("phn"), " operon richness"))

my_y_title <- expression(atop("functional gene and species richness","in relation to day 0 [%]"))
my_legend <- expression(paste(italic("phn"),"-operon gene richness"))

phn_operon_richness_mean_ci <- ggplot(phn_operon_richness_summary, aes(x = new_day, y = gene_richness2_relative, group = treatment, colour = treatment))+
geom_line(aes(linetype = group), size = 2)+
geom_ribbon(aes(ymax = gene_richness2_relative + ci, ymin = gene_richness2_relative - ci), alpha = 0.1, colour = NA)+
geom_line(data = tax_rich_omics, aes(y = richness_rel, linetype = group), size = 1.5)+
#geom_point(data = more_cell_counts_0, aes(x = new_new_day, y = cells_rel, group = treatment), size = 3, alpha = 0.8)+
scale_colour_manual(values = c("glyph" = "black", "control" = "grey70"),
						name = "  ",
						breaks = c("glyph", "control"),
						labels = c("Treatment", "Control"))+
scale_linetype_manual(values = c("group_func" = 1, "group_tax" = 3),
						name = "",
						breaks = c("group_func", "group_tax"),
						labels = c(my_legend, "16S rRNA gene richness"))+								
theme_bw()+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	theme(axis.title = element_text(size = 25))+
	theme(axis.text.x = element_text(angle = 0))+
	theme(axis.title.y = element_text(angle = 90, hjust = 0.5))+
	theme(axis.text=element_text(size = 20))+
	theme(legend.position="right")+
	theme(strip.background = element_blank(),
		strip.text.x = element_blank())+
	xlab("new_day after glyphosate addition")+
	ylab(my_y_title)+
	guides(lty = guide_legend(keywidth = 2, keyheight = 1))+
	coord_cartesian(ylim = c(50, 280)) 
phn_operon_richness_mean_ci
ggsave(file = "phn_operon_richness_and_tax.png", width = 12, height = 10)
ggsave(file = "phn_operon_richness_and_tax.pdf", width = 12, height = 10)
ggsave(file = "phn_operon_richness_and_tax.eps", width = 12, height = 10)
  
  
########

ggplot(monooxy, aes(x = new_day, y = gene_richness2, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(katalase_rich, aes(x = new_day, y = gene_richness2, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(more_phosphonate, aes(x = new_day, y = gene_richness2, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)
  
ggplot(p_starve_rich_subset, aes(x = new_day, y = gene_richness2, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(phnM_like, aes(x = new_day, y = gene_richness2, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)  
  
ggplot(housekeeping_richness, aes(x = new_day, y = gene_richness2_relative, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(nitrogen_rich_subset, aes(x = new_day, y = gene_richness2, group = gene, colour = gene)) +
  geom_line(size = 2) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)

ggplot(phn_operon_richness_subset_glyph, aes(x = new_day, y = gene_richness2_relative, group = gene, colour = gene)) +
  geom_line(size = 1.0) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3) +
  facet_grid(~ treatment)

ggplot(sox_richness_subset, aes(x = new_day, y = gene_richness2, group = gene, colour = gene)) +
  geom_line(size = 2) +
  geom_text(aes(label = gene), show.legend = FALSE) +
  facet_grid(~ treatment)
  
ggplot(ref_richness_subset_glyph, aes(x = new_day, y = gene_richness2, group = gene, colour = gene)) +
  geom_line(size = 1.2) +
  geom_text(aes(label = gene), show.legend = FALSE, size = 3.5) +
  geom_line(data = phn_operon_richness_subset_glyph, aes(x = new_day, y = gene_richness2_relative, group = gene),size = 0.8)
  geom_text(data = phn_operon_richness_subset_glyph, label = gene, show.legend = FALSE, size = 3.5)







