# data is collected in prokka_all in workspace data_merging_single_assemblies.RData
library(ggplot2)
library(data.table)

# create subset for richness and abundance calculation
all_products <- prokka_all[,c(4,9,11,14,15,19)]

# calculating richness data per gene
all_products_richness <- aggregate( . ~  adjusted_products + gene + new_day + treatment, data = all_products, length)
all_products_richness <- all_products_richness[with(all_products_richness, order(adjusted_products)), ]
all_products_richness <- all_products_richness[,c(1:5)]
colnames(all_products_richness)[colnames(all_products_richness)=="contig_id"] <- "gene_richness"
all_products_richness <- setDT(all_products_richness)
all_products_richness[,gene_richness_relative := gene_richness/gene_richness[new_day == 0]*100, by = .(adjusted_products, treatment)]
all_products_richness[is.na(all_products_richness)] <- 0
  
  
# calculate mean value for the richness data

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

# vectors for later grouping between functional richness and taxonomic richness
functional<- c("group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func")
# taxonomic <- c(rep("group_tax", 32))  # because we have taxonomy from 32 DNA samples, might be already present in tax_rich.csv

# read in taxonomic richness file, add variable for grouping
tax_rich_omics <- read.csv(file.choose(), header = TRUE)  # tax_rich.csv
tax_rich_omics <- tax_rich_omics[,-1]

# creating statistical data based on script summarySE, it does not understand "NA"! dataframe[is.na(dataframe)] <- 0
all_products_richness_summary <- summarySE(all_products_richness, measurevar = "gene_richness_relative", groupvars = c("new_day" ,"treatment")) 
all_products_richness_summary["group"] <- functional

my_y_title <- expression(atop("functional gene and species richness","in relation to day 0 [%]"))
my_legend <- expression("all genes richness")

all_products_richness_mean_ci <- ggplot(all_products_richness_summary, aes(x = new_day, y = gene_richness_relative, group = treatment, colour = treatment))+
geom_line(aes(linetype = group), size = 2)+
geom_ribbon(aes(ymax = gene_richness_relative + ci, ymin = gene_richness_relative - ci), alpha = 0.1, colour = NA)+
geom_line(data = tax_rich_omics, aes(x = days, y = richness_rel, linetype = group), size = 1.5)+
#geom_point(data = more_cell_counts_0, aes(x = new_day, y = cells_rel, group = treatment), size = 3, alpha = 0.8)+
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
	xlab("Days after glyphosate addition")+
	ylab(my_y_title)+
	guides(lty = guide_legend(keywidth = 2, keyheight = 1))+
	coord_cartesian(ylim = c(50, 280)) 
all_products_richness_mean_ci