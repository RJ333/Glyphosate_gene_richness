

#### problem: prokka output does not contain all products > 1, running samtools again



# what to do with merged output?
--> calculate normalized reads per contig and per product (for tax and gene abundance)
--> GC content biplot length of contigs, phn genes marked
--> script for abundance data
--> script for richness data
important genes and taxonomic annotation, abundance in 16S data


library(data.table)
library(argparse)
library(tidyverse)

# Rscript 01_funct_tax_data_merging.bash -i=/data/Rene/glyph -o=/data/Rene/glyph/merged_input

parser <- ArgumentParser()

parser$add_argument("-i", "--input_directory", default = NULL,
                    dest = "input", required = TRUE,
                    help = "directory containing input file")
					
# parser$add_argument("-t", "--threshold_richness", default = NULL,
                    dest = "threshold", required = TRUE, type= "integer",
                    help = "threshold defining on how many contigs a product \
							has to be found to be included in further analysis")
					
parser$add_argument("-o", "--output_dir", default = NULL,
                    dest = "output", required = TRUE,
                    help = "directory to store merged input data")					
args <- parser$parse_args()

# reading in prokka data

# noch probleme mit 1. spalte (data.table mag keine row,names aus dateien)
# to remove row names (should be done when generating file in first place)
# old <- read.table("./omics/merge_input/prokka_processed.tsv", sep = "\t", header = TRUE, row.names = 1)
# write.table(old, "./omics/merge_input/prokka_processed.tsv", sep = "\t", row.names = FALSE)


print("reading prokka data ...")
# with argparse
prokka <- fread(args$input/"prokka_processed.tsv",
  colClasses = c("factor", "factor", "integer", "factor", "factor", "factor", "factor"))
# or interactive
prokka <- fread("./omics/merge_input/prokka_processed.tsv",
  colClasses = c("factor", "factor", "integer", "factor", "factor", "factor", "factor"))
  
  

  # reading in kaiju data
# we can't set the colClasses directly, as we read info about 15 files and then drop the first
# this is fixed in kaiju script, row names are dropped in future
#old <- read.table("./omics/merge_input/filtered_kaiju_output.tsv", sep = "\t", header = TRUE, row.names = 1)
#write.table(old, "./omics/merge_input/filtered_kaiju_output.tsv", sep = "\t", row.names = FALSE)

print("reading kaiju data ...")
# with argparse
kaiju <- fread(args$input/"filtered_kaiju_output.tsv")
# or interactive
kaiju <- fread("./omics/merge_input/filtered_kaiju_output.tsv",
  colClasses = c("factor", "factor", "integer", "integer", "factor", "factor","factor", 
				 "factor","factor", "double", "double", "double", "double", "double"))

# reading in bbmap data
print("reading bbmap data ...")
# with argparse
bbmap <- fread(args$input/"contig_length_coverage.tsv",
  col.names = c("sample", "contig_id", "average_coverage", "contig_length"),
  colClasses = c("factor", "factor", "double", "integer"))
# or interactive
bbmap <- fread("./omics/merge_input/contig_length_coverage.tsv",
  col.names = c("sample", "contig_id", "average_coverage", "contig_length"),
  colClasses = c("factor", "factor", "double", "integer"))

# reading in reads per product data
print("reading reads per product data ...")
# with argparse
product_reads <- fread(args$input/"contig_gene_sample_reads.tsv", header = TRUE,
  col.names = c("contig_id", "product2", "sample", "reads_per_product"),
  colClasses = c("factor", "factor", "factor", "integer"))
  # or interactive
product_reads <- fread("./omics/merge_input/contig_gene_sample_reads.tsv", header = TRUE,
  col.names = c("contig_id", "product2", "sample", "reads_per_product"),
  colClasses = c("factor", "factor", "factor", "integer"))

# reading in concoct binning data
# print("reading concoct binning data ...")
# concoct <-

# reading in checkm binning quality data
# print("reading checkM binning QC data ...")
# checkm <-

# rechecking character substition before merging
product_reads$product2 <- gsub("_|-|'| |/|:|\\(|\\.|\\)|\\|", "@", product_reads$product2)
prokka$product2 <- 		  gsub("_|-|'| |/|:|\\(|\\.|\\)|\\|", "@", prokka$product2)
product_reads$product2 <- as.factor(product_reads$product2)
prokka$product2 <- as.factor(prokka$product2)
# ' was not replaced to @ in product_reads, resulted in 36673 NAs

# merging data.tables
prokka_bbmap <- merge(prokka, bbmap, 
  by.x = c("sample", "contig_id"), 
  by.y = c("sample", "contig_id"), all.x = TRUE)

prokka_bbmap_prodreads <- merge(prokka_bbmap, product_reads, 
  by.x = c("sample", "contig_id", "product2"), 
  by.y = c("sample", "contig_id", "product2"), all.x = TRUE)

prokka_bbmap_prodreads_kaiju <- merge(prokka_bbmap_prodreads, kaiju,
  by.x = c("sample", "contig_id"), 
  by.y = c("sample", "contig_id"), all.x = TRUE)

# now all data to calculate the tpm (RNA) or rpm (DNA) is present: 
# http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

# first generate rpk for each gene and each kaiju annotated contig
prokka_bbmap_prodreads_kaiju$product_reads_rpk <- prokka_bbmap_prodreads_kaiju$reads_per_product/(prokka_bbmap_prodreads_kaiju$gene_length/1000)
prokka_bbmap_prodreads_kaiju$kaiju_reads_rpk <- prokka_bbmap_prodreads_kaiju$total_reads/(prokka_bbmap_prodreads_kaiju$contig_length.y/1000)


# get scaling factors for each sample

# initialize vector
scaling_factor <- vector("numeric")

scaling_factor_kaiju <- vector("numeric")
# calculate scaling factor
for(i in c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10")) {
 value <- sum(subset(prokka_bbmap_prodreads_kaiju,sample == i)$product_reads_rpk, na.rm = TRUE)
 scaling_factor[i] <- value }
 
# calculate scaling factor for kaiju data
for(i in c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10")) {
 value <- sum(subset(prokka_bbmap_prodreads_kaiju,sample == i)$kaiju_reads_rpk, na.rm = TRUE)
 scaling_factor_kaiju[i] <- value }

 <- as.data.frame(scaling_factor_kaiju)
 
# merge with omics meta data
meta_omics_small <- read.csv("./omics/merge_input/meta_omics_small.csv", header = TRUE, sep = ";")

 # add vector to meta data 
meta_omics_small <- cbind(meta_omics_small, scaling_factor, scaling_factor_kaiju)

prokka_all <- merge(prokka_bbmap_prodreads_kaiju, meta_omics_small, by = "sample")

names(prokka_all)[names(prokka_all) == 'contig_length.x'] <- 'contig_length'
names(prokka_all)[names(prokka_all) == 'total_reads.x'] <- 'total_reads_kaiju'
names(prokka_all)[names(prokka_all) == 'total_reads.y'] <- 'library_size'

# check how sum of reads per kilobase correlate to sequencing depth (in read pairs)
library(ggplot2)
ggplot(meta_data_omics_small, aes( x = new_day, group = treatment, lty = treatment)) +
  geom_line(aes(y = total_reads * 2, colour = "total single reads"), size = 2)+
  geom_line(aes(y = scaling_factor, colour = "scale factor"), size = 2)

# generate reads per million for products and kaiju
prokka_all[,product_rpm := product_reads_rpk/(scaling_factor/1000000), by = sample]
prokka_all[,kaiju_rpm := kaiju_reads_rpk/(scaling_factor_kaiju/1000000), by = sample]

# check if each sample has a total of 1 million reads now for products
sum(subset(prokka_all, sample == "A1")$product_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "A2")$product_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "A3")$product_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "A4")$product_rpm, na.rm = TRUE)  # 1e+06  
sum(subset(prokka_all, sample == "A5")$product_rpm, na.rm = TRUE)  # 1e+06  
sum(subset(prokka_all, sample == "A6")$product_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "A7")$product_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "B8")$product_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "B9")$product_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "B10")$product_rpm, na.rm = TRUE)  # 1e+06
# check if each sample has a total of 1 million reads now for products
sum(subset(prokka_all, sample == "A1")$kaiju_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "A2")$kaiju_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "A3")$kaiju_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "A4")$kaiju_rpm, na.rm = TRUE)  # 1e+06  
sum(subset(prokka_all, sample == "A5")$kaiju_rpm, na.rm = TRUE)  # 1e+06  
sum(subset(prokka_all, sample == "A6")$kaiju_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "A7")$kaiju_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "B8")$kaiju_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "B9")$kaiju_rpm, na.rm = TRUE)  # 1e+06
sum(subset(prokka_all, sample == "B10")$kaiju_rpm, na.rm = TRUE)  # 1e+06

# remove some columns which are probably not important

write.table(prokka_all_tax, file = "prokka_tax_glyph.tsv", sep = "\t")