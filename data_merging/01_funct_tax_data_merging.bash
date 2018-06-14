# list prokka files with sample name

# what files do we need?

# taxonomy/processing_kaiju/quality_check_kaiju_output.sh
kaiju output : taxonomic output per contig : ready

# Glyphosate_gene_richness/data_merging/00_generate_list_contig_length_coverage_bbtools.sh
bbmap output : reads per contig, average coverage and contig length : ready

prokka output : prokka table with adjusted product names, fitting to samtools used names : in prep 

samtools gene richness output : reads per product per sample : ready

concoct binning output?

meta_data for samples --> ready

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
old <- read.table("./omics/merge_input/prokka_processed.tsv", sep = "\t", header = TRUE, row.names = 1)
write.table(old, "./omics/merge_input/prokka_processed.tsv", sep = "\t", row.names = FALSE)


print("reading prokka data ...")
prokka <- fread(args$input/"prokka_processed.tsv",
  colClasses = c("factor", "factor", "integer", "factor", "factor", "factor", "factor"))
prokka <- fread("./omics/merge_input/prokka_processed.tsv",
  colClasses = c("factor", "factor", "integer", "factor", "factor", "factor", "factor"))
  
  

  # reading in kaiju data
# we can't set the colClasses directly, as we read info about 15 files and then drop the first
# this is fixed in kaiju script, row names are dropped in future
old <- read.table("./omics/merge_input/filtered_kaiju_output.tsv", sep = "\t", header = TRUE, row.names = 1)
write.table(old, "./omics/merge_input/filtered_kaiju_output.tsv", sep = "\t", row.names = FALSE)

print("reading kaiju data ...")
kaiju <- fread(args$input/"filtered_kaiju_output.tsv")
kaiju <- fread("./omics/merge_input/filtered_kaiju_output.tsv",
  colClasses = c("factor", "factor", "integer", "integer", "factor", "factor","factor", 
				 "factor","factor", "double", "double", "double", "double", "double"))

# reading in bbmap data
print("reading bbmap data ...")
bbmap <- fread(args$input/"contig_length_coverage.tsv",
  col.names = c("sample", "contig_id", "average_coverage", "contig_length"),
  colClasses = c("factor", "factor", "double", "integer"))
bbmap <- fread("./omics/merge_input/contig_length_coverage.tsv",
  col.names = c("sample", "contig_id", "average_coverage", "contig_length"),
  colClasses = c("factor", "factor", "double", "integer"))

# reading in reads per product data
print("reading reads per product data ...")
product_reads <- fread(args$input/"contig_gene_sample_reads.tsv", header = TRUE)
  col.names = c("contig_id", "product2", "sample", "reads_per_product"),
  colClasses = c("factor", "factor", "factor", "integer"))
product_reads <- fread("./omics/merge_input/contig_gene_sample_reads.tsv", header = TRUE,
  col.names = c("contig_id", "product2", "sample", "reads_per_product"),
  colClasses = c("factor", "factor", "factor", "integer"))

# reading in concoct binning data
# print("reading concoct binning data ...")
# concoct <-

# reading in checkm binning quality data
# print("reading checkM binning QC data ...")
# checkm <-

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

# reorder columns


# add data from read mappings on genes/products

# merge with omics meta data
meta_data_omics_small <- read.csv("./omics/merge_input/meta_omics_small.csv", header = TRUE, sep = ";")
prokka_all <- merge(prokka_all, meta_omics_small, by = "sample")

# now all data to calculate the tpm (RNA) or rpm (DNA) is present: 
# http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

prokka_bbmap_prodreads_kaiju$product_reads_rpk <- prokka_bbmap_prodreads_kaiju$reads_per_product/(prokka_bbmap_prodreads_kaiju$gene_length/1000)
prokka_bbmap_prodreads_kaiju$kaiju_reads_rpk <- prokka_bbmap_prodreads_kaiju$total_reads/(prokka_bbmap_prodreads_kaiju$contig_length.y/1000)

# get scaling factors for each sample

# initialize vector
scaling_factor <- vector("numeric")

# calculate scaling factor
for(i in c("A1", "A2", "A3", "A4", "A5", "A6", "A7", "B8", "B9", "B10")) {
 value <- sum(subset(prokka_bbmap_prodreads_kaiju,sample == i)$product_reads_rpk, na.rm = TRUE)
 scaling_factor[i] <- value }

 # add vector to meta data 
meta_data_omics_small <- cbind(meta_data_omics_small, scaling_factor)

# check how sum of reads per kilobase correlate to sequencing depth (in read pairs)
library(ggplot2)
ggplot(meta_data_omics_small, aes( x = new_day, group = treatment, lty = treatment)) +
  geom_line(aes(y = total_reads * 2, colour = "total single reads"), size = 2)+
  geom_line(aes(y = scaling_factor, colour = "scale factor"), size = 2)

# generate reads per million  

# without merging both data.tables?

prokka_all <- merge(prokka_all, meta_omics_small[ , c("sample", "scale_factor")], by = "sample")
prokka_all$rpm <- prokka_all$rpk/(prokka_all$scale_factor/1000000)
 
# check if each sample has a total of 1 million reads now
sum(subset(prokka_all,sample == "A1")$rpm)  # 1e+06
sum(subset(prokka_all,sample == "A2")$rpm)  # 1e+06
sum(subset(prokka_all,sample == "A3")$rpm)  # 1e+06
sum(subset(prokka_all,sample == "A4")$rpm)  # 1e+06  
sum(subset(prokka_all,sample == "A5")$rpm)  # 1e+06  
sum(subset(prokka_all,sample == "A6")$rpm)  # 1e+06
sum(subset(prokka_all,sample == "A7")$rpm)  # 1e+06
sum(subset(prokka_all,sample == "B8")$rpm)  # 1e+06
sum(subset(prokka_all,sample == "B9")$rpm)  # 1e+06
sum(subset(prokka_all,sample == "B10")$rpm)  # 1e+06


# now only tax information and read normalization for that data is missing
kaiju_50_20 <- read.delim(file.choose(), header = TRUE)
# remove some columns which are probably not important
# relative values for order, class and phylum removed
kaiju_shorter <- kaiju_50_20[,c(1:11)]
# calculation columns for rpm and other unimportant columns removed
prokka_all_shorter <- prokka_all[,c(1,4,7,8,9,11,13,19)]

# merge shortened data
prokka_all_tax <- merge(prokka_all_shorter, kaiju_shorter, by = c("sample", "contig_id"), all.x = TRUE)

write.table(prokka_all_tax, file = "prokka_tax_glyph.tsv", sep = "\t")