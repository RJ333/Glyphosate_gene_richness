#!/usr/bin/env Rscript

##  TODO:

## compare with 04_get_reads_per... what is now überflüssig and can be removed from there

## when is table ready for data merging? make sure same characters were substitued by "@"
## in available data and change for future runs


library(argparse)
library(data.table)
library(tidyverse)
library(splitstackshape)

parser <- ArgumentParser()

parser$add_argument("-p", "--prokka_file", default = NULL,
                    dest = "prokka_file", required = TRUE,
                    help = "file containing modified and combined prokka files")
parser$add_argument("-t", "--threshold_richness", default = 1,
                    dest = "threshold", required = TRUE, type= "integer",
                    help = "threshold defining on how many contigs a product \
							has to be found to be included in further analysis")
					
parser$add_argument("-o", "--output_dir", default = NULL,
                    dest = "output_dir", required = TRUE,
                    help = "directory to store processed prokka data")					
args <- parser$parse_args()

# call like this: 
# conda activate Renv
# time Rscript 02_modify_prokka_products.r \
# -p="/data/Rene/glyph/prokka/prokka_all_modified.tsv" -o="/data/Rene/glyph/prokka/"


# read in contig_gene_sample_reads 
contig_gene_sample_reads <- fread("./omics/prokka/contig_gene_sample_reads.tsv")
contig_gene_sample_reads$gene <- gsub("_|-|'| |/|:|\\.|\\(|\\)|\\|", "@", contig_gene_sample_reads$gene)  


# read in the combined prokka data from all samples 
print("reading prokka data ...")
#prokka_tables <- fread(args$prokka_file, 
prokka_tables <- fread("./omics/prokka/prokka_all_modified.tsv", 
  col.names = c("sample", "contig_id", "gene_start", "gene_end", "annotation"),
  colClasses = c("factor", "factor", "integer", "integer", "character"))

# splitting annotation-column
prokka_neat <- prokka_tables %>% 
  cSplit("annotation", sep = ";") %>%
  gather(Key, value, -c(sample, contig_id, gene_start, gene_end)) %>% 
  separate(value, c("Col", "Value"), sep = "=") %>% 
  select(-Key) %>% 
  filter(!(is.na(Col) & is.na(Value))) %>% 
  spread(Col, Value)

# keep this during devel, if something fails 
full <- prokka_neat  
# prokka_neat <- full


# combine "note" and "product"
print("combining columns product and note")
prokka_neat$product2 <- ifelse(prokka_neat$product == "hypothetical protein" & !is.na(prokka_neat$note), prokka_neat$note, prokka_neat$product)

# calculate gene length
print("calculating gene length")
prokka_neat$gene_length <- prokka_neat$gene_end - prokka_neat$gene_start

# remove unneeded columns "note" and "product"
# gene start and end is needed for bed file generation!
print("selecting important columns")
prokka_select <- prokka_neat[, c("sample", "contig_id", "gene_start", "gene_end", "gene_length",
								 "eC_number", "gene", "product2", "locus_tag")]

## unify product names , all to lower
print("unifying product annotations")
prokka_select$product2 <- tolower(prokka_select$product2)

## replace white lines and special characters in product names with "@"
prokka_select$product2 <- gsub("_|-|'| |/|:|\\(|\\.|\\)|\\|", "@", prokka_select$product2)  


# unify gene numbering e.g. genE_01, genE_02 etc --> genE
print("unifying gene names")
prokka_select$gene <- gsub("_.*$","", prokka_select$gene)

# turning columns into factors
prokka_select$eC_number <- as.factor(prokka_select$eC_number)
prokka_select$product2 <- as.factor(prokka_select$product2)
prokka_select$gene <- as.factor(prokka_select$gene)
prokka_select$locus_tag <- as.factor(prokka_select$locus_tag)

# create table with all products across all samples and their global richness
all_unique_products_with_count <- as.data.frame(table(droplevels(prokka_select$product2)))

# generate vector of unique products with defined threshold
unique_products_selected <- all_unique_products_with_count %>%
filter(Freq > args$threshold) %>%
select(Var1)

unique_products_selected <- all_unique_products_with_count %>%
filter(Freq > 1) %>%
select(Var1)


# unique products for samtools script
print("writing unique products across all samples ...")
write.table(unique_products_selected , sep = "\t", 
  file = paste0(args$output_dir, "/unified_unique_prokka_products_greater_",args$threshold,".tsv"))
print("all done")

# unique products for samtools script
print("writing unique products across all samples ...")
write.table(unique_products_selected , sep = "\t", 
  file = paste0("./omics/prokka", "/unified_unique_prokka_products_greater_",1,".tsv"))
print("all done")


# generate prokka output for samtools to create bed.files
prokka_bed  <- setDT(prokka_select)
prokka_bed[, fwrite(.SD, sep = "\t", paste0(args$output_dir, "/prokka_for_bed_", sample,".tsv")), by = sample]

prokka_bed[, fwrite(.SD, sep = "\t", paste0("./omics/prokka", "/prokka_for_bed_", sample,".tsv")), by = sample]

# remove gene start/end from this file (column 3 and 4)
prokka_processed <- prokka_select[, c("sample", "contig_id", "gene_length",
								 "eC_number", "gene", "product2", "locus_tag")]

# write prokka files for further merging
print("writing processed prokka file ...")
write.table(prokka_processed , sep = "\t", file = paste0(args$output_dir, "/prokka_processed.tsv"))
write.table(prokka_processed , sep = "\t", file = paste0("./omics/prokka", "/prokka_processed.tsv"))
print("all done")

sessionInfo()