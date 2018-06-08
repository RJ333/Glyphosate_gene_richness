#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(tidyverse)
library(splitstackshape)

parser <- ArgumentParser()

parser$add_argument("-p", "--prokka_file", default = NULL,
                    dest = "prokka_file", required = TRUE,
                    help = "file containing modified and combined prokka files")
parser$add_argument("-o", "--output_dir", default = NULL,
                    dest = "output_dir", required = TRUE,
                    help = "directory to store processed prokka data")					
args <- parser$parse_args()

# call like this: 
# conda activate Renv
# time Rscript 02_modify_prokka_products.r \
# -p="/data/Rene/glyph/prokka/prokka_all_modified.tsv" -o="/data/Rene/glyph/prokka/"

# read in the combined prokka data from all samples 
print("reading prokka data ...")
#prokka_tables <- fread(args$prokka_file, 
prokka_tables <- fread("./omics/prokka/prokka_all_modified.tsv", 
  col.names = c("sample", "contig_id", "gene_start", "gene_end", "annotation"),
  colClasses = c("factor", "factor", "integer", "integer", "character"))

#splitting annotation-column
prokka_neat <- prokka_tables %>% 
  cSplit("annotation", sep = ";") %>%
  gather(Key, value, -c(sample, contig_id, gene_start, gene_end)) %>% 
  separate(value, c("Col", "Value"), sep = "=") %>% 
  select(-Key) %>% 
  filter(!(is.na(Col) & is.na(Value))) %>% 
  spread(Col, Value)

  
  
## combine "note" and "product"

full <- prokka_neat
full$product2 <- ifelse(full$product == "hypothetical protein" & !is.na(full$note), full$note, full$product)







  
# generate output
# prokka files for further merging
print("writing processed prokka table as xxx ...")
write.table(xxx , sep = "\t", file = paste0(args$output_dir, "/xxx_prok"))
print("all done")
# unique products for samtools script
print("writing unique products across all samples as xxx ...")
write.table(xxx , sep = "\t", file = paste0(args$output_dir, "/xxx_uniq"))
print("all done")

sessionInfo()
  
##  TODO:

## Test so far, if successful, then:

## name new columns correctly



## calculate gene length

## unify product names , all to lower

## replace white lines and special characters in product names with "@"
#  "/", "\", "-", "_","'", """ " ( ) [ ] 

## compare with 01_get_reads_per... what is now überflüssig and can be removed from there

## turn product names into factors?

## extract unique products for samtools script

## when is table ready for data merging?

## unify gene numbering e.g. genE_01, genE_02 etc --> genE
# in excel
#add headers to the columns: sample contig_id start_pos end_pos gene_length eC_number gene product 
#(also including note=...) products_adjusted
#add column sample_contig_id to uniquely identify contigs from different samples
#add sample_contig_id_product for same reason

#delete columns with prodigal, CDS, "0", "+/-", inference, locus tag, PROKKA id
#add column gene length
#remove counts from genE_01, genE_02 etc --> genE
#replace in product column > adjusted_products: 

