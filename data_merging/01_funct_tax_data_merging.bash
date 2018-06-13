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


prokka_all <- fread(...)
kaiju <-
bbmap <-
gene_reads <-
concoct <- 

names (file) <- c("contig_id", "contig_length", "sample")  # order depending on the script you used
contig_length$sample_contig_id <- do.call(paste, c(contig_length[c("sample", "contig_id")], sep = "_")) 
contig_length <- contig_length[c(4,2)]
prokka_all <- merge(prokka_all, contig_length, by = "sample_contig_id")


# add data from read mappings on genes/products

# after string adjustments for products, some read entries are duplicates (same product on same contig in same sample).
# they have to be summed up using plyr
library(plyr)
uniq_product_map2 <- ddply(product_map2,"sample_contig_id_product",numcolwise(sum))

# merge the reads with the big file
prokka_all <- merge(prokka_all, product_map2, by = "sample_contig_id_product")

# merge with omics meta data
meta_omics_small <- read.csv(file.choose(), sep = ";")
prokka_all <- merge(prokka_all, meta_omics_small, by = "sample")

# now all data to calculate the tpm (RNA) or rpm (DNA) is present: http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

prokka_all$rpk <- prokka_all$reads_per_gene/(prokka_all$gene_length/1000)

# get scaling factors for each sample
sum(subset(prokka_all,sample == "A1")$rpk)  # 104241798
sum(subset(prokka_all,sample == "A2")$rpk)  # 88531948
sum(subset(prokka_all,sample == "A3")$rpk)  # 123148955
sum(subset(prokka_all,sample == "A4")$rpk)  # 113719110
sum(subset(prokka_all,sample == "A5")$rpk)  # 76706929
sum(subset(prokka_all,sample == "A6")$rpk)  # 93840142
sum(subset(prokka_all,sample == "A7")$rpk)  # 142089935
sum(subset(prokka_all,sample == "B8")$rpk)  # 134019320
sum(subset(prokka_all,sample == "B9")$rpk)  # 87025821
sum(subset(prokka_all,sample == "B10")$rpk)  # 126025571

# add column to meta_omics_small
meta_omics_small$scale_factor <- c(104241798, 88531948, 123148955, 113719110, 76706929, 93840142, 142089935, 134019320, 87025821, 126025571)

# checkhow sum of reads per kilobase correlate to sequencing depth (in read pairs)
library(ggplot2)
ggplot(meta_omics_small, aes( x = new_day, group = treatment, lty = treatment)) +
  geom_line(aes(y = total_reads, colour = "total paired reads"))+
  geom_line(aes(y = total_reads * 2, colour = "total single reads"))+
  geom_line(aes(y = scale_factor, colour = "scale factor"))

# generate reads per million  
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

# another value, the average coverage, is created by  #35, bbmap
# it can be found under /data/jwerner/glyphosate/IMP/${i}/metagenome_coverage/covstats.txt 
# and needs to be merged by contig_id and sample 




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