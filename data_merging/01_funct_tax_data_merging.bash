# list prokka files with sample name

# path to prokka files:
/data/jwerner/glyphosate/IMP/*/output_IMP/Analysis/annotation

# create output folder
mkdir -p /data/Rene/glyph

# go to folder with the different samples
cd /data/jwerner/glyphosate/IMP

# put for loop in a function to simply copy into the terminal
# function adds a column with the sample name and saves the modified file somewhere else
function hihi {
for sample in *
do
  echo $sample
  cd /data/jwerner/glyphosate/IMP/$sample/output_IMP/Analysis/annotation
  awk -v sample=$sample '{a=sample;}{print a"\t"$0}' annotation.filt.gff > /data/Rene/glyph/${sample}_prokka.gff
done
}


# run function
hihi

# append the 10 modified prokka files to one list

cd /data/Rene/glyph

function hoho {
for prokkas in *_prokka.gff
do	
  cat $prokkas >> all_prokkas.gff
  echo "appending $prokkas done"
done
}

# run function
hoho



# resulting prokka file with more than 1.200.000 lines is to big for excel. for convenience, some cleaning tasks should be tested in excel
# before writing code. So, the file is split roughly in half, creating two new files (xaa and xab)

split -l 700000 all_prokkas.gff

mv xaa all_prokkas_part1.gff
mv xab all_prokkas_part2.gff

# in excel
add headers to the columns: sample contig_id start_pos end_pos gene_length eC_number gene product (also including note=...) products_adjusted
add column sample_contig_id to uniquely identify contigs from different samples
add sample_contig_id_product for same reason
split single column with annotations based on ";" as sep
order annotation columns, use "gene" etc as header
delete columns with prodigal, CDS, "0", "+/-", inference, locus tag, PROKKA id
add column gene length
remove counts from genE_01, genE_02 etc --> genE
replace in product column > adjusted_products: 
  "/", "-", "_","'", """ " and whitespace with @, 
  everything lower case (=klein())
  (or later in R, if forgotten in excel)
  
# preparing gene richness data to be merged within R (single_assembly_gene_richness... repo)
mkdir /data/Rene/glyph/unnamed_reads_per_gene
mkdir /data/Rene/glyph/named_reads_per_gene

for file in ../test_folder/results/*; do cp "$file" unnamed_reads_per_gene; done

# IMPORTANT: substitute "_" with "@" in product names, which will be used as separator later
# use rename.sh
# https://stackoverflow.com/questions/50039463/replace-substitute-a-list-of-substrings-in-a-list-of-filenames-in-bash-with-awk

# if you want to run it only on e.g. the files that have more than the two "_", move those to another folder "renamed"
find ../unnamed_reads_per_gene/ -mindepth 1 -maxdepth 1 -type f -name '*_*_*_*'


running 01_combine_samtools_view_files.sh  # starting from folder "renamed"
running 02_append_and_clean_samtools_view_output.sh  # output file contig_gene_sample_reads.tsv

# in R
# put both large prokka files together
prokka_all_1 <- read.delim(file.choose(), header = TRUE)  # all_prokkas_part1_cleaned.gff
prokka_all_2 <- read.delim(file.choose(), header = TRUE)  # all_prokkas_part2_cleaned.gff
prokka_all <- rbind(prokka_all_1, prokka_all_2)
nrow(prokka_all_1)  # 700.000
nrow(prokka_all_2)  # 544.955
nrow(prokka_all)  # 1.244.955


# use contig length file for merging (version already exists, "get_contig_length_from_prokkafna.sh" in taxonomy analysis repo)
# and adjust for merging
contig_length <- read.delim(file.choose(), header = FALSE)
names (contig_length) <- c("contig_id", "contig_length", "sample")
contig_length$sample_contig_id <- do.call(paste, c(contig_length[c("sample", "contig_id")], sep = "_")) 
contig_length <- contig_length[c(4,2)]
prokka_all <- merge(prokka_all, contig_length, by = "sample_contig_id")

# if substitution forgotten:
# prokka_all$adjusted_products <- as.factor(gsub("-|'| |/", "@", prokka_all$adjusted_products))
# prokka_all$sample_contig_id_product <- as.factor(gsub("-|'| |/", "@", prokka_all$sample_contig_id_product))

# check if adjusted product names reduced ambiguity
str(prokka_all)
prokka_all$adjusted_products <- as.factor(prokka_all$adjusted_products)
length(levels(prokka_all$adjusted_products))  # 9810
length(levels(prokka_all$product))  # 10234

# add data from read mappings on genes/products
product_map <- read.delim(file.choose(), header = TRUE)  # contig_gene_sample_reads.tsv
nrow(product_map)  # 1191982
# adjust product names similar to adjusted_products
product_map$gene <- as.factor(gsub("-|'| |/", "@", product_map$gene))
product_map$gene <- tolower(product_map$gene)
product_map$sample_contig_id_product <- do.call(paste, c(product_map[c("sample", "contig","gene")], sep = "_"))
product_map2 <- product_map[c(5,4)]
product_map2$sample_contig_id_product <- as.factor(product_map2$sample_contig_id_product)

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

# check how sum of reads per kilobase correlate to sequencing depth (in read pairs)
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

# now only tax information and read normalization for that data is missing