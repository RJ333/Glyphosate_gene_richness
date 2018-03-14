# how to read in multiple files into a list in R

## import_multiple_cov_files_to_R
# Purpose: Import multiple cov files to the Global Environment in R
 
# set working directory
setwd("D:/1_Fachliches/gene_richness_selected_genes")
 
# list all cov files from the current directory and the subdirectories
# use the pattern argument to define a common pattern of import files with regex
list.files(pattern = "counts.cov$", recursive = TRUE) 
 
# create a list from these files
list.filenames <- list.files(pattern = "counts.cov$", recursive = TRUE)
list.filenames
 
# create an empty list that will serve as a container to receive the incoming files
list.gene_richness <- list()
 
# create a loop to read in your data
for (i in 1:length(list.filenames))
{
list.gene_richness[[i]] <- read.delim(list.filenames[i],header=FALSE)
}
 
# add the names of your data to the list
names(list.gene_richness) <- basename(list.files(pattern = "counts.cov$", recursive = TRUE))
 
# now you can index one of your tables like this
list.gene_richness$phnm_counts.cov
 
# or this
list.gene_richness[1]
 
# you can make a function out of this
import.multiple.table.files <- function(mypath, mypattern, sep, header_read, ...)
{
tmp.list.1 <- list.files(mypath, pattern = mypattern, recursive = TRUE)
tmp.list.2 <- list(length = length(tmp.list.1))
  for (i in 1:length(tmp.list.1)){tmp.list.2[[i]]<-read.table(tmp.list.1[i], sep = sep, header = header_read, ...)}
names(tmp.list.2)<-basename(list.files(mypath, pattern = mypattern, recursive = TRUE))
tmp.list.2
}
 
# use it like this
list.gene_richness <- import.multiple.table.files("D:/1_Fachliches/gene_richness_selected_genes", "counts.cov$", sep = " ", header_read = FALSE )
list.gene_length <- import.multiple.table.files("D:/1_Fachliches/gene_richness_selected_genes", "trimmed.cov$", sep = " ", header_read = FALSE )
# note: with ... we enable the function to refine the import with parameters from read.cov.
# here we define the separator of entries in the cov files to be comma.
 
# save it to the folder with your custom functions
save(import.multiple.cov.files, file = "~/R/functions/import.multiple.cov.files.RData")
 
# load it like this whenever you need it in another script with
load("D:/1_Fachliches/gene_richness_selected_genes")
 
# end


list.files(pattern = "trimmed.cov$", recursive = TRUE)