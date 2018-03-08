# this is for phnM

# rename column headers
colnames_richness <- c("contig_id", "base_pos", "coverage", "sample", "gene")
list.gene_richness <- lapply(list.gene_richness, setNames, colnames_richness)

colnames_length <- c("contig_id", "gene_length", "gene", "contig_length")
list.gene_length <- lapply(list.gene_length, setNames, colnames_length)

str(list.gene_richness)
str(list.gene_length)

# remove everything before "="
# this approach only works if you modify existing columns
list.gene_length <- lapply(rapply(list.gene_length, function(x) 
  sub('.*\\=', '', x), how = "list"), 
  as.data.frame)

# add new column with contig_ids without appended fraction number
# this approach is suited to ADD a new column
list.gene_richness <- lapply(list.gene_richness, function(x) 
{x$full_contig_id <- sub('\\..*', '', x$contig_id);
x})


# merge coverage and length information
# without SIMPLIFY it creates a mess of data inside the new object
list.gene_merged <- list()
list.gene_merged <- mapply(function(x,y) {
  as.data.frame(merge(x, y, by.x = "full_contig_id", by.y = "contig_id"))
  }, x = list.gene_richness, y = list.gene_length, SIMPLIFY = FALSE)

# check the outcome, [13] is phnM here, don't forget the head.list function  
  
head.list(list.gene_merged)
lapply(list.gene_richness,nrow) 
lapply(list.gene_richness[13],nrow)  # 218253
lapply(list.gene_length[13],nrow)  # 43
sapply(list.gene_merged[13],nrow)  # 251580
sapply(list.gene_merged,str)


############################################ not finished below

# turn character into factors
str(phnM_cov_length)
phnM_cov_length$full_contig_id<-as.factor(phnM_cov_length$full_contig_id)
phnM_cov_length$gene.y<-as.factor(phnM_cov_length$gene.y)
str(phnM_cov_length)
# check how many contigs have been found, note that not all from the length file must be present as 
# some might be below the coverage threshold (k141_136374,k141_79399)
levels(phnM_cov_length$full_contig_id)
# 36
# add sequence value
phnM_cov_length$ID <- seq.int(nrow(phnM_cov_length))
cast_phnM_cov_length<-dcast(phnM_cov_length,ID+full_contig_id+contig_id+contig_length+gene_length+base_pos+gene.x+gene.y~sample, value.var="coverage")
head(cast_phnM_cov_length)
str(cast_phnM_cov_length)
# turning NAs into 0s
cast_phnM_cov_length[is.na(cast_phnM_cov_length)] <- 0

# count per base
colSums(cast_phnM_cov_length[,c(9:18)] != 0)
plot(colSums(cast_phnM_cov_length[,c(9:18)] != 0))
# count per contig
A1_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,9)]),]
A1_subset2_phnM<-subset(A1_subset_phnM, A1 != 0)
length(table(droplevels(A1_subset2_phnM$full_contig_id)))
# 13
A2_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,10)]),]
A2_subset2_phnM<-subset(A2_subset_phnM, A2 != 0)
length(table(droplevels(A2_subset2_phnM$full_contig_id)))
# 11
A3_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,11)]),]
A3_subset2_phnM<-subset(A3_subset_phnM, A3 != 0)
length(table(droplevels(A3_subset2_phnM$full_contig_id)))
# 16
A4_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,12)]),]
A4_subset2_phnM<-subset(A4_subset_phnM, A4 != 0)
length(table(droplevels(A4_subset2_phnM$full_contig_id)))
# 22
A5_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,13)]),]
A5_subset2_phnM<-subset(A5_subset_phnM, A5 != 0)
length(table(droplevels(A5_subset2_phnM$full_contig_id)))
# 16
A6_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,14)]),]
A6_subset2_phnM<-subset(A6_subset_phnM, A6 != 0)
length(table(droplevels(A6_subset2_phnM$full_contig_id)))
# 17
A7_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,15)]),]
A7_subset2_phnM<-subset(A7_subset_phnM, A7 != 0)
length(table(droplevels(A7_subset2_phnM$full_contig_id)))
# 19
B8_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,17)]),]
B8_subset2_phnM<-subset(B8_subset_phnM, B8 != 0)
length(table(droplevels(B8_subset2_phnM$full_contig_id)))
# 16
B9_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,18)]),]
B9_subset2_phnM<-subset(B9_subset_phnM, B9 != 0)
length(table(droplevels(B9_subset2_phnM$full_contig_id)))
# 11
B10_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,16)]),]
B10_subset2_phnM<-subset(B10_subset_phnM, B10 != 0)
length(table(droplevels(B10_subset2_phnM$full_contig_id)))
# 17