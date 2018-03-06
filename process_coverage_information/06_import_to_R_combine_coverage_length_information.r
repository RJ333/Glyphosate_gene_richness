library(reshape2)

# this is for phnM

# read in coverage and length file
phnM_cov<-read.csv(file.choose(),sep=" ",header=FALSE)  # phnM_counts
length_phnM<-read.csv(file.choose(),sep=" ",header=FALSE)  # phnM_contig_and_gene_length.txt_trimmed

# rename column headers
colnames(phnM_cov)<- c("contig_id","base_pos","coverage","sample","gene")
colnames(length_phnM)<- c("contig_id","gene_length","gene","contig_length")
head(phnM_cov)
head(length_phnM)

# remove everything before "="
length_phnM$gene<-sub('.*\\=', '', length_phnM$gene)

# add new column with contig_ids without appended fraction number
phnM_cov$full_contig_id<-sub('\\..*', '', phnM_cov$contig_id)

# merge coverage and length information
phnM_cov_length<-merge(phnM_cov,length_phnM,by.x="full_contig_id", by.y="contig_id")
head(phnM_cov_length)
nrow(phnM_cov)  # 218253
nrow(length_phnM)  # 43
nrow(phnM_cov_length)  # 251580

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


# count per gene
A1_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,9)]),]
A1_subset2_phnM<-subset(A1_subset_phnM, A1 != 0)
length(table(droplevels(A1_subset2_phnM$gene.x)))
# 16
A2_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,10)]),]
A2_subset2_phnM<-subset(A2_subset_phnM, A2 != 0)
length(table(droplevels(A2_subset2_phnM$gene.x)))
# 14
A3_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,11)]),]
A3_subset2_phnM<-subset(A3_subset_phnM, A3 != 0)
length(table(droplevels(A3_subset2_phnM$gene.x)))
# 20
A4_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,12)]),]
A4_subset2_phnM<-subset(A4_subset_phnM, A4 != 0)
length(table(droplevels(A4_subset2_phnM$gene.x)))
# 26
A5_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,13)]),]
A5_subset2_phnM<-subset(A5_subset_phnM, A5 != 0)
length(table(droplevels(A5_subset2_phnM$gene.x)))
# 20
A6_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,14)]),]
A6_subset2_phnM<-subset(A6_subset_phnM, A6 != 0)
length(table(droplevels(A6_subset2_phnM$gene.x)))
# 22
A7_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,15)]),]
A7_subset2_phnM<-subset(A7_subset_phnM, A7 != 0)
length(table(droplevels(A7_subset2_phnM$gene.x)))
# 26
B8_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,17)]),]
B8_subset2_phnM<-subset(B8_subset_phnM, B8 != 0)
length(table(droplevels(B8_subset2_phnM$gene.x)))
# 17
B9_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,18)]),]
B9_subset2_phnM<-subset(B9_subset_phnM, B9 != 0)
length(table(droplevels(B9_subset2_phnM$gene.x)))
# 12
B10_subset_phnM<-cast_phnM_cov_length[!duplicated(cast_phnM_cov_length[,c(3,7,16)]),]
B10_subset2_phnM<-subset(B10_subset_phnM, B10 != 0)
length(table(droplevels(B10_subset2_phnM$gene.x)))
# 19





# this is for phnE

# read in coverage and length file
phnE_cov<-read.csv(file.choose(),sep=" ",header=FALSE)  # phnE_counts
length_phnE<-read.csv(file.choose(),sep=" ",header=FALSE)  # phnE_contig_and_gene_length.txt_trimmed

# rename column headers
colnames(phnE_cov)<- c("contig_id","base_pos","coverage","sample","gene")
colnames(length_phnE)<- c("contig_id","gene_length","gene","contig_length")
head(phnE_cov)
head(length_phnE)

# remove everything before "="
length_phnE$gene<-sub('.*\\=', '', length_phnE$gene)

# add new column with contig_ids without appended fraction number
phnE_cov$full_contig_id<-sub('\\..*', '', phnE_cov$contig_id)

# merge coverage and length information
phnE_cov_length<-merge(phnE_cov,length_phnE,by.x="full_contig_id", by.y="contig_id")
head(phnE_cov_length)
nrow(phnE_cov)  # 661771
nrow(length_phnE)  # 94
nrow(phnE_cov_length)  # 654498

# turn character into factors
str(phnE_cov_length)
phnE_cov_length$full_contig_id<-as.factor(phnE_cov_length$full_contig_id)
phnE_cov_length$gene.y<-as.factor(phnE_cov_length$gene.y)
str(phnE_cov_length)
# check how many contigs have been found, note that not all from the length file must be present as 
# some might be below the coverage threshold (k141_136374,k141_79399)
levels(phnE_cov_length$full_contig_id)
# 62
# add sequence value
phnE_cov_length$ID <- seq.int(nrow(phnE_cov_length))
cast_phnE_cov_length<-dcast(phnE_cov_length,ID+full_contig_id+contig_id+contig_length+gene_length+base_pos+gene.x+gene.y~sample, value.var="coverage")
head(cast_phnE_cov_length)
str(cast_phnE_cov_length)
# turning NAs into 0s
cast_phnE_cov_length[is.na(cast_phnE_cov_length)] <- 0

# count per base
colSums(cast_phnE_cov_length[,c(9:18)] != 0)
plot(colSums(cast_phnE_cov_length[,c(9:18)] != 0))
# count per contig
A1_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,9)]),]
A1_subset2_phnE<-subset(A1_subset_phnE, A1 != 0)
length(table(droplevels(A1_subset2_phnE$full_contig_id)))
# 24
A2_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,10)]),]
A2_subset2_phnE<-subset(A2_subset_phnE, A2 != 0)
length(table(droplevels(A2_subset2_phnE$full_contig_id)))
# 19
A3_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,11)]),]
A3_subset2_phnE<-subset(A3_subset_phnE, A3 != 0)
length(table(droplevels(A3_subset2_phnE$full_contig_id)))
# 25
A4_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,12)]),]
A4_subset2_phnE<-subset(A4_subset_phnE, A4 != 0)
length(table(droplevels(A4_subset2_phnE$full_contig_id)))
# 35
A5_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,13)]),]
A5_subset2_phnE<-subset(A5_subset_phnE, A5 != 0)
length(table(droplevels(A5_subset2_phnE$full_contig_id)))
# 30
A6_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,14)]),]
A6_subset2_phnE<-subset(A6_subset_phnE, A6 != 0)
length(table(droplevels(A6_subset2_phnE$full_contig_id)))
# 28
A7_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,15)]),]
A7_subset2_phnE<-subset(A7_subset_phnE, A7 != 0)
length(table(droplevels(A7_subset2_phnE$full_contig_id)))
# 30
B8_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,17)]),]
B8_subset2_phnE<-subset(B8_subset_phnE, B8 != 0)
length(table(droplevels(B8_subset2_phnE$full_contig_id)))
# 33
B9_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,18)]),]
B9_subset2_phnE<-subset(B9_subset_phnE, B9 != 0)
length(table(droplevels(B9_subset2_phnE$full_contig_id)))
# 25
B10_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,16)]),]
B10_subset2_phnE<-subset(B10_subset_phnE, B10 != 0)
length(table(droplevels(B10_subset2_phnE$full_contig_id)))
# 32


# count per gene
A1_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,9)]),]
A1_subset2_phnE<-subset(A1_subset_phnE, A1 != 0)
length(table(droplevels(A1_subset2_phnE$gene.x)))
# 34
A2_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,10)]),]
A2_subset2_phnE<-subset(A2_subset_phnE, A2 != 0)
length(table(droplevels(A2_subset2_phnE$gene.x)))
# 28
A3_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,11)]),]
A3_subset2_phnE<-subset(A3_subset_phnE, A3 != 0)
length(table(droplevels(A3_subset2_phnE$gene.x)))
# 38
A4_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,12)]),]
A4_subset2_phnE<-subset(A4_subset_phnE, A4 != 0)
length(table(droplevels(A4_subset2_phnE$gene.x)))
# 55
A5_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,13)]),]
A5_subset2_phnE<-subset(A5_subset_phnE, A5 != 0)
length(table(droplevels(A5_subset2_phnE$gene.x)))
# 48
A6_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,14)]),]
A6_subset2_phnE<-subset(A6_subset_phnE, A6 != 0)
length(table(droplevels(A6_subset2_phnE$gene.x)))
# 46
A7_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,15)]),]
A7_subset2_phnE<-subset(A7_subset_phnE, A7 != 0)
length(table(droplevels(A7_subset2_phnE$gene.x)))
# 50
B8_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,17)]),]
B8_subset2_phnE<-subset(B8_subset_phnE, B8 != 0)
length(table(droplevels(B8_subset2_phnE$gene.x)))
# 53
B9_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,18)]),]
B9_subset2_phnE<-subset(B9_subset_phnE, B9 != 0)
length(table(droplevels(B9_subset2_phnE$gene.x)))
# 44
B10_subset_phnE<-cast_phnE_cov_length[!duplicated(cast_phnE_cov_length[,c(3,7,16)]),]
B10_subset2_phnE<-subset(B10_subset_phnE, B10 != 0)
length(table(droplevels(B10_subset2_phnE$gene.x)))
# 54


# this is for gyrA

# read in coverage and length file
gyrA_cov<-read.csv(file.choose(),sep=" ",header=FALSE)  # gyrA_counts
length_gyrA<-read.csv(file.choose(),sep=" ",header=FALSE)  # gyrA_contig_and_gene_length.txt_trimmed

# rename column headers
colnames(gyrA_cov)<- c("contig_id","base_pos","coverage","sample","gene")
colnames(length_gyrA)<- c("contig_id","gene_length","gene","contig_length")
head(gyrA_cov)
head(length_gyrA)

# remove everything before "="
length_gyrA$gene<-sub('.*\\=', '', length_gyrA$gene)

# add new column with contig_ids without appended fraction number
gyrA_cov$full_contig_id<-sub('\\..*', '', gyrA_cov$contig_id)

# merge coverage and length information
gyrA_cov_length<-merge(gyrA_cov,length_gyrA,by.x="full_contig_id", by.y="contig_id")
head(gyrA_cov_length)
nrow(gyrA_cov)  # 5499208
nrow(length_gyrA)  # 100
nrow(gyrA_cov_length)  # 3385220

# turn character into factors
str(gyrA_cov_length)
gyrA_cov_length$full_contig_id<-as.factor(gyrA_cov_length$full_contig_id)
gyrA_cov_length$gene.y<-as.factor(gyrA_cov_length$gene.y)
str(gyrA_cov_length)
# check how many contigs have been found, note that not all from the length file must be present as 
# some might be below the coverage threshold (k141_136374,k141_79399)
length(levels(gyrA_cov_length$full_contig_id))
# 91

# add sequence value
gyrA_cov_length$ID <- seq.int(nrow(gyrA_cov_length))
cast_gyrA_cov_length<-dcast(gyrA_cov_length,ID+full_contig_id+contig_id+contig_length+gene_length+base_pos+gene.x+gene.y~sample, value.var="coverage")
head(cast_gyrA_cov_length)
str(cast_gyrA_cov_length)
# turning NAs into 0s
cast_gyrA_cov_length[is.na(cast_gyrA_cov_length)] <- 0

# count per base
colSums(cast_gyrA_cov_length[,c(9:18)] != 0)
plot(colSums(cast_gyrA_cov_length[,c(9:18)] != 0))
# count per contig
A1_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,9)]),]
A1_subset2_gyrA<-subset(A1_subset_gyrA, A1 != 0)
length(table(droplevels(A1_subset2_gyrA$full_contig_id)))
# 41
A2_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,10)]),]
A2_subset2_gyrA<-subset(A2_subset_gyrA, A2 != 0)
length(table(droplevels(A2_subset2_gyrA$full_contig_id)))
# 37
A3_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,11)]),]
A3_subset2_gyrA<-subset(A3_subset_gyrA, A3 != 0)
length(table(droplevels(A3_subset2_gyrA$full_contig_id)))
# 45
A4_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,12)]),]
A4_subset2_gyrA<-subset(A4_subset_gyrA, A4 != 0)
length(table(droplevels(A4_subset2_gyrA$full_contig_id)))
# 55
A5_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,13)]),]
A5_subset2_gyrA<-subset(A5_subset_gyrA, A5 != 0)
length(table(droplevels(A5_subset2_gyrA$full_contig_id)))
# 44
A6_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,14)]),]
A6_subset2_gyrA<-subset(A6_subset_gyrA, A6 != 0)
length(table(droplevels(A6_subset2_gyrA$full_contig_id)))
# 44
A7_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,15)]),]
A7_subset2_gyrA<-subset(A7_subset_gyrA, A7 != 0)
length(table(droplevels(A7_subset2_gyrA$full_contig_id)))
# 56
B8_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,17)]),]
B8_subset2_gyrA<-subset(B8_subset_gyrA, B8 != 0)
length(table(droplevels(B8_subset2_gyrA$full_contig_id)))
# 48
B9_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,18)]),]
B9_subset2_gyrA<-subset(B9_subset_gyrA, B9 != 0)
length(table(droplevels(B9_subset2_gyrA$full_contig_id)))
# 39
B10_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,16)]),]
B10_subset2_gyrA<-subset(B10_subset_gyrA, B10 != 0)
length(table(droplevels(B10_subset2_gyrA$full_contig_id)))
# 55


# count per gene
A1_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,9)]),]
A1_subset2_gyrA<-subset(A1_subset_gyrA, A1 != 0)
length(table(droplevels(A1_subset2_gyrA$gene.x)))
# 46
A2_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,10)]),]
A2_subset2_gyrA<-subset(A2_subset_gyrA, A2 != 0)
length(table(droplevels(A2_subset2_gyrA$gene.x)))
# 41
A3_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,11)]),]
A3_subset2_gyrA<-subset(A3_subset_gyrA, A3 != 0)
length(table(droplevels(A3_subset2_gyrA$gene.x)))
# 50
A4_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,12)]),]
A4_subset2_gyrA<-subset(A4_subset_gyrA, A4 != 0)
length(table(droplevels(A4_subset2_gyrA$gene.x)))
# 60
A5_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,13)]),]
A5_subset2_gyrA<-subset(A5_subset_gyrA, A5 != 0)
length(table(droplevels(A5_subset2_gyrA$gene.x)))
# 48
A6_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,14)]),]
A6_subset2_gyrA<-subset(A6_subset_gyrA, A6 != 0)
length(table(droplevels(A6_subset2_gyrA$gene.x)))
# 49
A7_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,15)]),]
A7_subset2_gyrA<-subset(A7_subset_gyrA, A7 != 0)
length(table(droplevels(A7_subset2_gyrA$gene.x)))
# 61
B8_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,17)]),]
B8_subset2_gyrA<-subset(B8_subset_gyrA, B8 != 0)
length(table(droplevels(B8_subset2_gyrA$gene.x)))
# 52
B9_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,18)]),]
B9_subset2_gyrA<-subset(B9_subset_gyrA, B9 != 0)
length(table(droplevels(B9_subset2_gyrA$gene.x)))
# 43
B10_subset_gyrA<-cast_gyrA_cov_length[!duplicated(cast_gyrA_cov_length[,c(3,7,16)]),]
B10_subset2_gyrA<-subset(B10_subset_gyrA, B10 != 0)
length(table(droplevels(B10_subset2_gyrA$gene.x)))
# 59