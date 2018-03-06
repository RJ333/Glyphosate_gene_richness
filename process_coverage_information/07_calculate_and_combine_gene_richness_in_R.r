# all genes used here are longer than 400 bp, the contigs they are on are longer than 1000 bp and the minimum coverage per base is 5 x 
# this script brings the quality checked data into the right format to display the gene richness per gene per sample

# no clue how to do this with a loop or lapply for all my genes, so I copy the code for each gene
library(reshape2)
library(ggplot2)

# this is for gyrA

test_gyrA<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_gyrA)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_gyrA<-test_gyrA[,c(-2,-3)]

# check for unique contigs
unique(test_gyrA)
test_gyrA_uniq<-test_gyrA[!duplicated(test_gyrA[,c(1:3)]),]
cast_uniq_gyrA<-dcast(test_gyrA_uniq,contig_id+gene~sample)

# turn NAs into 0 and gyrA genes into 1
cast_uniq_gyrA[is.na(cast_uniq_gyrA)] <- 0
cast_uniq_gyrA[,c(3:12)]<-data.frame(lapply(cast_uniq_gyrA[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_gyrA[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_gyrA[,c(3:12)])

# create new dataset with extra column gene name
gyrA_richness<-as.data.frame(colSums(cast_uniq_gyrA[,c(3:12)]))
colnames(gyrA_richness)<-"richness"

# prepare vectors containing meta data
gyrA_gene_vector<-"gyrA"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
gyrA_richness$gene<-gyrA_gene_vector
gyrA_richness$days<-day_vector
gyrA_richness$treatment<-treatment_vector

# this is for gyrB

test_gyrB<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_gyrB)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_gyrB<-test_gyrB[,c(-2,-3)]

# check for unique contigs
unique(test_gyrB)
test_gyrB_uniq<-test_gyrB[!duplicated(test_gyrB[,c(1:3)]),]
cast_uniq_gyrB<-dcast(test_gyrB_uniq,contig_id+gene~sample)

# turn NAs into 0 and gyrB genes into 1
cast_uniq_gyrB[is.na(cast_uniq_gyrB)] <- 0
cast_uniq_gyrB[,c(3:12)]<-data.frame(lapply(cast_uniq_gyrB[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_gyrB[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_gyrB[,c(3:12)])

# create new dataset with extra column gene name
gyrB_richness<-as.data.frame(colSums(cast_uniq_gyrB[,c(3:12)]))
colnames(gyrB_richness)<-"richness"

# prepare vectors containing meta data
gyrB_gene_vector<-"gyrB"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
gyrB_richness$gene<-gyrB_gene_vector
gyrB_richness$days<-day_vector
gyrB_richness$treatment<-treatment_vector

# this is for phnC

test_phnC<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnC)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnC<-test_phnC[,c(-2,-3)]

# check for unique contigs
unique(test_phnC)
test_phnC_uniq<-test_phnC[!duplicated(test_phnC[,c(1:3)]),]
cast_uniq_phnC<-dcast(test_phnC_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnC genes into 1
cast_uniq_phnC[is.na(cast_uniq_phnC)] <- 0
cast_uniq_phnC[,c(3:12)]<-data.frame(lapply(cast_uniq_phnC[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnC[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnC[,c(3:12)])

# create new dataset with extra column gene name
phnC_richness<-as.data.frame(colSums(cast_uniq_phnC[,c(3:12)]))
colnames(phnC_richness)<-"richness"

# prepare vectors containing meta data
phnC_gene_vector<-"phnC"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnC_richness$gene<-phnC_gene_vector
phnC_richness$days<-day_vector
phnC_richness$treatment<-treatment_vector

# this is for phnD

test_phnD<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnD)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnD<-test_phnD[,c(-2,-3)]

# check for unique contigs
unique(test_phnD)
test_phnD_uniq<-test_phnD[!duplicated(test_phnD[,c(1:3)]),]
cast_uniq_phnD<-dcast(test_phnD_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnD genes into 1
cast_uniq_phnD[is.na(cast_uniq_phnD)] <- 0
cast_uniq_phnD[,c(3:12)]<-data.frame(lapply(cast_uniq_phnD[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnD[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnD[,c(3:12)])

# create new dataset with extra column gene name
phnD_richness<-as.data.frame(colSums(cast_uniq_phnD[,c(3:12)]))
colnames(phnD_richness)<-"richness"

# prepare vectors containing meta data
phnD_gene_vector<-"phnD"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnD_richness$gene<-phnD_gene_vector
phnD_richness$days<-day_vector
phnD_richness$treatment<-treatment_vector

# this is for phnE

test_phnE<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnE)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnE<-test_phnE[,c(-2,-3)]

# check for unique contigs
unique(test_phnE)
test_phnE_uniq<-test_phnE[!duplicated(test_phnE[,c(1:3)]),]
cast_uniq_phnE<-dcast(test_phnE_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnE genes into 1
cast_uniq_phnE[is.na(cast_uniq_phnE)] <- 0
cast_uniq_phnE[,c(3:12)]<-data.frame(lapply(cast_uniq_phnE[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnE[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnE[,c(3:12)])

# create new dataset with extra column gene name
phnE_richness<-as.data.frame(colSums(cast_uniq_phnE[,c(3:12)]))
colnames(phnE_richness)<-"richness"

# prepare vectors containing meta data
phnE_gene_vector<-"phnE"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnE_richness$gene<-phnE_gene_vector
phnE_richness$days<-day_vector
phnE_richness$treatment<-treatment_vector

# this is for phnF

test_phnF<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnF)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnF<-test_phnF[,c(-2,-3)]

# check for unique contigs
unique(test_phnF)
test_phnF_uniq<-test_phnF[!duplicated(test_phnF[,c(1:3)]),]
cast_uniq_phnF<-dcast(test_phnF_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnF genes into 1
cast_uniq_phnF[is.na(cast_uniq_phnF)] <- 0
cast_uniq_phnF[,c(3:12)]<-data.frame(lapply(cast_uniq_phnF[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnF[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnF[,c(3:12)])

# create new dataset with extra column gene name
phnF_richness<-as.data.frame(colSums(cast_uniq_phnF[,c(3:12)]))
colnames(phnF_richness)<-"richness"

# prepare vectors containing meta data
phnF_gene_vector<-"phnF"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnF_richness$gene<-phnF_gene_vector
phnF_richness$days<-day_vector
phnF_richness$treatment<-treatment_vector

# this is for phnG

test_phnG<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnG)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnG<-test_phnG[,c(-2,-3)]

# check for unique contigs
unique(test_phnG)
test_phnG_uniq<-test_phnG[!duplicated(test_phnG[,c(1:3)]),]
cast_uniq_phnG<-dcast(test_phnG_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnG genes into 1
cast_uniq_phnG[is.na(cast_uniq_phnG)] <- 0
cast_uniq_phnG[,c(3:12)]<-data.frame(lapply(cast_uniq_phnG[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnG[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnG[,c(3:12)])

# create new dataset with extra column gene name
phnG_richness<-as.data.frame(colSums(cast_uniq_phnG[,c(3:12)]))
colnames(phnG_richness)<-"richness"

# prepare vectors containing meta data
phnG_gene_vector<-"phnG"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnG_richness$gene<-phnG_gene_vector
phnG_richness$days<-day_vector
phnG_richness$treatment<-treatment_vector

# this is for phnH

test_phnH<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnH)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnH<-test_phnH[,c(-2,-3)]

# check for unique contigs
unique(test_phnH)
test_phnH_uniq<-test_phnH[!duplicated(test_phnH[,c(1:3)]),]
cast_uniq_phnH<-dcast(test_phnH_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnH genes into 1
cast_uniq_phnH[is.na(cast_uniq_phnH)] <- 0
cast_uniq_phnH[,c(3:12)]<-data.frame(lapply(cast_uniq_phnH[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnH[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnH[,c(3:12)])

# create new dataset with extra column gene name
phnH_richness<-as.data.frame(colSums(cast_uniq_phnH[,c(3:12)]))
colnames(phnH_richness)<-"richness"

# prepare vectors containing meta data
phnH_gene_vector<-"phnH"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnH_richness$gene<-phnH_gene_vector
phnH_richness$days<-day_vector
phnH_richness$treatment<-treatment_vector

# this is for phnI

test_phnI<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnI)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnI<-test_phnI[,c(-2,-3)]

# check for unique contigs
unique(test_phnI)
test_phnI_uniq<-test_phnI[!duplicated(test_phnI[,c(1:3)]),]
cast_uniq_phnI<-dcast(test_phnI_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnI genes into 1
cast_uniq_phnI[is.na(cast_uniq_phnI)] <- 0
cast_uniq_phnI[,c(3:12)]<-data.frame(lapply(cast_uniq_phnI[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnI[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnI[,c(3:12)])

# create new dataset with extra column gene name
phnI_richness<-as.data.frame(colSums(cast_uniq_phnI[,c(3:12)]))
colnames(phnI_richness)<-"richness"

# prepare vectors containing meta data
phnI_gene_vector<-"phnI"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnI_richness$gene<-phnI_gene_vector
phnI_richness$days<-day_vector
phnI_richness$treatment<-treatment_vector

# this is for phnJ

test_phnJ<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnJ)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnJ<-test_phnJ[,c(-2,-3)]

# check for unique contigs
unique(test_phnJ)
test_phnJ_uniq<-test_phnJ[!duplicated(test_phnJ[,c(1:3)]),]
cast_uniq_phnJ<-dcast(test_phnJ_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnJ genes into 1
cast_uniq_phnJ[is.na(cast_uniq_phnJ)] <- 0
cast_uniq_phnJ[,c(3:12)]<-data.frame(lapply(cast_uniq_phnJ[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnJ[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnJ[,c(3:12)])

# create new dataset with extra column gene name
phnJ_richness<-as.data.frame(colSums(cast_uniq_phnJ[,c(3:12)]))
colnames(phnJ_richness)<-"richness"

# prepare vectors containing meta data
phnJ_gene_vector<-"phnJ"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnJ_richness$gene<-phnJ_gene_vector
phnJ_richness$days<-day_vector
phnJ_richness$treatment<-treatment_vector

# this is for phnK

test_phnK<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnK)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnK<-test_phnK[,c(-2,-3)]

# check for unique contigs
unique(test_phnK)
test_phnK_uniq<-test_phnK[!duplicated(test_phnK[,c(1:3)]),]
cast_uniq_phnK<-dcast(test_phnK_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnK genes into 1
cast_uniq_phnK[is.na(cast_uniq_phnK)] <- 0
cast_uniq_phnK[,c(3:12)]<-data.frame(lapply(cast_uniq_phnK[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnK[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnK[,c(3:12)])

# create new dataset with extra column gene name
phnK_richness<-as.data.frame(colSums(cast_uniq_phnK[,c(3:12)]))
colnames(phnK_richness)<-"richness"

# prepare vectors containing meta data
phnK_gene_vector<-"phnK"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnK_richness$gene<-phnK_gene_vector
phnK_richness$days<-day_vector
phnK_richness$treatment<-treatment_vector

# this is for phnL

test_phnL<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnL)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnL<-test_phnL[,c(-2,-3)]

# check for unique contigs
unique(test_phnL)
test_phnL_uniq<-test_phnL[!duplicated(test_phnL[,c(1:3)]),]
cast_uniq_phnL<-dcast(test_phnL_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnL genes into 1
cast_uniq_phnL[is.na(cast_uniq_phnL)] <- 0
cast_uniq_phnL[,c(3:12)]<-data.frame(lapply(cast_uniq_phnL[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnL[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnL[,c(3:12)])

# create new dataset with extra column gene name
phnL_richness<-as.data.frame(colSums(cast_uniq_phnL[,c(3:12)]))
colnames(phnL_richness)<-"richness"

# prepare vectors containing meta data
phnL_gene_vector<-"phnL"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnL_richness$gene<-phnL_gene_vector
phnL_richness$days<-day_vector
phnL_richness$treatment<-treatment_vector

# this is for phnM

test_phnM<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnM)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnM<-test_phnM[,c(-2,-3)]

# check for unique contigs
unique(test_phnM)
test_phnM_uniq<-test_phnM[!duplicated(test_phnM[,c(1:3)]),]
cast_uniq_phnM<-dcast(test_phnM_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnM genes into 1
cast_uniq_phnM[is.na(cast_uniq_phnM)] <- 0
cast_uniq_phnM[,c(3:12)]<-data.frame(lapply(cast_uniq_phnM[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnM[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnM[,c(3:12)])

# create new dataset with extra column gene name
phnM_richness<-as.data.frame(colSums(cast_uniq_phnM[,c(3:12)]))
colnames(phnM_richness)<-"richness"

# prepare vectors containing meta data
phnM_gene_vector<-"phnM"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnM_richness$gene<-phnM_gene_vector
phnM_richness$days<-day_vector
phnM_richness$treatment<-treatment_vector

# this is for phnN

test_phnN<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnN)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnN<-test_phnN[,c(-2,-3)]

# check for unique contigs
unique(test_phnN)
test_phnN_uniq<-test_phnN[!duplicated(test_phnN[,c(1:3)]),]
cast_uniq_phnN<-dcast(test_phnN_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnN genes into 1
cast_uniq_phnN[is.na(cast_uniq_phnN)] <- 0
cast_uniq_phnN[,c(3:12)]<-data.frame(lapply(cast_uniq_phnN[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnN[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnN[,c(3:12)])

# create new dataset with extra column gene name
phnN_richness<-as.data.frame(colSums(cast_uniq_phnN[,c(3:12)]))
colnames(phnN_richness)<-"richness"

# prepare vectors containing meta data
phnN_gene_vector<-"phnN"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnN_richness$gene<-phnN_gene_vector
phnN_richness$days<-day_vector
phnN_richness$treatment<-treatment_vector

# this is for phnP

test_phnP<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_phnP)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_phnP<-test_phnP[,c(-2,-3)]

# check for unique contigs
unique(test_phnP)
test_phnP_uniq<-test_phnP[!duplicated(test_phnP[,c(1:3)]),]
cast_uniq_phnP<-dcast(test_phnP_uniq,contig_id+gene~sample)

# turn NAs into 0 and phnP genes into 1
cast_uniq_phnP[is.na(cast_uniq_phnP)] <- 0
cast_uniq_phnP[,c(3:12)]<-data.frame(lapply(cast_uniq_phnP[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_phnP[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_phnP[,c(3:12)])

# create new dataset with extra column gene name
phnP_richness<-as.data.frame(colSums(cast_uniq_phnP[,c(3:12)]))
colnames(phnP_richness)<-"richness"

# prepare vectors containing meta data
phnP_gene_vector<-"phnP"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
phnP_richness$gene<-phnP_gene_vector
phnP_richness$days<-day_vector
phnP_richness$treatment<-treatment_vector

# this is for recA

test_recA<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_recA)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_recA<-test_recA[,c(-2,-3)]

# check for unique contigs
unique(test_recA)
test_recA_uniq<-test_recA[!duplicated(test_recA[,c(1:3)]),]
cast_uniq_recA<-dcast(test_recA_uniq,contig_id+gene~sample)

# turn NAs into 0 and recA genes into 1
cast_uniq_recA[is.na(cast_uniq_recA)] <- 0
cast_uniq_recA[,c(3:12)]<-data.frame(lapply(cast_uniq_recA[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_recA[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_recA[,c(3:12)])

# create new dataset with extra column gene name
recA_richness<-as.data.frame(colSums(cast_uniq_recA[,c(3:12)]))
colnames(recA_richness)<-"richness"

# prepare vectors containing meta data
recA_gene_vector<-"recA"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
recA_richness$gene<-recA_gene_vector
recA_richness$days<-day_vector
recA_richness$treatment<-treatment_vector

# this is for rpoC

test_rpoC<-read.csv(file.choose(),sep=" ",header=FALSE)

# rename column headers
colnames(test_rpoC)<- c("contig_id","base_pos","coverage","sample","gene")

# remove base_pos and coverage
test_rpoC<-test_rpoC[,c(-2,-3)]

# check for unique contigs
unique(test_rpoC)
test_rpoC_uniq<-test_rpoC[!duplicated(test_rpoC[,c(1:3)]),]
cast_uniq_rpoC<-dcast(test_rpoC_uniq,contig_id+gene~sample)

# turn NAs into 0 and rpoC genes into 1
cast_uniq_rpoC[is.na(cast_uniq_rpoC)] <- 0
cast_uniq_rpoC[,c(3:12)]<-data.frame(lapply(cast_uniq_rpoC[,c(3:12)], function(x) as.numeric(x!="0")))

# check if each gene is present in at least one sample rowsums > 0
rowSums(cast_uniq_rpoC[,c(3:12)])
# sum gene occurrences per samples, this is the gene richness
colSums(cast_uniq_rpoC[,c(3:12)])

# create new dataset with extra column gene name
rpoC_richness<-as.data.frame(colSums(cast_uniq_rpoC[,c(3:12)]))
colnames(rpoC_richness)<-"richness"

# prepare vectors containing meta data
rpoC_gene_vector<-"rpoC"
day_vector<-(c(0,3,7,14,22,43,71,71,0,22))
treatment_vector<-rep(c("glyph","control"),times=c(7,3))

# and add them to new data.frame
rpoC_richness$gene<-rpoC_gene_vector
rpoC_richness$days<-day_vector
rpoC_richness$treatment<-treatment_vector

# turn gene data.frames per rbind into a melt type data.frame for ggplo2
richness_combined<-rbind(gyrA_richness,gyrB_richness,phnC_richness,phnD_richness,phnE_richness,phnF_richness,phnG_richness,phnH_richness,phnI_richness,phnJ_richness,phnK_richness,phnL_richness,phnM_richness,phnN_richness,phnP_richness,recA_richness,rpoC_richness)
# row names are messed up due to being not unique, so rbind changed them. Doesn't matter as each sample is defined per day and treatment

phn_subset<-subset(richness_combined,grepl("phn", gene))
housekeeping_subset<-subset(richness_combined,grepl("gyr|rec|rpo", gene))

ggplot(phn_subset,aes(x=days,y=richness))+
  geom_line(aes(colour=gene))+
  facet_wrap(~treatment)
  
ggplot(housekeeping_subset,aes(x=days,y=richness))+
  geom_line(aes(colour=gene))+
  facet_wrap(~treatment)
  
  
#turn into relative numbers for better comparison
cast_richness_combined<-dcast(richness_combined,days+treatment~gene,value.var="richness")
write.csv(cast_richness_combined,file="cast_richness_combined.csv")

# in excel relative values calculated
cast_relative_richness<-read.csv(row.names=1, file.choose(),sep=";")
melt_relative_richness<-melt(cast_relative_richness,id=c("treatment", "days"))

phn_subset2<-subset(melt_relative_richness,grepl("phn", variable))
housekeeping_subset2<-subset(melt_relative_richness,grepl("gyr|rec|rpo", variable))

ggplot(phn_subset2,aes(x=days,y=value))+
  geom_line(aes(colour=variable))+
  facet_wrap(~treatment)
  
ggplot(housekeeping_subset2,aes(x=days,y=value))+
  geom_line(aes(colour=variable))+
  facet_wrap(~treatment)


ggplot(melt_relative_richness,aes(x=days,y=value))+
  geom_line(aes(colour=variable))+
  facet_wrap(~treatment)
