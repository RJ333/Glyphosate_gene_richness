library(analyzeMock)

ggplot2::theme_set(ggplot2::theme_bw())

# set up folders and files
projectPath <- "/media/chrisgaby/DATA/IOW/mock-community-analysis/r_package"
testdataPath <- file.path(projectPath, "testdata")
rawPath <- file.path(testdataPath, "raw")

trainingFasta <- file.path(projectPath,
                           "tax/silva_nr_v132_train_set.fa.gz")
speciesFasta <- file.path(projectPath,
                          "tax/silva_species_assignment_v132.fa.gz")

readPattern <- c("_R1_001.fastq", "_R2_001.fastq")
fwReads <- sort(list.files(rawPath, pattern=readPattern[1], full.names=TRUE))
rvReads <- sort(list.files(rawPath, pattern=readPattern[2], full.names=TRUE))
sampleNames <- sapply(strsplit(basename(fwReads), "_"), `[`, 1)

filteredPattern <- c("_F_filt.fastq.gz", "_R_filt.fastq.gz")
fwReadsFiltered <- file.path(testdataPath, "filtered",
                             paste0(sampleNames, filteredPattern[1]))
rvReadsFiltered <- file.path(testdataPath, "filtered",
                             paste0(sampleNames, filteredPattern[2]))


createDirs(basefolder=testdataPath)
createQualityPlots(input=c(fwReads, rvReads),
                   output=file.path(testdataPath, "plots",
                                    "quality_profile_untrimmed.pdf"))
filteredReads <- dada2::filterAndTrim(fwReads, fwReadsFiltered,
                               rvReads, rvReadsFiltered,
                               truncLen=c(300, 280),
                               trimLeft=c(17, 21),
                               maxN=0,
                               maxEE=c(4, 4),
                               truncQ=2,
                               rm.phix=FALSE,
                               compress=TRUE,
                               multithread=TRUE)
checkPrimers(fwReads, rvReads,
             fwReadsFiltered, rvReadsFiltered,
             fwPrimer="CCTACGGGNBGCASCAG",
             rvPrimer="GACTACNVGGGTATCTAATCC",
             plotPath=file.path(testdataPath, "plots"))
createQualityPlots(input=c(fwReadsFiltered, rvReadsFiltered),
                   file.path(testdataPath, "plots",
                             "quality_profile_trimmed.pdf"))
seqTable <- runDada2(fwReads, fwReadsFiltered, rvReadsFiltered,
                     filteredReads, sampleNames,
                     output=c(file.path(testdataPath, "plots",
                                        "error_profile.pdf"),
                              file.path(testdataPath, "plots",
                                        "tracked_reads.pdf")))
taxa <- assignTaxInformation(seqTable, trainingFasta, speciesFasta)
# samdf <- data.frame(sample_names=sampleNames, bacterial=TRUE)
ps <- phyloseq::phyloseq(phyloseq::otu_table(seqTable, taxa_are_rows=FALSE), phyloseq::tax_table(taxa))
ps <- modifyPhyloseq(ps)

plot_bar(ps, "Genus", "Abundance", "Genus")
ggplot2::ggsave(file.path(testdataPath, "plots", "ps_barplot.pdf"))

sessionInfo()
