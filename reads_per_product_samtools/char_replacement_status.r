what does the current version of the prokka process substitute?

## unify product names , all to lower
print("unifying product annotations")
prokka_select$product2 <- tolower(prokka_select$product2)

## replace white lines and special characters in product names with "@"
prokka_select$product2 <- gsub("-|'| |/|:", "@", prokka_select$product2)  
# is "|" a problem?
prokka_select$product2 <- gsub("|", "@", fixed = TRUE)

# unify gene numbering e.g. genE_01, genE_02 etc --> genE
print("unifying gene names")
prokka_select$gene <- gsub("_.*$","", prokka_select$gene)




what was done to the processed data?

replaced () | "" / ; - _
still present [] '' :

due the same to contig_gene_sample_reads.tsv

. _ was missing