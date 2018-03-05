#!/usr/bin/env R

get_variable_name <- function(x) {
  cat(deparse(substitute(x)))
}

setwd("/data/projects/scripts/Glyphosate_gene_richness/read_abundance")
a1 <- read.table("data/A1.txt", sep = "", header = FALSE)
a2 <- read.table("data/A2.txt", sep = "", header = FALSE)
a3 <- read.table("data/A3.txt", sep = "", header = FALSE)
a4 <- read.table("data/A4.txt", sep = "", header = FALSE)
a5 <- read.table("data/A5.txt", sep = "", header = FALSE)
a6 <- read.table("data/A6.txt", sep = "", header = FALSE)
a7 <- read.table("data/A7.txt", sep = "", header = FALSE)

b8 <- read.table("data/B8.txt", sep = "", header = FALSE)
b9 <- read.table("data/B9.txt", sep = "", header = FALSE)
b10 <- read.table("data/B10.txt", sep = "", header = FALSE)


expanded_grid <- expand.grid(c("a1", "a2", "a3", "a4", "a5", "a6", 'a7'),c("b8", "b9", "b10"))
index <- 1
for (i in list(a1, a2, a3, a4, a5, a6, a7)) {
  for (j in list(b8, b9, b10)) {
    cat("-------------------\n", expanded_grid[index, ])
    index <- index + 1
    print(wilcox.test(i[["V1"]], j[["V1"]], exact = FALSE))
  } 
}

