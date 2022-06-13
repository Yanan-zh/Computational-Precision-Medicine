library(ggplot2)
expr <- read.table("data/expr.txt", header = TRUE, sep = "\t")
pheno <- read.table("data/pheno.txt", header = TRUE, sep = "\t")
PCAplot(expr, color = paste(pheno$tissue,pheno$timepoint,sep=": "),shape = pheno$treatment_response)
