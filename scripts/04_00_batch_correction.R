BiocManager::install("sva")
library(sva)


expr = read.table('data/gene_expression.tsv',header = TRUE,row.names = 1)
pheno = read.table('data/pheno.txt',header = TRUE)



pheno = pheno %>% mutate(age = str_remove_all(age, "y")) %>% unite(type, tissue:timepoint) %>% mutate(type = gsub("_pretreatment", "", type)) %>% mutate(type = gsub("eeks after start of", "", type))

pheno$age = sapply(pheno$age, as.numeric)

#__________________________________________________________________

# batch correction
#__________________________________________________________________

batch = pheno$batch

batch_expr = ComBat(dat=expr, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

## Recalculate PCA
pca_batch <- prcomp(t(batch_expr))

## Check out how much variance is explained
summary(pca_batch)

df <- data.frame(PC1 = pca_batch$x[,1], PC2 = pca_batch$x[,2])
df$batch <- as.factor(pheno$batch)
df$shape <- pheno$type
# show batch corrected PCA
ggplot(df, aes(x = PC1, y = PC2, color=batch, shape = shape )) +
  geom_point()


#write out batch corrected
write.table(batch_expr, "data/batch_expr.txt", sep = "\t", row.names = TRUE, quote = FALSE)
