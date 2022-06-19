library(tidyverse)
tab <- read.table("data/GSE60331_series_matrix.txt", sep = "\t", skip = 31, fill = TRUE)
div <- which(tab[,1] == "!series_matrix_table_begin") # This line divides expr and pheno
#microarray_platform <- read.table("data/GPL15207-17536.txt", sep="\t", header = T, fill = T)

# Index of duplicate samples
duplicates <- pheno %>% 
    duplicated() %>% 
    which()

# Expression data
expr <- tab[(div+1):(nrow(tab)-1),]
colnames(expr) <- expr[1,]
rownames(expr) <- expr[,1]
expr <- expr[-1,-1] %>% select(-duplicates)

# Meta data
pheno <- tab[1:(div-1),] %>% 
    filter(V1 == "!Sample_characteristics_ch1") %>% 
    select(-V1) %>% t %>% 
    as.data.frame()
colnames(pheno) <- apply(pheno, 2, function(x) {str_extract(x,"(.+)(?=:)")[1]})
pheno <- pheno %>% 
    mutate(across(.fns = ~str_extract(.x, "(?<=:\\s)(.+)"))) %>% 
    rename(treatment_response = "treatment response") %>% 
    slice(-duplicates)


probes_genes <- read_tsv("data/GPL15207-17536.txt",
                         skip = 36) %>%
    dplyr::select(ID,`Gene Symbol`,`Entrez Gene`) %>%
    write_tsv(file = "data/check.tsv") %>%
    mutate(across(.fns = ~str_split(.x,
                                    pattern = "\\/\\/\\/"))) %>%
    unnest(cols = everything()) %>%
    mutate(across(.fns = str_trim)) %>%
    rename(c(gene_symbol = `Gene Symbol`,
             entrez_gene = `Entrez Gene`)) %>% 
    filter(!is.na(gene_symbol)) %>% 
    filter(!(gene_symbol == "---"))


# Save data
write.table(expr, "data/expr.txt", sep = "\t", row.names = TRUE, quote = FALSE)
write.table(pheno, "data/pheno.txt", sep = "\t", row.names = FALSE, quote = TRUE)
#write.table(probeID_GeneSymbol, "data/probeID_GeneSymbol.txt", sep = "\t", row.names = F)
write_tsv(x = probes_genes, "data/probeID_gene.tsv")
