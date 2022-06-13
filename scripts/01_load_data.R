library(tidyverse)
tab <- read.table("data/GSE60331_series_matrix.txt", sep = "\t", skip = 31, fill = TRUE)
div <- which(tab[,1] == "!series_matrix_table_begin") # This line divides expr and pheno



# Expression data
expr <- tab[(div+1):(nrow(tab)-1),]
colnames(expr) <- expr[1,]
rownames(expr) <- expr[,1]
expr <- expr[-1,-1]

# Meta data
pheno <- tab[1:(div-1),] %>% 
    filter(V1 == "!Sample_characteristics_ch1") %>% 
    select(-V1) %>% t %>% 
    as.data.frame()
colnames(pheno) <- apply(pheno, 2, function(x) {str_extract(x,"(.+)(?=:)")[1]})
pheno <- pheno %>% 
    mutate(across(.fns = ~str_extract(.x, "(?<=:\\s)(.+)"))) %>% 
    rename(treatment_response = "treatment response")

# Microarray platform data

probes_genes <- read_tsv("data/GPL15207-17536.txt",
                         skip = 36) %>% 
    dplyr::select(ID,`Gene Symbol`,`Entrez Gene`) %>% 
    mutate(across(.fns = ~str_split(.x, 
                                    pattern = "\\/\\/\\/"))) %>% 
    unnest(cols = everything()) %>% 
    mutate(across(.fns = str_trim)) %>% 
    rename(c(gene_symbol = `Gene Symbol`,
             entrez_gene = `Entrez Gene`))


# Save data
write.table(expr, "data/expr.txt", sep = "\t", row.names = TRUE, quote = FALSE)
write.table(pheno, "data/pheno.txt", sep = "\t", row.names = FALSE, quote = TRUE)
write_tsv(x = probes_genes, "data/probeID_gene.tsv")
