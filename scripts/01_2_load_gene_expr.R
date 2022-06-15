probes_genes <- read_tsv("data/probeID_gene.tsv") %>% 
    select(!entrez_gene) %>%
    filter(!is.na(gene_symbol)) %>% 
    filter(!(gene_symbol == "---")) %>% 
    distinct(.keep_all = T) 

expr <- read.table("data/expr.txt",
                              header = TRUE)

expr_tib <- expr %>% 
    as_tibble(rownames = "ID") %>% 
    left_join(y = probes_genes) %>% 
    relocate(gene_symbol) %>% 
    filter(!is.na(gene_symbol))

#probes that map to multiple genes
expr_tib %>% 
    select(ID,gene_symbol) %>% 
    distinct(gene_symbol, .keep_all = T) %>% 
    add_count(ID) %>% 
    filter(n > 1) %>% 
    distinct(ID)

#gene-level expression
gene_expr <- 
    expr_tib  %>% 
    mutate(m = rowMeans(across(where(is.numeric)))) %>% 
    group_by(gene_symbol) %>% 
    filter(m == max(m)) %>% 
    relocate(m) %>% 
    ungroup() %>% 
    distinct(gene_symbol, .keep_all = T) 


### genes with multiple genes mapping to one probe
gene_expr %>%  
    select(gene_symbol, ID) %>%
    distinct(gene_symbol, .keep_all = T) %>% 
    add_count(ID) %>% 
    filter(n >1) %>% 
    distinct(ID)

### genes with multiple probes mapping to one gene
gene_expr %>%  
    select(gene_symbol, ID, m) %>%
    distinct(ID, .keep_all = T) %>% 
    add_count(gene_symbol) %>% 
    filter(n >1) %>% 
    distinct(gene_symbol)

# genes <- expr_tib %>% pull(gene_symbol)
# probes <- expr_tib %>% pull(ID)
    
# WGCNA::collapseRows(datET = expr_matrix,
#                     rowGroup = genes,
#                     rowID = probes,
#                     method = "MaxMean"
#                     )

write_tsv(gene_expr %>% 
              select(!m & !ID), file = "data/gene_expression.tsv")


