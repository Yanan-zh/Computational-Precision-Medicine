library(estimate)
library(tidyverse)

# load probe level data
start <- Sys.time()

expr <- read.table("data/expr.txt",
                   header = TRUE)

probes_genes <- read_tsv("data/probeID_gene.tsv") %>% 
    select(!entrez_gene)


# gene_names <- read_tsv(file = "data/gene_names.tsv") %>% 
#     mutate(name_abbr = str_extract(string = Name,
#         pattern = "(?<=\\()[\\w\\-]+?(?=\\))")) 

# write_tsv(gene_names, file = "data/gene_names.tsv")

# rownames(expr) <- gene_names

# associate genes with probes
# 
# probes_genes <- read_tsv("data/GPL15207-17536.txt",
#                          skip = 36) %>% 
#     dplyr::select(ID,`Gene Symbol`,`Entrez Gene`) %>% 
#     mutate(across(.fns = ~str_split(.x, 
#                                     pattern = "\\/\\/\\/"))) %>% 
#     unnest(cols = everything()) %>% 
#     mutate(across(.fns = str_trim)) %>% 
#     rename(c(gene_symbol = `Gene Symbol`,
#              entrez_gene = `Entrez Gene`))
# 
# write_tsv(x = probes_genes,
#           "data/probeID_gene.tsv")

# merge gene symbols with probe expression data

expr_tib <- expr %>% 
    as_tibble(rownames = "ID") %>% 
    left_join(y = probes_genes) %>% 
    relocate(gene_symbol) %>% 
    filter(!is.na(gene_symbol)) %>% 
    filter(!(gene_symbol == "---"))

# prepare data for collapsing rows 

#rowID <- rownames(expr)

# genes <- probes_genes %>%
#     select(ID,gene_symbol) %>%
#     distinct(ID, .keep_all = T) %>% 
#     pull()
## how to choose which probes to use? 

# sanity check
# length(rowID)
# length(genes)

# create gene level expression data
# gene_level_expr <- WGCNA::collapseRows(
#     datET = expr,
#     rowGroup = genes,
#     rowID = rowID,
# )

expr_tib <- expr_tib %>% 
    group_by(gene_symbol) %>% 
    summarise(across(.cols = !ID, 
                     .fns = ~mean(.x, na.rm = T),
                     .names = "{.col}"))


gene_expr <- as.data.frame(expr_tib %>% 
                               select(!gene_symbol))

rownames(gene_expr) <- expr_tib %>% pull(gene_symbol)


#gene_expr <- as.data.frame(gene_level_expr$datETcollapsed)

write.table(x = gene_expr, 
            file = "data/gene_expr.txt", 
            sep = "\t",
            quote = F)

#Intersect genes with common genes from ESTIMATE as per protocol

filterCommonGenes(input.f = "data/gene_expr.txt",
                  output.f = "data/gene_expr_comm.gct")


# read.table("data/gene_expr.txt", header = TRUE, row.names = 1, 
#            sep = "\t", quote = "", stringsAsFactors = FALSE)

#Create ESTIMATE scores
estimateScore("data/gene_expr_comm.gct",
              output.ds = "data/estimate_scores.gct",
              platform = "affymetrix")


ES <- read_delim(file = "data/estimate_scores.gct",
                 skip = 2) 
# transposed data to fit pheno
ES_2 <- ES %>% 
    select(where(is.numeric)) %>% 
    as.matrix() %>% 
    t() %>% 
    as_tibble() 

colnames(ES_2) <- ES %>% 
    pull(1)


# Some plots of distribution of scores 
tumor_density <- ES_2 %>% 
    ggplot(mapping = aes(x = TumorPurity)) +
    geom_density(fill = "gray") 
    
ggsave(filename = "images/03_purity_density.jpg",
       plot = tumor_density,
       device = "jpeg",
       width = 20,
       height = 10
       )

score_densities <- ES_2 %>% 
    select(!TumorPurity) %>% 
    pivot_longer(cols = everything(),
                 names_to = "score_type",
                 values_to = "score") %>% 
    ggplot(aes(x = score,
               fill = score_type)) +
    geom_density(alpha = 0.5) +
    labs(fill = element_blank()) 
    
ggsave(filename = "images/03_score_densities.jpg",
       plot = score_densities,
       device = "jpeg",
       width = 20,
       height = 10)

write_tsv(x = ES_2,
          file = "data/tumor_purity.tsv")

end <- Sys.time()

print(end - start)