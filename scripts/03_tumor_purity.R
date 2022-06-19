library(estimate)
library(tidyverse)

# load probe level data
start <- Sys.time()

gene_expr_tib <- read_tsv("data/gene_expression.tsv")

gene_expr <- as.data.frame(gene_expr_tib %>% 
                               select(!gene_symbol))

rownames(gene_expr) <- gene_expr_tib %>% pull(gene_symbol)


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
