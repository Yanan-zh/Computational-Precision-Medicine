# Load packages ----------------------------------------------------------------
library("tidyverse")
library("fgsea")
library("writexl")

# Load data --------------------------------------------------------------------

## File with logFC for responders pre- and post-treatment
ranks_R_pre_vs_post <- read.table("results/DEG_R_pre_vs_post_purity.tsv", header = T)

## File with logFC for non-responders pre- and post-treatment
ranks_NR_pre_vs_post <- read.table("results/DEG_NR_pre_vs_post_purity.tsv", header = T)

## File with logFC for pre-treatment responders vs non-responders
ranks_pre_R_vs_NR <- read.table("results/DEG_pre_R_vs_NR_purity.tsv", header = T)

## File with GO BP gene sets
GO_BP_genesets <- gmtPathways("data/GeneSymbols.gmt")



# Clean data -------------------------------------------------------------------

# Ranks for responders pre- and post-treatment
ranks_R_pre_vs_post <- setNames(ranks_R_pre_vs_post$logFC, 
                                ranks_R_pre_vs_post$gene_symbol)
# Ranks for non-responders pre- and post-treatment
ranks_NR_pre_vs_post <- setNames(ranks_NR_pre_vs_post$logFC, 
                                 ranks_NR_pre_vs_post$gene_symbol)
# Ranks for pre-treatment responders vs non-responders 
ranks_pre_R_vs_NR <- setNames(ranks_pre_R_vs_NR$logFC, 
                              ranks_pre_R_vs_NR$gene_symbol)


# Run fGSEA analysis -----------------------------------------------------------

fgseaRes_R_pre_vs_post <- fgsea(pathways = GO_BP_genesets, # list of gene sets to check
                                stats    = ranks_R_pre_vs_post, # ranked named vector
                                minSize  = 15, # minimum gene set size to include
                                maxSize  = 500) 

fgseaRes_NR_pre_vs_post <- fgsea(pathways = GO_BP_genesets, # list of gene sets to check
                                 stats    = ranks_NR_pre_vs_post, # ranked named vector
                                 minSize  = 15, # minimum gene set size to include
                                 maxSize  = 500) 

fgseaRes_pre_R_vs_NR <- fgsea(pathways = GO_BP_genesets, # list of gene sets to check
                              stats    = ranks_pre_R_vs_NR, # ranked named vector
                              minSize  = 15, # minimum gene set size to include
                              maxSize  = 500) # maximum gene set size to include

## View most significant GO terms
head(fgseaRes_R_pre_vs_post[order(pval), ], 10)
head(fgseaRes_NR_pre_vs_post[order(pval), ], 10)
head(fgseaRes_pre_R_vs_NR[order(pval), ], 10)


# Plot data --------------------------------------------------------------------
## Enrichment plot for a specific GO term 
plotEnrichment(GO_BP_genesets[["GOBP_BLOOD_VESSEL_MORPHOGENESIS"]],
               ranks_R_pre_vs_post) + labs(title="GOBP_BLOOD_VESSEL_MORPHOGENESIS")


## Table plot of top up- and down-regulated GO BP
topUpGO_R_pre_post <- fgseaRes_R_pre_vs_post[ES > 0][head(order(pval), n=10), pathway]
topDownGO_R_pre_post <- fgseaRes_R_pre_vs_post[ES < 0][head(order(pval), n=10), pathway]
topGO_R_pre_post <- c(topUpGO_R_pre_post, rev(topDownGO_R_pre_post))
plotGseaTable(GO_BP_genesets[topGO_R_pre_post], 
              ranks_R_pre_vs_post, 
              fgseaRes_R_pre_vs_post, 
              gseaParam=0.5)

## Table plot of top up- and down regulated GO BP
topUpGO_NR_pre_post <- fgseaRes_NR_pre_vs_post[ES > 0][head(order(pval), n=10), pathway]
topDownGO_NR_pre_post <- fgseaRes_NR_pre_vs_post[ES < 0][head(order(pval), n=10), pathway]
topGO_NR_pre_post <- c(topUpGO_NR_pre_post, rev(topDownGO_NR_pre_post))
plotGseaTable(GO_BP_genesets[topGO_NR_pre_post], 
              ranks_NR_pre_vs_post, 
              fgseaRes_NR_pre_vs_post, 
              gseaParam=0.5)

## Table plot of top up- and down regulated GO BP
topUpGO_pre_R_NR <- fgseaRes_pre_R_vs_NR[ES > 0][head(order(pval), n=10), pathway]
topDownGO_pre_R_NR <- fgseaRes_pre_R_vs_NR[ES < 0][head(order(pval), n=10), pathway]
topGO_pre_R_NR <- c(topUpGO_pre_R_NR, rev(topDownGO_pre_R_NR))
plotGseaTable(GO_BP_genesets[topGO_pre_R_NR], 
              ranks_pre_R_vs_NR, 
              fgseaRes_pre_R_vs_NR, 
              gseaParam=0.5)


# Write files ------------------------------------------------------------------

## To tsv
write_tsv(x = fgseaRes_R_pre_vs_post, file = "results/fgsea_R_pre_vs_post.tsv")
write_tsv(x = fgseaRes_NR_pre_vs_post, file = "results/fgseaRes_NR_pre_vs_post.tsv")
write_tsv(x = fgseaRes_pre_R_vs_NR, file = "results/fgseaRes_pre_R_vs_NR.tsv")

## To xlsx
fgseaRes_R_pre_vs_post <- fgseaRes_R_pre_vs_post %>% 
    rowwise() %>% mutate(leadingEdge = paste(unlist(leadingEdge), collapse = ','))
fgseaRes_NR_pre_vs_post <- fgseaRes_NR_pre_vs_post %>% 
    rowwise() %>% mutate(leadingEdge = paste(unlist(leadingEdge), collapse = ','))
fgseaRes_pre_R_vs_NR <- fgseaRes_pre_R_vs_NR %>% 
    rowwise() %>% mutate(leadingEdge = paste(unlist(leadingEdge), collapse = ','))

write_xlsx(x = fgseaRes_R_pre_vs_post, path = "results/fgsea_R_pre_vs_post.xlsx")
write_xlsx(x = fgseaRes_NR_pre_vs_post, path = "results/fgsea_NR_pre_vs_post.xlsx")
write_xlsx(x = fgseaRes_pre_R_vs_NR, path = "results/fgsea_pre_R_vs_NR.xlsx")
