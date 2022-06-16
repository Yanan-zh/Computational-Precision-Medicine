# Load packages ----------------------------------------------------------------
library("tidyverse")
library("fgsea")


# Load data --------------------------------------------------------------------

## File with logFC for all genes
#ranks_response_all <- read.table("data/response_geneset_all.txt", header = T)

## File with logFC for pre-treatment responder vs non-responder
#ranks_pre_R_vs_NR <- read.table("data/response_geneset_pre.tsv", header = T)

## File with logFC for responders pre- and post-treatment
ranks_R_pre_vs_post <- read.table("data/DEG_R_pre_vs_post.tsv", header = T)

## File with logFC for non-responders pre- and post-treatment
ranks_NR_pre_vs_post <- read.table("data/DEG_NR_pre_vs_post.tsv", header = T)

## File with logFC for pre-treatment responders vs non-responders
ranks_pre_R_vs_NR <- read.table("data/DEG_pre_R_vs_NR.tsv", header = T)

## File with GO BP gene sets
GO_BP_genesets <- gmtPathways("data/GeneSymbols.gmt")
#GO_entrez <- gmtPathways("data/Entrez.gmt")



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
                              maxSize  = 400) 

fgseaRes_NR_pre_vs_post <- fgsea(pathways = GO_BP_genesets, # list of gene sets to check
                                 stats    = ranks_NR_pre_vs_post, # ranked named vector
                                 minSize  = 15, # minimum gene set size to include
                                 maxSize  = 400) 

fgseaRes_pre_R_vs_NR <- fgsea(pathways = GO_BP_genesets, # list of gene sets to check
                              stats    = ranks_pre_R_vs_NR, # ranked named vector
                              minSize  = 15, # minimum gene set size to include
                              maxSize  = 400) # maximum gene set size to include

## View most significant GO terms
head(fgseaRes_R_pre_vs_post[order(pval), ], 10)
head(fgseaRes_NR_pre_vs_post[order(pval), ], 10)
head(fgseaRes_pre_R_vs_NR[order(pval), ], 10)


# Plot data --------------------------------------------------------------------
## Enrichment plot for a specific GO term for responders
plotEnrichment(GO_BP_genesets[["GOBP_CELL_MIGRATION"]],
               ranks_pre_R_vs_NR) + labs(title="GOBP_CELL_MIGRATION")


## Table plot for some selected GO terms for responders
topPathwaysUp <- fgseaRes_pre_R_vs_NR[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes_pre_R_vs_NR[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(GO_BP_genesets[topPathways], ranks_pre_R_vs_NR, fgseaRes_pre_R_vs_NR, 
              gseaParam=0.5)

# Write files ------------------------------------------------------------------
write_tsv(x = fgseaRes_R_pre_vs_post, file = "results/fgsea_R_pre_vs_post.tsv")
write_tsv(x = fgseaRes_NR_pre_vs_post, file = "results/fgsea_NR_pre_vs_post.tsv")
write_tsv(x = fgseaRes_pre_R_vs_NR, file = "results/fgsea_pre_R_vs_NR.tsv")
