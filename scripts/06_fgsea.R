# Load packages ----------------------------------------------------------------
library("tidyverse")
library("fgsea")


# Load data --------------------------------------------------------------------
## File for mapping probe ID to gene ID
probeID_gene <- read_tsv(file = "data/probeID_gene.tsv")

## File with logFC for differentially expressed genes
ranks_response_DEG <- read.table("data/response_geneset.txt", header = T)

## File with logFC for all genes
ranks_response_all <- read.table("data/response_geneset_all.txt", header = T)

## File with GO BP gene sets
GO_BP_genesets <- gmtPathways("data/GeneSymbols.gmt")



# Clean data -------------------------------------------------------------------

## Prepare rank files

### DEG
ranks_response_DEG <- setNames(ranks_response_DEG$logFC, 
                               ranks_response_DEG$GeneSymbol) # named vector

### All
ranks_response_all <- ranks_response_all %>% 
    select(-adj.P.Val) %>% 
    arrange(desc(logFC)) 
ranks_response_all <- setNames(ranks_response_all$logFC, 
                               ranks_response_all$GeneSymbol) # named vector


# Run fGSEA analysis -----------------------------------------------------------
fgseaRes_DEG <- fgsea(pathways = GO_BP_genesets, # list of gene sets to check
                      stats    = ranks_response_DEG, # ranked named vector
                      minSize  = 15, # minimum gene set size to include
                      maxSize  = 500) # maximum gene set size to include
fgseaRes_all <- fgsea(pathways = GO_BP_genesets, # list of gene sets to check
                      stats    = ranks_response_all, # ranked named vector
                      minSize  = 15, # minimum gene set size to include
                      maxSize  = 500) # maximum gene set size to include


## View most significant GO terms
head(fgseaRes_all[order(pval), ])



# Plot data --------------------------------------------------------------------
## Enrichment plot for a specific GO term
plotEnrichment(GO_BP_genesets[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]],
               ranks_response_all) + labs(title="GOBP_ADAPTIVE_IMMUNE_RESPONSE")


## Table plot for some selected GO terms
topPathwaysUp <- fgseaRes_all[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes_all[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(GO_BP_genesets[topPathways], ranks_response_all, fgseaRes_all, 
              gseaParam=0.5)
