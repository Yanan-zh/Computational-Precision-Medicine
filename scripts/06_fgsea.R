# Load packages ----------------------------------------------------------------
library("tidyverse")
library("fgsea")


# Load data --------------------------------------------------------------------
expr <- read.table("data/expr.txt", header = TRUE, sep = "\t")


## File with ranks: cols = c(probeID, log2FC)
ranks <- read.table(rnk.file,
                    header = TRUE, 
                    colClasses = c("character", "numeric"))
ranks <- setNames(ranks$t, ranks$ID) # Convert to named numeric vector
str(ranks)

## File with GO BP gene sets
GO_BP_genesets <- gmtPathways("GeneSymbols.gmt")
str(head(pathways))

## File for mapping probe ID to gene ID
probeID_GeneSymbol <- read.table(file = "data/probeID_GeneSymbol.txt",
                                 header = T)


# Clean data -------------------------------------------------------------------

## Map probe IDs to gene IDs
# probes <- ranks$probes %>% as_tibble %>% rename(probeID = value)
probes <- as.vector(rownames(expr)[1:30]) %>% as_tibble %>% rename(probeID = value)
left_join(probes, probeID_GeneSymbol,
          by = "probeID")

## Ensure all gene names of pathways are present in ranks - else remove them 



# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# probe_to_entrez <- getBM(attributes = c('affy_hg_u133_plus_2', 'entrezgene_id'),
#                          filters = 'affy_hg_u133_plus_2',
#                          values = probes, 
#                          mart = ensembl)



# Run fGSEA analysis -----------------------------------------------------------
fgseaRes <- fgsea(pathways = examplePathways, # list of gene sets to check
                  stats    = exampleRanks, # ranked named vector
                  minSize  = 15, # minimum gene set size to include
                  maxSize  = 500) # maximum gene set size to include

## View most significant GO terms
head(fgseaRes[order(pval), ])


# Plot data --------------------------------------------------------------------
## Enrichment plot for a specific GO term
plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
               exampleRanks) + labs(title="Programmed Cell Death")


## Table plot for some selected GO terms
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam=0.5)
