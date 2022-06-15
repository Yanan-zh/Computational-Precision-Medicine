# Load packages
library(tidyverse)
library(limma)
library(data.table)

# Read data
expr <- as.matrix(fread("data/expr.txt", sep = "\t"), rownames=1) %>% 
    data.frame()
pheno <- read.table("data/pheno.txt", header = TRUE, sep = "\t")
probeID_gene <- read_tsv("data/gene_expr_probe_id.tsv") %>% 
    select(gene_symbol, ID) %>% 
    rename(probeID = ID)

# DEG from paper
true_R <- read_tsv("data/test.txt", col_names = F) %>% 
    mutate(X1 = str_split(X1, " /// ")) %>% unnest()
true_NR <- read_tsv("data/test_NR.txt", col_names = F) %>% 
    mutate(X1 = str_split(X1, " /// ")) %>% unnest(X1)

# All probes + genes
probeID_gene_all <- read_tsv("data/probeID_gene.tsv") %>% 
    rename(probeID = ID)


# Divide tumor data in responders and non-responders
R_index <- which(pheno$treatment_response == "Responder (R)" &
                 pheno$tissue == "Tumor" &
                 !pheno$patient %in% c("MA", "NK", "DF"))
NR_index <- which(pheno$treatment_response == "Non-responder (NR)" &
                  pheno$tissue == "Tumor" &
                  !pheno$patient %in% c("RVS", "IP", "WR", "VJ", "AL", "DM"))
expr_R <- expr[ ,R_index]
expr_NR <- expr[ ,NR_index]
pheno_R <- pheno[R_index, ]
pheno_NR <- pheno[NR_index, ]

# Factor variable for grouping
R_timepoint <- factor(pheno_R$timepoint)
R_timepoint <- relevel(R_timepoint, ref = "pretreatment")
NR_timepoint <- factor(pheno_NR$timepoint)
NR_timepoint <- relevel(NR_timepoint, ref = "pretreatment")
R_patient <- factor(pheno_R$patient)
NR_patient <- factor(pheno_NR$patient)

# Design matrix
R_design <- model.matrix(~ 0 + R_patient + R_timepoint)
colnames(R_design) <- gsub("R_patient", "", colnames(R_design))
colnames(R_design) <- gsub("R_timepoint", "", colnames(R_design))
rownames(R_design) <- colnames(expr_R)

NR_design <- model.matrix(~ 0 + NR_patient + NR_timepoint)
colnames(NR_design) <- gsub("NR_patient", "", colnames(NR_design))
colnames(NR_design) <- gsub("NR_timepoint", "", colnames(NR_design))
rownames(NR_design) <- colnames(expr_NR)

# Differential expression analysis
fit_R <- lmFit(expr_R, R_design)
fit_R <- eBayes(fit_R)
fit_NR <- lmFit(expr_NR, NR_design)
fit_NR <- eBayes(fit_NR)

# Extract results and convert to gene names
topTable_R <- topTable(fit_R, coef = "3 weeks after start of treatment",
                       number = Inf, adjust.method = "BH") %>% 
    rownames_to_column("probeID") %>% 
    left_join(probeID_gene_all, by = "probeID") %>% 
    drop_na()
topTable_NR <- topTable(fit_NR, coef = "3 weeks after start of treatment",
                        number = Inf, adjust.method = "BH") %>% 
    rownames_to_column("probeID") %>% 
    left_join(probeID_gene, by = "probeID") %>% 
    drop_na()

# Differential expressed genes (according to logFC)
DEG_R <- topTable_R %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(gene_symbol, logFC) %>% 
    arrange(desc(logFC))
DEG_NR <- topTable_NR %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(gene_symbol, logFC) %>% 
    arrange(desc(logFC))

# Significant probes
# sig_probes_R <- topTable_R %>% 
#     filter(adj.P.Val < 0.05)
# topTable_NR %>% 
#     filter(adj.P.Val < 0.05)



# Write file
write_tsv(x = DEG_R, file = "data/DEG_R_pre_vs_post.tsv")
write_tsv(x = DEG_NR, file = "data/DEG_NR_pre_vs_post.tsv")

