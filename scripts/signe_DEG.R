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

# All probes + genes
probeID_gene_all <- read_tsv("data/probeID_gene.tsv") %>% 
    rename(probeID = ID)

# Tumor purity
purity <- read_tsv(file = "data/tumor_purity.tsv") %>% 
    select(TumorPurity)

############################################################################
###                   Responder pre- vs post-treatment                   ###
###                                and                                   ###
###                Non-responder pre- vs post-treatment                  ###
############################################################################


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
purity_R <- purity[R_index, ] %>% pull()
purity_NR <- purity[NR_index, ] %>% pull()


# Factor variable for grouping
R_timepoint <- factor(pheno_R$timepoint)
R_timepoint <- relevel(R_timepoint, ref = "pretreatment")
NR_timepoint <- factor(pheno_NR$timepoint)
NR_timepoint <- relevel(NR_timepoint, ref = "pretreatment")
R_patient <- factor(pheno_R$patient)
NR_patient <- factor(pheno_NR$patient)


# Design matrix
R_design <- model.matrix(~ R_patient + R_timepoint + purity_R)
colnames(R_design) <- gsub("R_patient", "", colnames(R_design))
colnames(R_design) <- gsub("R_timepoint", "", colnames(R_design))
rownames(R_design) <- colnames(expr_R)

NR_design <- model.matrix(~ NR_patient + NR_timepoint + purity_NR)
colnames(NR_design) <- gsub("NR_patient", "", colnames(NR_design))
colnames(NR_design) <- gsub("NR_timepoint", "", colnames(NR_design))
rownames(NR_design) <- colnames(expr_NR)


# Differential expression analysis
fit_R <- lmFit(combat_expr_R, R_design)
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
    left_join(probeID_gene_all, by = "probeID") %>% 
    drop_na()


# Differential expressed genes (according to logFC)
DEG_R <- topTable_R %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(gene_symbol, logFC) %>% 
    group_by(gene_symbol) %>%
    slice(which.max(abs(logFC))) %>% # select gene with largest absolute logFC
    arrange(desc(logFC))
DEG_NR <- topTable_NR %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(gene_symbol, logFC) %>% 
    group_by(gene_symbol) %>%
    slice(which.max(abs(logFC))) %>%  
    arrange(desc(logFC))

# Significant genes
topTable_R %>% 
    filter(adj.P.Val < 0.05) %>% 
    group_by(gene_symbol) %>% 
    slice(which.min(adj.P.Val))
# 4
topTable_NR %>% 
    filter(adj.P.Val < 0.05)
# 0


# Write file
write_tsv(x = DEG_R, file = "results/DEG_R_pre_vs_post_purity.tsv")
write_tsv(x = DEG_NR, file = "results/DEG_NR_pre_vs_post_purity.tsv")



############################################################################
###              Pre-treatment responder vs non-responder                ###
############################################################################


# Extract tumor pre-treatment data
pre_index <- which(pheno$tissue == "Tumor" &
                   pheno$timepoint == "pretreatment") 
expr_pre <- expr[ ,pre_index]
pheno_pre <- pheno[pre_index, ]
purity_pre <- purity[pre_index, ] %>% pull()
    

# Factor variable for grouping
pre_treatment_response <- factor(pheno_pre$treatment_response)
pre_treatment_response <- relevel(pre_treatment_response, ref = "Non-responder (NR)")
pre_batch <- factor(pheno_pre$batch)


# Design matrix
pre_design <- model.matrix(~ pre_treatment_response + purity_pre + pre_batch)
colnames(pre_design) <- gsub("pre_treatment_response", "", colnames(pre_design))
colnames(pre_design) <- gsub("pre_batch", "", colnames(pre_design))
rownames(pre_design) <- colnames(expr_pre)


# Differential expression analysis
fit_pre <- lmFit(expr_pre, pre_design)
fit_pre <- eBayes(fit_pre)

# Extract results and convert to gene names
topTable_pre <- topTable(fit_pre, coef = "Responder (R)",
                         number = Inf, adjust.method = "BH") %>% 
    rownames_to_column("probeID") %>% 
    left_join(probeID_gene_all, by = "probeID") %>% 
    drop_na()

# Differential expressed genes (according to logFC)
DEG_pre <- topTable_pre %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(gene_symbol, logFC) %>% 
    group_by(gene_symbol) %>% 
    slice(which.max(abs(logFC))) %>% # select gene with largest absolute logFC
    arrange(desc(logFC))

# Significant genes
topTable_pre %>% 
     filter(adj.P.Val < 0.05)
# 0


# Write file
write_tsv(x = DEG_pre, file = "results/DEG_pre_R_vs_NR_purity.tsv")
