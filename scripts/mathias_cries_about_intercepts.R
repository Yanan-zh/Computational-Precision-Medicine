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

# Extract tumor pre-treatment data
pre_index <- which(pheno$tissue == "Tumor" &
                       pheno$timepoint == "pretreatment") 
expr_pre <- expr[ ,pre_index]
pheno_pre <- pheno[pre_index, ]
purity_pre <- purity[pre_index, ] %>% pull()


# Factor variable for grouping
pre_treatment_response <- factor(pheno_pre$treatment_response)
pre_treatment_response <- relevel(pre_treatment_response, ref = "Non-responder (NR)")
patients <- pheno$patient[pre_index]
dat <- fread("data/GSE60331_series_matrix.txt", fill = TRUE)
batch <- dat[57,] %>% str_extract('(?<=\"\")(\\d+)') %>% as.integer() %>% .[-1]
batch <- batch[!str_detect(dat[32,], "rep2")[-1]][pre_index]

# Design matrix
pre_design <- model.matrix(~pre_treatment_response)

# Differential expression analysis
fit_pre <- lmFit(expr_pre, pre_design)
fit_pre <- eBayes(fit_pre)

#### Find DEGs ####
topTable(fit_pre, 
         coef=2, 
         adjust="BH", 
         number = Inf, 
         p.value=0.05,
         lfc = 1)
#####


fit_pre$coefficients[1,]

dx <- pre_design[,2] %>% as.factor()
NR <- pre_design[,1] %>% as.integer()
R <- pre_design[,2] %>% as.integer()

dy <- expr_pre[1,] %>% as.numeric()
ggplot(mapping = aes(x = dx,
                     y = dy)) +
    geom_point() +
    geom_hline(aes(yintercept = 4.281548*NR+4.160127*R,
                   color = dx)) +
    ggtitle("No Intercept")
