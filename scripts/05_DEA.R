#### Differential Expression Analysis ####
library(tidyverse)
library(limma)

# Read data
expr <- read.table("data/expr.txt", header = TRUE, sep = "\t")
# expr <- read.table("data/GeneSymbol_expr.txt", header = TRUE, sep = "\t")
pheno <- read.table("data/pheno.txt", header = TRUE, sep = "\t") %>% 
    mutate(group = apply(pheno,1,function(x) {
        if (x[4] == "Non-responder (NR)"){a <- "NR"}
        else {a <- "R"}
        if (x[5] == "Tumor"){
            if (x[6] == "pretreatment"){b <- "T1"}
            else{b <- "W3"}
        }else{b <- "M"}
        str_c(b,a,sep="_")
    }))

design <- model.matrix(~0+group+patient, data = pheno)
keep <- str_detect(colnames(design),"group")
colnames(design)[keep] <- str_extract(colnames(design)[keep], "(?<=group)(.+)")
fit <- lmFit(expr, design)

contrast.matrix <- makeContrasts((M_R+T1_R+W3_R)/3-(M_NR+T1_NR+W3_NR)/3,
                                 (M_R+M_NR)/2-(T1_R+T1_NR+W3_R+W3_NR)/4,
                                 levels = design)
colnames(contrast.matrix) <- c("R_vs_NR", 
                               "Muc_vs_Tumor")
fit2 <-  contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tab <- topTable(fit2, coef=2, adjust="BH")
results <- decideTests(fit2)
vennDiagram(results)

decideTests(fit2, adjust.method = "BH",
            p.value = 0.05,
            lfc = 1) %>% 
    as.data.frame() %>% 
    filter(R_vs_NR != 0,
           Muc_vs_Tumor == 0) %>% 
    rownames()


geneset <- topTable(fit2, 
         coef=1, 
         adjust="BH", 
         number = nrow(expr), 
         p.value=0.05) %>% 
    select(logFC) %>%
    rownames_to_column("probeID") %>% 
    filter(!probeID %in% (topTable(fit2, 
                                coef=2, 
                                adjust="BH", 
                                number = nrow(expr), 
                                p.value=0.05) %>% rownames())) %>% 
    arrange(desc(logFC)) %>% 
    inner_join(read.table("data/probeID_GeneSymbol.txt", 
                          header = TRUE)) %>% 
    select(-probeID) %>% 
    arrange(GeneSymbol) %>% 
    distinct(GeneSymbol, .keep_all = TRUE) %>% 
    relocate(GeneSymbol) %>% 
    arrange(desc(logFC))

write.table(geneset, "data/response_geneset.txt", 
            row.names = FALSE)

