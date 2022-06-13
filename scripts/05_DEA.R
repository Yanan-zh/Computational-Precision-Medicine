#### Differential Expression Analysis ####
library(tidyverse)
library(limma)

# Read data
expr <- read.table("data/expr.txt", header = TRUE, sep = "\t")
pheno <- read.table("data/pheno.txt", header = TRUE, sep = "\t")

group <- apply(pheno,1,function(x) {
    if (x[4] == "Non-responder (NR)"){
        a <- "NR"
    }else{
        a <- "R"
    }
    
    if (x[5] == "Tumor"){
        if (x[6] == "pretreatment"){
            b <- "T1"
        }else{
            b <- "W3"
        }
    }else{
        b <- "M"
    }
    
    str_c(b,a,sep="_")
})

pheno <- pheno %>% 
    mutate(group = group)

design <- model.matrix(~0+group+patient, data = pheno)
keep <- str_detect(colnames(design),"group")
colnames(design)[keep] <- str_extract(colnames(design)[keep], "(?<=group)(.+)")
fit <- lmFit(expr, design)

contrast.matrix <- makeContrasts((M_R+T1_R+W3_R)/3-(M_NR+T1_NR+W3_NR)/3,
                                 (M_R+M_NR)/2-(T1_R+T1_NR+W3_R+W3_NR)/4,
                                 levels = design)
colnames(contrast.matrix) <- c("R_vs_NR", "Muc_vs_Tumor")
fit2 <-  contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

topTable(fit2, coef=2, adjust="BH") %>% View()
results <- decideTests(fit2)
vennDiagram(results)

