library(tidyverse)
library(data.table)

### Function for PCA ###

PCAplot <- function(dat, color, shape = NULL,legend=TRUE,plot.title=""){
    # perform a PCA on the data
    pca <- prcomp(t(dat))
    
    # the contribution to the total variance for each component
    percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
    
    # assembly the data for the plot
    d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Color=color, name=colnames(dat))
    if(is.null(shape)){
        Aes <- aes_string(x="PC1", y="PC2",color= "color")
    }else{
        d<-data.frame(d,Shape=shape)
        Aes <- aes_string(x="PC1", y="PC2",color="color", shape = "shape")
    }
    legend.pos <- if(legend) {"right"} else {"none"}
    
    ggplot(data=d) +
        Aes +
        geom_point(size=3) +
        labs(x = paste0("PC1: ",round(percentVar[1] * 100),"% variance"),
             y = paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
        ggtitle(plot.title) +
        theme(legend.position = legend.pos)
}

### Load data ###

expr <- as.matrix(fread("data/expr.txt", sep = "\t"),rownames=1) %>% 
    data.frame()
pheno <- read.table("data/pheno.txt", header = TRUE, sep = "\t")


# PCA for mucosal and tumor pre-treatment samples
index_M_T0 <- which(pheno$timepoint == "pretreatment") 
expr_M_T0 <- expr[ ,index_M_T0]
pheno_M_t0 <- pheno[index_M_T0, ]

PCAplot(expr_M_T0, 
        color = paste(pheno_M_t0$tissue,pheno_M_t0$timepoint,sep=": "),
        shape = pheno_M_t0$treatment_response)

# PCA for responders
index_R <- which(pheno$treatment_response == "Responder (R)") 
expr_R <- expr[ ,index_R]
pheno_R <- pheno[index_R, ]

PCAplot(expr_R, 
        shape = paste(pheno_R$tissue,pheno_R$timepoint,sep=": "),
        color = pheno_R$patient)

# PCA for all samples
PCAplot(expr, 
        color = paste(pheno$tissue,pheno$timepoint,sep=": "),
        shape = pheno$treatment_response) +
    labs(shape = "Treatment Response",
         color = "Tissue and Timepoint",
         title = "Tissue and Responders") +
    theme_minimal() +
    ggthemes::scale_color_colorblind()


# PCA showing batches
PCAplot(expr,
        color = pheno$batch) +
    labs(color = "Batch", 
         title = "Batch Effect") + 
    theme_minimal() +
    ggthemes::scale_color_colorblind()

