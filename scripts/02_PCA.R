library(ggplot2)
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

# DO IT
expr <- read.table("data/expr.txt", header = TRUE, sep = "\t")
pheno <- read.table("data/pheno.txt", header = TRUE, sep = "\t")
PCAplot(expr, color = paste(pheno$tissue,pheno$timepoint,sep=": "),shape = pheno$treatment_response)
