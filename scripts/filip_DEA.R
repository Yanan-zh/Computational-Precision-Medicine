library(tidyverse)
library(limma)

# load the data 
pheno <- read_delim("data/pheno.txt")

gene_symbols <- read_tsv(file = "data/probeID_gene.tsv") %>% 
    select(ID, gene_symbol)

pheno_trans <- pheno %>% 
    mutate(age = str_extract(age, 
                             pattern = "[\\d]+") %>% 
               as.numeric()) %>% 
    mutate(treatment_response = str_extract(treatment_response,
                                            pattern = "(?<=\\().+(?=\\))")) %>% 
    mutate(timepoint = case_when(str_detect(timepoint,
                                            pattern = "(pre)") ~ "T0",
                                 T ~ "T1")) %>% 
    mutate(response_timepoint = str_c(treatment_response, 
                                      timepoint,
                                      sep = "_") %>% 
               as_factor()) %>% 
    filter(tissue == "Tumor") 

i <- which(pheno$tissue == "Tumor")

purity <- read_tsv(file = "data/tumor_purity.tsv") %>% 
    select(TumorPurity)
purity <- purity[i, ] %>% pull()


MCP <- read.table(file = "data/MCPcounter.txt",
                  sep = "\t", 
                  header = T) %>% 
    t() %>% 
    as_tibble()

MCP <- MCP[i, ]
CD8 <- MCP %>% pull(`CD8 T cells`)
B_lin <- MCP %>% pull(`B lineage`)


### diff expression on gene_level

expr_gene <- read.table("data/gene_expression.tsv",
                   header = T,
                   row.names = 1) 
expr_gene <- expr_gene[, i]

expr_probe <- read.table("data/expr.txt",
                         header = T)
expr_probe <- expr_probe[, i]




patient_1 <- as.factor(pheno_trans$patient)
response_timepoint <- relevel(pheno_trans$response_timepoint,
                              ref = "NR_T0")

### Design for multilevel model
design <- model.matrix(~0 + response_timepoint + purity + B_lin + CD8)
design %>% view()

keep <- str_detect(colnames(design),"response_timepoint")
# keep_2 <- str_detect(colnames(design), "patient")

colnames(design)[keep] <- str_extract(colnames(design)[keep], 
                                      "(?<=response_timepoint)(.+)")
# colnames(design)[1] <- "NR_T0"

# colnames(design)[keep_2] <- str_extract(colnames(design)[keep_2],
#                                         "(?<=patient_)(.+)")


### within patient correlation of gene expression model
# corfit_gene <- duplicateCorrelation(expr_gene,
#                                design,
#                                block = pheno_trans$patient,
#                                )
# ### within patient correlation of probe expression model
# corfit_probe <- duplicateCorrelation(expr_probe,
#                                design,
#                                block = pheno_trans$patient)
# 
?lmFit

corfit_gene <- duplicateCorrelation(expr_gene,
                                    design = design,
                                    block = pheno_trans$patient)
corfit_probe <- duplicateCorrelation(expr_gene,
                                     design = design,
                                     block = pheno_trans$patient)

fit_gene <- lmFit(expr_gene, 
             design = design,
             block = pheno_trans$patient,
             correlation = corfit_gene$consensus)

fit_probe <- lmFit(expr_probe,
                   design = design,
                   block = pheno_trans$patient,
                   correlation = corfit_probe$consensus)

# contrasts 
cont <- makeContrasts(
    NR_pre_treat = NR_T1-NR_T0,
    R_pre_treat = R_T1-R_T0,
    Treat_NR_R = R_T1-NR_T1,
    Pre_NR_R = R_T0-NR_T0,
    levels = design
)


fit_2_gene <- contrasts.fit(fit_gene,
                            contrasts = cont)
fit_2_gene <- eBayes(fit_2_gene)
fit_2_probe <- contrasts.fit(fit_probe,
                             contrasts = cont)
fit_2_probe <- eBayes(fit_2_probe)



### Responder pre and after treatment
top_table_gene_1 <- topTable(fit_2_gene, 
                             coef = "R_pre_treat",
                             adjust.method = "BH",
                             number = Inf)
fit_2_gene$correlation

R_pre_treat_gene <- top_table_gene_1 %>% 
    as_tibble(rownames = "gene_symbol") %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(gene_symbol, logFC)

top_table_probe_1 <- topTable(fit_2_probe,
                              coef = "R_pre_treat",
                              adjust.method = "BH",
                              number = Inf)
R_pre_treat_probe <- top_table_probe_1 %>% 
    as_tibble(rownames = "ID") %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(ID, logFC)
### Non-Responder pre and after treatment
top_table_gene_2 <- topTable(fit_2_gene, 
                               coef = "NR_pre_treat",
                               adjust.method = "BH",
                               number = Inf)

NR_pre_treat_gene <- top_table_gene_2 %>% 
    as_tibble(rownames = "gene_symbol") %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(gene_symbol, logFC)
top_table_probe_2 <- topTable(fit_2_probe, 
                               coef = "NR_pre_treat",
                               adjust.method = "BH",
                               number = Inf)
NR_pre_treat_probe <- top_table_probe_2 %>% 
    as_tibble(rownames = "ID") %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(ID, logFC)
### Difference in non versus responders after treatment
top_table_gene_3 <- topTable(fit_2_gene, 
                                coef = "Treat_NR_R",
                                adjust.method = "BH",
                                number = Inf)
Treat_NR_R_gene <- top_table_gene_3 %>% 
    as_tibble(rownames = "gene_symbol") %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(gene_symbol, logFC)

top_table_probe_3 <- topTable(fit_2_probe, 
                                coef = "Treat_NR_R",
                                adjust.method = "BH",
                                number = Inf)
Treat_NR_R_probe <- top_table_probe_3 %>% 
    as_tibble(rownames = "ID") %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(ID, logFC)
### Difference in non versus responders before treatment
top_table_gene_4 <- topTable(fit_2_gene, 
                               coef = "Pre_NR_R",
                               adjust.method = "BH",
                               number = Inf)
Pre_NR_R_gene <- top_table_gene_4 %>% 
    as_tibble(rownames = "gene_symbol") %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(gene_symbol, logFC)

top_table_probe_4 <- topTable(fit_2_probe, 
                                coef = "Pre_NR_R",
                                adjust.method = "BH",
                                number = Inf)

Pre_NR_R_probe <- top_table_probe_4 %>% 
    as_tibble(rownames = "ID") %>% 
    filter(logFC > 1 | logFC < -1) %>% 
    select(ID, logFC)



### results


result_tables_text <- c(
        "R_pre_treat_gene",
        "R_pre_treat_probe",
        "NR_pre_treat_gene",
        "NR_pre_treat_probe",
        "Treat_NR_R_gene",
        "Treat_NR_R_probe",
        "Pre_NR_R_gene",
        "Pre_NR_R_probe"
    )

result_tables <- list(
    R_pre_treat_gene,
    R_pre_treat_probe,
    NR_pre_treat_gene,
    NR_pre_treat_probe,
    Treat_NR_R_gene,
    Treat_NR_R_probe,
    Pre_NR_R_gene,
    Pre_NR_R_probe
)

signi_tables <- list(
    top_table_gene_1,
    top_table_probe_1,
    top_table_gene_2,
    top_table_probe_2,
    top_table_gene_3,
    top_table_probe_3,
    top_table_gene_4,
    top_table_probe_4
)
## save significant genes and probes tables 
walk2(.x = signi_tables,
      .y = result_tables_text,
      .f = ~{
                if (str_detect(.y, 
                             pattern = "(probe)")) {
                  .x %>% 
                  as_tibble(rownames = "ID") %>% 
                      left_join(y = gene_symbols,
                                by = "ID") %>% 
                        relocate(gene_symbol) %>% 
                      filter(adj.P.Val < 0.05) %>% 
                        arrange(desc(logFC)) %>% 
                      write_tsv(., 
                                file = 
                                    str_c("data/tables/",
                                          .y,
                                          "_signi.tsv")) }
          else {
          .x %>% as_tibble(rownames = "gene") %>%
              filter(adj.P.Val < 0.05) %>% 
                  arrange(desc(logFC)) %>% 
              write_tsv(., 
                        file = 
                            str_c("data/tables/",
                                  .y,
                                  "_signi.tsv")) }
                  
      })

# save logFC sorted tables
walk2(.x = result_tables,
       .y = result_tables_text,
       .f = ~{
           if (str_detect(.y, 
                          pattern = "(probe)")) {
               .x %>%
                   left_join(y = gene_symbols,
                             by = "ID") %>% 
                   relocate(gene_symbol) %>%   
                   arrange(desc(logFC)) %>% 
           write_tsv(x = ., 
                    file = str_c("data/tables/", .y ,".tsv"))
               
               }
           else {
               .x %>% arrange(desc(logFC)) %>% 
               write_tsv(x = ., 
                         file = str_c("data/tables/", .y ,".tsv"))
           }
           })
