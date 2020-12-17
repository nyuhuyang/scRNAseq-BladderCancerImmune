library(Seurat)
library(DT)
library(tibble)
library(SummarizedExperiment)
library(ggplot2)
library(TCGAbiolinks)
library(immunedeconv)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

#==========================================
query.exp.BLCA <- GDCquery(project = "TCGA-BLCA",
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type  = "results", 
                      #sample.type = "Primary solid Tumor",
                      legacy = TRUE)

GDCdownload(query.exp.BLCA)
exp.BLCA <- GDCprepare(query.exp.BLCA, save = FALSE)
exp.BLCA[1:4,1:4]
genes = gsub("\\|.*","",rownames(assays(exp.BLCA)$raw_count))
FPKM.BLCA <- assays(exp.BLCA)$raw_count[!duplicated(genes),]
rownames(FPKM.BLCA) = gsub("\\|.*","",rownames(FPKM.BLCA))

FPKM.BLCA %<>% as.data.frame %>% rownames_to_column
colnames(FPKM.BLCA)[1] = "gene_symbol"
FPKM.BLCA[1:4,1:4];class(FPKM.BLCA)
data.table::fwrite(FPKM.BLCA,file =paste0(path,"TCGA-BLCA_FPKM.txt"),
                   sep = "\t",row.names = FALSE, col.names = TRUE)


#=======================
# expression
query.exp.BLCA.LAML <- GDCquery(project = c("TCGA-BLCA","TCGA-LAML"),
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type  = "results", 
                      #sample.type = "Primary solid Tumor",
                      legacy = TRUE)

GDCdownload(query.exp.BLCA.LAML)
exp.BLCA.LAML <- GDCprepare(query.exp.BLCA.LAML, save = FALSE)
assays(exp.BLCA.LAML)[[1]][1:4,1:4]
genes = gsub("\\|.*","",rownames(assays(exp.BLCA.LAML)$raw_count))
FPKM.BLCA.LAML <- assays(exp.BLCA.LAML)$raw_count[!duplicated(genes),]
rownames(FPKM.BLCA.LAML) = gsub("\\|.*","",rownames(FPKM.BLCA.LAML))

FPKM.BLCA.LAML %<>% as.data.frame %>% rownames_to_column
colnames(FPKM.BLCA.LAML)[1] = "gene_symbol"
FPKM.BLCA.LAML[1:4,1:4];class(FPKM.BLCA.LAML)
data.table::fwrite(FPKM.BLCA.LAML,file =paste0(path,"TCGA-BLCA_LAML_FPKM.txt"),
                   sep = "\t",row.names = FALSE, col.names = TRUE)
FPKM.BLCA.LAML = read.delim(file ="output/20201215/TCGA-BLCA_LAML_FPKM.txt")
# =========== cibersort=====================

URL <- "https://raw.githubusercontent.com/veltenlab/rnamagnet/master/R/CIBERSORT.R"
download.file(URL, destfile = "R/archive/CIBERSORT.R", method="curl")
# download LM22 from https://cibersort.stanford.edu/download.php
set_cibersort_binary("R/archive/CIBERSORT.R")
set_cibersort_mat("data/LM22.txt")

cibersort = immunedeconv::deconvolute(gene_expression = FPKM, method = "cibersort_abs")
cibersort.BLCA = read.table("output/20201215/CIBERSORTx_Job7_Results.txt",sep = "\t",
                       row.names = 1,header = TRUE)
table(rowSums(cibersort.BLCA[,1:22]))
cibersort.BLCA[,20:25] %>% kable %>% kable_styling()

cibersort.BLCA.LAML = read.table("output/20201215/CIBERSORTx_Job8_Results.txt",sep = "\t",
                            row.names = 1,header = TRUE)
table(rowSums(cibersort.BLCA.LAML[,1:22]))
#============ Convert to Seurat ===================
colnames(FPKM.BLCA)[1] = "rowname"
FPKM.BLCA %<>% column_to_rownames()
BLCA_seurat <- CreateSeuratObject(counts = FPKM.BLCA, project = "TCGA-BLCA",
                                  meta.data = as.data.frame(colData(exp.BLCA)))
BLCA_seurat@meta.data %<>% cbind(cibersort.BLCA)
table(BLCA_seurat$definition, BLCA_seurat$P.value < 0.104)

colnames(FPKM.BLCA.LAML)[1] = "rowname"
FPKM.BLCA.LAML %<>% column_to_rownames()

BLCA_LAML_seurat <- CreateSeuratObject(counts = FPKM.BLCA.LAML, project = "TCGA-BLCA_LAML",
                                  meta.data = as.data.frame(colData(exp.BLCA.LAML)))
BLCA_LAML_seurat = readRDS("data/20201010_TCGA-BLCA&LAML.rds")
BLCA_LAML_seurat@meta.data %<>% cbind(cibersort.BLCA.LAML)

meta_data = BLCA_seurat@meta.data
Cell_types = colnames(cibersort.BLCA)[1:22]
meta_data = meta_data[meta_data$P.value < 0.104 ,]
calculate_survive_curve(meta_data = meta_data, Cell_types, day = "days_to_death", 
                        save.path = paste0(path,"P_0.104_","days_to_death/"))

Remove = is.na(meta_data[,"paper_Combined.days.to.last.followup.or.death"]) |
    meta_data[,"paper_Combined.days.to.last.followup.or.death"] %in% c("[Not Available]","[Discrepancy]") |
    meta_data[,"paper_Combined.days.to.last.followup.or.death"] < 0
meta_data = meta_data[!Remove,]
meta_data$paper_Combined.days.to.last.followup.or.death %<>% as.integer()
calculate_survive_curve(meta_data = meta_data, Cell_types, day = "paper_Combined.days.to.last.followup.or.death", 
                        save.path = paste0(path,"P_0.104_","paper_Combined.days.to.last.followup.or.death/"))


ggscatter(x = "P.value", y = "RMSE",data = meta_data) + 
    geom_vline(xintercept=0.05, linetype="dashed",
               color = "red", size=0.5)+
    geom_vline(xintercept=0.104, linetype="dashed",
               color = "blue", size=0.5)



VlnPlot(BLCA_LAML_seurat, features = c("P.value","Correlation","RMSE"),group.by = "orig.ident")
meta_data = BLCA_LAML_seurat@meta.data
meta_data = meta_data[meta_data$orig.ident %in% "TCGA-LAML",]









days = c("sample","definition","days_to_last_follow_up","days_to_death","subtype_Combined.days.to.last.followup.or.death") 
sapply(days,function(x) head(meta.data[,x])) %>% kable %>% kable_styling()

dup_sample = meta.data$sample[duplicated(meta.data$sample)]
dup_sample1 = meta.data[meta.data$sample %in% dup_sample[1],]

Idents(BLCA_seurat) = "definition"
Normal = meta.data[meta.data$definition %in% "Solid Tissue Normal",]
sapply(days,function(x) head(Normal[,x])) %>% kable %>% kable_styling()
BLCA_seurat %<>% subset(idents = "Primary solid Tumor")

BLCA_seurat@meta.data[,days[5]] %<>% as.numeric()
sapply(days[3:5],function(x) table(BLCA_seurat@meta.data[,x] >= 0 )) 

hist(BLCA_seurat@meta.data[,days[3]], main = paste("Histogram of",days[3]))
hist(BLCA_seurat@meta.data[,days[4]], main = paste("Histogram of",days[4]))
BLCA_seurat@meta.data[,days[5]] %<>% as.numeric()
hist(BLCA_seurat@meta.data[,days[5]], main = paste("Histogram of",days[5]))

Remove = is.na(BLCA_seurat@meta.data[,"subtype_Combined.days.to.last.followup.or.death"]) |
    BLCA_seurat@meta.data[,"subtype_Combined.days.to.last.followup.or.death"] %in% c("[Not Available]","[Discrepancy]") |
    BLCA_seurat@meta.data[,"subtype_Combined.days.to.last.followup.or.death"] < 0
BLCA_seurat = BLCA_seurat[,!Remove]

hist(BLCA_seurat@meta.data[,days[5]], main = paste("Histogram of",days[5]))
