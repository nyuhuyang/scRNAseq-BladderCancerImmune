library(Seurat)
library(DT)
library(ggplot2)
library(magrittr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path))dir.create(path, recursive = T)

object = readRDS("data/20200920_TCGA-BLCA.rds")
BLCA.prop_main = readRDS("data/20200929_TCGA-BLCA.prop_main.rds")
table(rownames(object@meta.data) == rownames(BLCA.prop_main$Est.prop.weighted))
object@meta.data %<>% cbind(BLCA.prop_main$Est.prop.weighted)
object$HSC_level = factor(object$HSC > 0 , labels =c("HSC% ==0", "HSC% >0"))

table(object$HSC_level)
death = c("days_to_death",
          "subtype_Days.until.death",
          "subtype_Combined.days.to.last.followup.or.death")[3]
Remove = is.na(BLCA_seurat@meta.data[,death]) |
    BLCA_seurat@meta.data[,death] %in% c("[Not Available]","[Discrepancy]","NA") |
    BLCA_seurat@meta.data[,death] < 0

table(!Remove)
object = object[,!Remove]

Idents(object) = "HSC_level"
Idents(object) %<>% factor(levels = c("HSC% >0","HSC% ==0"))
table(Idents(object))
DEG <- FindAllMarkers.UMI(object, logfc.threshold = 0.1,test.use = "MAST",
                       only.pos = TRUE,latent.vars = "nFeature_RNA",
                       p.adjust.methods = "fdr")
DEG$gene %<>% paste0(" ",.)
colnames(DEG) %<>% gsub("UMI","FPKM",.)
colnames(DEG)[8] ="group"
write.csv(DEG, file = paste0(path,"BLCA_HSC_markers.csv"))
