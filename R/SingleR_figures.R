library(SingleR)
library(Seurat)
library(SeuratDisk)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================
# conda activate r4.0.2
Convert("output/20200914_BladderCancerImmune.h5ad", 
        dest = "output/20200914_BladderCancerImmune.h5seurat", overwrite = F)
object <- LoadH5Seurat("output/20200914_BladderCancerImmune.h5seurat")


##############################
# create singleR data frame
###############################
singlerDF = readRDS("output/20200915_singleR_pred.rds")
as.data.frame(table(singlerDF$labels)) %>% kable %>% kable_styling()

singlerDF$main_labels = singlerDF$labels
singlerDF$main_labels %<>% gsub("T cells, CD4\\+","CD4T",.)
singlerDF$main_labels %<>% gsub("T cells, CD8\\+","CD8T",.)
singlerDF$main_labels %<>% gsub(",.*","",.)

main_labels_summary <- as.data.frame(table(singlerDF$main_labels))
main_labels_summary$prec = main_labels_summary$Freq/sum(main_labels_summary$Freq)*100 
main_labels_summary$prec %<>% round(digits = 1) %>% as.character %>% paste0("%")


colnames(main_labels_summary)[1:2]= c("cell_type", "cell_number")
main_labels_summary[,"abb_cell_type"] = substr(main_labels_summary$cell_type, 1, 4) %>% gsub(" .*","",.)
main_labels_summary[,"cell_type_per"] = paste0(main_labels_summary$abb_cell_type,":", main_labels_summary$prec)

jpeg(paste0(path,"main_labels_summary.jpeg"), units="in", width=10, height=7,res=300)
ggpie(main_labels_summary, "cell_number", label = "cell_type_per",
      fill = "cell_type", color = "white",lab.pos = "out",
      palette = Singler.colors)
dev.off()

T_NK <- grepl("Macrophages|B cells|NK|T cells|Monocytes", singlerDF$labels)
table(grepl("Macrophages|B cells|NK|T cells|Monocytes", singlerDF$labels))
T_NK_df <- singlerDF[T_NK,]
T_NK_df$labels %<>% gsub("T cells, CD4\\+","CD4T",.)
T_NK_df$labels %<>% gsub("T cells, CD8\\+","CD8T",.)
#T_NK_df$labels %<>% gsub(",.*","",.)

T_NK_summary <- as.data.frame(table(T_NK_df$labels))
T_NK_summary$prec = T_NK_summary$Freq/sum(T_NK_summary$Freq)*100 
T_NK_summary$prec %<>% round(digits = 1) %>% as.character %>% paste0("%")
colnames(T_NK_summary)[1:2]= c("cell_type", "cell_number")
T_NK_summary[,"abb_cell_type"] = substr(T_NK_summary$cell_type, 1, 4) %>% gsub(" .*","",.)
T_NK_summary[,"cell_type_per"] = paste0(T_NK_summary$abb_cell_type,":", T_NK_summary$prec)

jpeg(paste0(path,"NK_T_summary.jpeg"), units="in", width=10, height=7,res=300)
ggpie(T_NK_summary, "cell_number", label = "cell_type_per",
      fill = "cell_type", color = "white",lab.pos = "out",
      palette = Singler.colors)
dev.off()

# GSE149652
object = readRDS(file= "data/20201117_GSE149652.rds")

labels_summary <- as.data.frame(table(object$labels))
labels_summary$prec = labels_summary$Freq/sum(labels_summary$Freq)*100 
labels_summary$prec %<>% round(digits = 1) %>% as.character %>% paste0("%")
labels_summary$Var1 %<>% gsub("T cells, CD4\\+","4",.)
labels_summary$Var1 %<>% gsub("T cells, CD8\\+","8",.)

colnames(labels_summary)[1:2]= c("cell_type", "cell_number")
labels_summary[,"cell_type_per"] = paste0(labels_summary$cell_type,":", labels_summary$prec)

jpeg(paste0(path,"GSE149652_labels_summary.jpeg"), units="in", width=10, height=7,res=300)
ggpie(labels_summary, "cell_number", label = "cell_type_per",
      fill = "cell_type", color = "white",lab.pos = "out",
      palette = Singler.colors)
dev.off()
