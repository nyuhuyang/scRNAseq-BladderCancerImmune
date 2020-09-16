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

source("R/util.R")
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
as.data.frame(table(singlerDF$labels))
T_NK <- grepl("NK|T cells", singlerDF$labels)
table(grepl("NK|T cells", singlerDF$labels))
T_NK <- grepl("NK|T cells", singlerDF$labels)
T_NK_df <- singlerDF[T_NK,]
T_NK_summary <- as.data.frame(table(T_NK_df$labels))
colnames(T_NK_summary)= c("cell_type", "cell_number")
T_NK_summary[,"cell_type_cell_number"] = paste(T_NK_summary$cell_type,":", T_NK_summary$cell_number)

options(repr.plot.width=10, repr.plot.height=7)

ggpie(T_NK_summary, "cell_number", label = "cell_type_cell_number",
      fill = "cell_type", color = "white",
      palette = singler.colors)