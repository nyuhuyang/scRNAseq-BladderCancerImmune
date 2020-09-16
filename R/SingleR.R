#devtools::install_github('dviraran/SingleR')
library(Seurat)
library(SeuratDisk)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================
# conda activate r3.6.2
Convert(paste0("output/20200728_SCT.h5ad"), 
        dest = paste0("output/20200728_SCT.h5seurat"), overwrite = F)
object <- LoadH5Seurat("output/20200728_SCT.h5seurat")

# conda activate r4.0
library(SingleR)
object <- LoadH5Seurat("output/SCT.h5seurat")

object %<>% subset(subset = n_genes >= 100)
object_data <- object[["RNA"]]@data
remove(object); GC()
ref = readRDS(file='data/ref_blue_encode_GSE107011_20200609.rds')

singler <- CreateBigSingleRObject.1(object_data, annot = NULL, project.name="Lung",
                N = 5000, min.genes = 100, technology = "10X",
                species = "Human", citation = "", ref.list = list(ref),
                normalize.gene.length = F, variable.genes = "de", fine.tune = T,
                reduce.file.size = F, do.signatures = F, do.main.types = T,
                temp.dir = getwd(), numCores = SingleR.numCores)
saveRDS(singler,file="output/singlerT_BladderCancer_20200722.rds")
