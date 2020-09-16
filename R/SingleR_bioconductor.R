#devtools::install_github('dviraran/SingleR')
library(Seurat)
library(SeuratDisk)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.1 Create Singler Object  ==========================================

# conda activate r4.0.2
library(SingleR)
library(SingleCellExperiment)
library(magrittr)
#meta_data = data.table::fread("output/20200729_meta_data.csv", showProgress = T)
#counts = data.table::fread("output/20200729_SCT_rawCounts.csv", showProgress = T)
#library(mltools)
# ====== load single cell =============
meta_data = readRDS("output/20200914_metadata.rds")
log_counts = readRDS("output/20200914_rawCounts.rds")
genes = readRDS("output/20200914_genes.rds")
genes = genes[,1]
rownames(meta_data) = colnames(log_counts)

log_counts <- as.data.frame(sapply(log_counts, as.numeric))
rownames(log_counts) = genes

sce <- SingleCellExperiment(list(logcounts=log_counts),
                                colData=DataFrame(meta_data))
# ====== load reference =============
blue_encode <- BlueprintEncodeData()
remove = grepl("CD4|CD8|Tregs|B-cells",blue_encode$label.fine)
blue_encode = blue_encode[,!remove]

immue_exp <- DatabaseImmuneCellExpressionData()

common <- Reduce(intersect, list(rownames(sce), 
                                 rownames(blue_encode),
                                 rownames(immue_exp)))
length(common)
combine_ref = do.call("cbind", list(blue_encode[common,],
                                    immue_exp[common,]))
table(combine_ref$label.fine)
system.time(trained <- trainSingleR(ref = combine_ref,
                                    labels=combine_ref$label.fine))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/20200915_singleR_pred.rds")


meta_data$barcode = rownames(meta_data)
meta_data %<>% cbind(pred)
saveRDS(x = pred, file = "output/20200915_singleR_pred.rds")

tbl <- table(pred$labels)
tbl[order(tbl,decreasing = T)]