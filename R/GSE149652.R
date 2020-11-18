# GSE149652
library(Seurat)
library(magrittr)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
library(data.table)
library(R.utils)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

# download https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149652
CD4_meta = fread("~/Downloads/GSE149652_CD4_TIL_droplet_cellinfo_matrice.csv.gz") %>% 
    as.data.frame
CD4_counts = fread("~/Downloads/GSE149652_CD4_TIL_droplet_count_matrice.csv.gz")
CD8_meta = fread("~/Downloads/GSE149652_CD8_TIL_droplet_cellinfo_matrice.csv.gz") %>% 
    as.data.frame
CD8_counts = fread("~/Downloads/GSE149652_CD8_TIL_droplet_count_matrice.csv.gz")

# https://www.kaggle.com/cartographic/data-table-to-sparsematrix
to_sparse <- function(d_table){
    index = d_table$index
    d_table = d_table[,-1]
    i_list <- lapply(d_table, function(x) which(x != 0))
    counts <- unlist(lapply(i_list, length), use.names = F)
    
    sparseMatrix <- sparseMatrix(
                    i = unlist(i_list, use.names = F),
                    j = rep(1:ncol(d_table), counts),
                    x = unlist(lapply(d_table, function(x) x[x != 0]), use.names = F),
                    dims = dim(d_table),
                    dimnames = list(NULL, names(d_table)))
    rownames(sparseMatrix) = index
    return(sparseMatrix)
}

system.time(CD4_counts %<>% to_sparse)
system.time(CD8_counts %<>% to_sparse)

# merge
colnames(CD4_counts) %<>% paste0("CD4_",.)
colnames(CD8_counts) %<>% paste0("CD8_",.)
CD4_meta$index %<>% paste0("CD4_",.)
CD8_meta$index %<>% paste0("CD8_",.)

system.time(counts <- RowMergeSparseMatrices(CD4_counts, CD8_counts))
meta.data = rbind(CD4_meta,CD8_meta[,colnames(CD4_meta)])
rownames(meta.data) = meta.data$index
object = CreateSeuratObject(counts = counts,
                            project = "GSE149652",
                            meta.data = meta.data)
object$labels = object$cell_types
object$labels %<>% gsub("CD4","T cells, CD4+, ",.)
object$labels %<>% gsub("CD8","T cells, CD8+, ",.)

object$main_labels = object$labels
object$main_labels %<>% gsub("T cells, CD4\\+,.*","T cells CD4+",.)
object$main_labels %<>% gsub("T cells, CD8\\+,.*","T cells CD8+",.)

saveRDS(object,file= "data/20201117_GSE149652.rds")
