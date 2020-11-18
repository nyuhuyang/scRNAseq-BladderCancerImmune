# conda activate r4.0.2
library(Seurat)
library(DT)
library(SummarizedExperiment)
library(ggplot2)
library(scater)
library(MuSiC)
library(xbioc)
library(Biobase)
library(BisqueRNA)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# 3 Deconvolution by MuSiC
# read BladderCancerImmune data
BLCA_seurat = readRDS("data/20200920_TCGA-BLCA.rds")
BLCA_LAML_seurat = readRDS("data/20201010_TCGA-BLCA&LAML.rds")

# convert TCGA-BLCA counts to a ExpressionSet object
    
BLCA.sce <- new("ExpressionSet", exprs = as.matrix(BLCA_seurat[["RNA"]]@counts),
                phenoData = AnnotatedDataFrame(BLCA_seurat@meta.data))
BLCA_LAML.sce <- new("ExpressionSet", exprs = as.matrix(BLCA_LAML_seurat[["RNA"]]@counts),
                phenoData = AnnotatedDataFrame(BLCA_LAML_seurat@meta.data))
# read BladderCancerImmune data
scBLCA_seurat = readRDS("data/20200914_BladderCancerImmune.rds")
scT_seurat = readRDS("data/20201117_GSE149652.rds")

scBLCA.sce <- new("ExpressionSet", exprs = as.matrix(scBLCA_seurat[["RNA"]]@counts),
                  phenoData = AnnotatedDataFrame(scBLCA_seurat@meta.data))
scT.sce <- new("ExpressionSet", exprs = as.matrix(scT_seurat[["RNA"]]@counts),
                  phenoData = AnnotatedDataFrame(scT_seurat@meta.data))
scT.sce$orig.ident = scT.sce$Sample_ID
common <- Reduce(intersect, list(rownames(scBLCA.sce), 
                                 rownames(scT.sce),
                                 #rownames(BLCA.sce),
                                 rownames(BLCA_LAML.sce)))
length(common)
scBLCA_T.sce = do.call("cbind", list(scBLCA.sce[common,],
                                    scT.sce[common,]))

table(scBLCA_T.sce$labels)
table(scBLCA_T.sce$orig.ident)
BLCA.prop = music_prop(bulk.eset = BLCA.sce, sc.eset = scBLCA_T.sce,
                       clusters = "labels",
                       samples = 'orig.ident', 
                       verbose = T)
saveRDS(BLCA.prop, file = "data/20201118_TCGA-BLCA.prop.rds")

str(BLCA.prop)
names(BLCA.prop)

# Jitter plot of estimated all cell type proportions
cell.labels = sort(colnames(BLCA.prop$Est.prop.weighted))
BLCA.prop$Est.prop.weighted = BLCA.prop$Est.prop.weighted[,cell.labels]
jitter.fig = Jitter_Est(list(data.matrix(BLCA.prop$Est.prop.weighted),
                             data.matrix(BLCA.prop$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jpeg(paste0(path,"MuSiC_NNLS_Sub_celltype.jpeg"), units="in", width=10, height=7,res=600)
print(jitter.fig)
dev.off()

# Estimate major cell type proportions
scBLCA_T.sce$main_labels %<>% gsub("T cells,", "T cells",. )
table(scBLCA_T.sce$main_labels)

BLCA.prop.main = music_prop(bulk.eset = BLCA.sce, sc.eset = scBLCA.sce,
                       clusters = "main_labels",
                       samples = 'orig.ident', 
                       verbose = T)
saveRDS(BLCA.prop.main, file = "data/20201118_TCGA-BLCA.prop_main.rds")

# Jitter plot of estimated Epi cell type proportions
jitter.fig = Jitter_Est(list(data.matrix(BLCA.prop.main$Est.prop.weighted),
                             data.matrix(BLCA.prop.main$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jpeg(paste0(path,"MuSiC_NNLS_main_celltype.jpeg"), units="in", width=10, height=7,res=600)
print(jitter.fig)
dev.off()

DoHeatmap.matrix(data.use = t(BLCA.prop.main$Est.prop.weighted),
                 save.path = "output/20200929/Est.prop.weighted", do.print = T)



BLCA_LAML.prop = music_prop(bulk.eset = BLCA_LAML.sce, sc.eset = scBLCA_T.sce,
                       clusters = "labels",
                       samples = 'orig.ident', 
                       verbose = T)
saveRDS(BLCA_LAML.prop, file = "data/20201118_BLCA_LAML.prop.rds")

str(BLCA_LAML.prop)
names(BLCA_LAML.prop)
