library(TCGAWorkflowData)
library(TCGAbiolinks)

library(Seurat)
library(DT)
library(SummarizedExperiment)
library(ggplot2)
library(scater)
library(MuSiC)
library(xbioc)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
# expression
query.exp <- GDCquery(project = "TCGA-BLCA",
                      data.category = "Gene expression",
                      data.type = "Gene expression quantification",
                      platform = "Illumina HiSeq", 
                      file.type  = "results", 
                      #sample.type = "Primary solid Tumor",
                      legacy = TRUE)

GDCdownload(query.exp)
exp <- GDCprepare(query.exp, save = FALSE)
exp
assays(exp)[[1]][1:4,1:4]
genes = gsub("\\|.*","",rownames(assays(exp)$raw_count))
raw_count <- assays(exp)$raw_count[!duplicated(genes),]
rownames(raw_count) = gsub("\\|.*","",rownames(raw_count))

# 2 Process TCGA data and visualization
BLCA_seurat <- CreateSeuratObject(counts = raw_count, project = "TCGA-BLCA",
                                  meta.data = as.data.frame(colData(exp))).
BLCA_seurat <- NormalizeData(BLCA_seurat)
BLCA_seurat <- FindVariableFeatures(BLCA_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BLCA_seurat)
BLCA_seurat <- ScaleData(BLCA_seurat, features = all.genes)

# 2.1 Perform linear dimensional reduction
BLCA_seurat <- RunPCA(BLCA_seurat, features = VariableFeatures(object = BLCA_seurat))
DimPlot(BLCA_seurat, reduction = "pca",pt.size = 3)
DimHeatmap(BLCA_seurat, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(BLCA_seurat)

# 2.2 Clustering and Visualization
BLCA_seurat <- FindNeighbors(BLCA_seurat, dims = 1:20)
BLCA_seurat <- FindClusters(BLCA_seurat, resolution = 0.5)

BLCA_seurat <- RunUMAP(BLCA_seurat, dims = 1:20)
DimPlot(BLCA_seurat, reduction = "umap",pt.size = 3)
DimPlot(BLCA_seurat, reduction = "umap", group.by = "tissue_or_organ_of_origin",
        pt.size = 3)
FeaturePlot(BLCA_seurat, features = "subtype_Neoantigen.load",
            pt.size = 3)
saveRDS(BLCA_seurat, file = "data/20200920_TCGA-BLCA.rds")

# 3 Deconvolution by MuSiC
# read BladderCancerImmune data
BLCA_seurat = readRDS("data/20200920_TCGA-BLCA.rds")
# convert TCGA-BLCA counts to a ExpressionSet object
BLCA.sce <- new("ExpressionSet", exprs = as.matrix(BLCA_seurat[["RNA"]]@counts),
                phenoData = AnnotatedDataFrame(BLCA_seurat@meta.data))
# read BladderCancerImmune data
scBLCA_seurat = readRDS("data/20200914_BladderCancerImmune.rds")
scBLCA.sce <- new("ExpressionSet", exprs = as.matrix(scBLCA_seurat[["RNA"]]@counts),
                  phenoData = AnnotatedDataFrame(scBLCA_seurat@meta.data))
table(scBLCA.sce$labels)
BLCA.prop = music_prop(bulk.eset = BLCA.sce, sc.eset = scBLCA.sce,
                       clusters = "labels",
                       samples = 'orig.ident', 
                       verbose = T)
saveRDS(BLCA.prop, file = "data/20200929_TCGA-BLCA.prop.rds")

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
BLCA.prop.main = music_prop(bulk.eset = BLCA.sce, sc.eset = scBLCA.sce,
                       clusters = "main_labels",
                       samples = 'orig.ident', 
                       verbose = T)

# Jitter plot of estimated Epi cell type proportions
jitter.fig = Jitter_Est(list(data.matrix(BLCA.prop.main$Est.prop.weighted),
                             data.matrix(BLCA.prop.main$Est.prop.allgene)),
                        method.name = c('MuSiC', 'NNLS'), title = 'Jitter plot of Est Proportions')
jpeg(paste0(path,"MuSiC_NNLS_main_celltype.jpeg"), units="in", width=10, height=7,res=600)
print(jitter.fig)
dev.off()

DoHeatmap.matrix(data.use = t(BLCA.prop.main$Est.prop.weighted),
                 save.path = "output/20200929/Est.prop.weighted", do.print = T)