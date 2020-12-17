library(Seurat)
library(DT)
library(SummarizedExperiment)
library(ggplot2)
library(scater)
library(MuSiC)
library(xbioc)
library(Biobase)
library(BisqueRNA)
library(TCGAbiolinks)
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


#=======================
# expression
query.exp <- GDCquery(project = c("TCGA-BLCA","TCGA-LAML"),
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

raw_count <- as.data.frame(raw_count) %>% rownames_to_column
colnames(raw_count)[1] = "gene_symbol"
raw_count[1:4,1:4];class(raw_count)
data.table::fwrite(raw_count,file =paste0(path,"TCGA-BLCA_LAML_FPKM.txt"),
                   sep = "\t",row.names = FALSE, col.names = TRUE)

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
saveRDS(BLCA_seurat, file = "data/20201010_TCGA-BLCA&LAML.rds")