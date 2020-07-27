############################################
# combine hpca and blueprint_encode
############################################
#devtools::install_github('dviraran/SingleR')
library(SingleR)
library(genefilter)
library(dplyr)
library(magrittr)
library(kableExtra)
library(GEOquery)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")
path <- paste0("./output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
####functions===========
#remove duplicate rownames with lower rowsumns
#' @param mat input as data.frame with gene name
#' @export mat matrix with gene as rownames, no duplicated genes
RemoveDup <- function(mat){
    gene_id <- as.matrix(mat[,1])
    mat <- mat[,-1]
    if(!is.matrix(mat)) mat <- sapply(mat,as.numeric)
    rownames(mat) <- 1:nrow(mat)
    mat[is.na(mat)] = 0
    mat <- cbind(mat, "rowSums" = rowSums(mat))
    mat <- mat[order(mat[,"rowSums"],decreasing = T),]
    gene_id <- gene_id[as.numeric(rownames(mat))]
    remove_index <- duplicated(gene_id)
    mat <- mat[!remove_index,]
    rownames(mat) <- gene_id[!remove_index]
    return(mat[,-ncol(mat)])
}


# 1. check and filter blueprint_encode data==============================
data("blueprint_encode")
names(blueprint_encode)
unique(blueprint_encode$types)
dim(blueprint_encode$data)
head(blueprint_encode$types)
length(blueprint_encode$types)
length(unique(blueprint_encode$types))
length(unique(blueprint_encode$main_types))
dim(blueprint_encode$data)
head(blueprint_encode$data[,1:5])
table(is.na(blueprint_encode$data))
blueprint_encode$data[is.na(blueprint_encode$data)] = 0
head(colSums(blueprint_encode$data))
#testMMM(blueprint_encode$data)
#boxplot(blueprint_encode$data, main="blueprint_encode")#slow!
#boxplot(blueprint_encode$data[,1:100])#slow!
# remove low quanlity blueprint_encode data
#par(mfrow=c(2,1))
#hist(colMeans(blueprint_encode$data),breaks=ncol(blueprint_encode$data))
quantile_75 <- apply(blueprint_encode$data,2,function(x) quantile(x,probs =0.75))
#hist(quantile_75, breaks=ncol(blueprint_encode$data))
rm_samples <- names(quantile_75)[quantile_75<1]
(rm_index <- which(colnames(blueprint_encode$data) %in% rm_samples))
blueprint_encode_rm <- blueprint_encode$data[,-rm_index]
quantile_75_new <- apply(blueprint_encode_rm,2,function(x) quantile(x,probs =0.75))
#hist(quantile_75, breaks=ncol(blueprint_encode$data))
#hist(quantile_75_new, breaks=ncol(blueprint_encode_rm),xlim = c(0,4.1))
#par(mfrow=c(1,1))
#boxplot(blueprint_encode_rm)#slow!
#title(main="blueprint_encode_rm")

# replace cell type names
Mapvalues = data.frame(c("Preadipocytes",   "Adipocytes Pre"),
                       c("CD4+ T-cells",   "T-cells CD4"),
                       c("CD4+ Tcm",       "T-cells CD4 central.memory"),
                       c("CD4+ Tem",       "T-cells CD4 effector.memory"),
                       c("CD8+ T-cells",   "T-cells CD8"),
                       c("CD8+ Tcm",       "T-cells CD8 central.memory"),
                       c("CD8+ Tem",       "T-cells CD8 effector.memory"),
                       c("Class-switched memory B-cells","B-cells class-switched.memory"),
                       c("DC",             "Dendritic-cells myeloid"),
                       c("Endothelial cells","Endothelial-cells"),
                       c("Epithelial cells","Epithelial-cells"),
                       c("Memory B-cells",  "B-cells memory"),
                       c("Mesangial cells", "Mesangial-cells"),
                       c("mv Endothelial cells","Endothelial-cells microvascular"),
                       c("naive B-cells",   "B-cells naive"),
                       c("NK cells",        "NK-cells"),
                       c("Plasma cells",    "B-cells plasma-cells"),
                       c("Skeletal muscle", "Skeletal-muscle"),
                       c("Smooth muscle",   "Smooth-muscle"),
                       c("Tregs",           "T-cells CD4 Tregs")
                       ) %>% t %>% as.data.frame(stringsAsFactors = F)
colnames(Mapvalues) = c("from","to")
blue_encode_types <- data.frame("original_names" = blueprint_encode$types[-rm_index], # didn't use main_types
                                "colnames" = colnames(blueprint_encode_rm),
                                stringsAsFactors = F)

blue_encode_types$main_types = plyr::mapvalues(blue_encode_types$original_names,
                                                  from  = Mapvalues$from, 
                                                  to = Mapvalues$to)
blue_encode_types$sub_types <- gsub("\\.[[:digit:]]+$","",blue_encode_types$colnames) %>% gsub("\\.$","",.)
# correct mistake in dataset
blue_encode_types$main_types[blue_encode_types$sub_types == "Megakaryocytes"] = "Megakaryocytes"
blue_encode_types %>% kable() %>% kable_styling()
#Blueprint_encode = CreateSinglerReference(name = 'Blueprint_encode',
#                                          expr = blueprint_encode_rm,
#                                          types = blue_encode_types$sub_types, 
#                                          main_types = blue_encode_types$main_types)
#save(Blueprint_encode,file='../SingleR/data/Blueprint_encode.RData')

# 2. check and prepare GSE107011 data==============================
# get TPM, it takes some time
getGEOSuppFiles("GSE107011",baseDir = "data/") 
GSE107011_data = read.table(gzfile("data/GSE107011/GSE107011_Processed_data_TPM.txt.gz"),
                  header = T)
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
ens <- rownames(GSE107011_data)
ensLookup <- gsub("\\.[0-9]*$", "", ens)
annotLookup <- getBM(mart=mart,
                     attributes=c("ensembl_gene_id","external_gene_name"),
                     filter="ensembl_gene_id",
                     values=ensLookup, uniqueRows=TRUE)
annotLookup <- data.frame(ens[match(annotLookup$ensembl_gene_id, ensLookup)],
                          annotLookup)
colnames(annotLookup) <- c("original_id","ensembl_gene_id","external_gene_name")
ens %<>% plyr::mapvalues(from = annotLookup$original_id, to = annotLookup$external_gene_name)
GSE107011_data = data.frame("gene_id" = ens,
                            GSE107011_data)
GSE107011_data <- RemoveDup(GSE107011_data)

# get cell types
gse <- getGEO('GSE107011') %>% .[[1]] 
GSE107011_types = gse@phenoData@data[,c("cell type:ch1","description.2")]
colnames(GSE107011_types) = c("original_names","sub_types")
GSE107011_types$colnames = colnames(GSE107011_data)
GSE107011_types$sub_types %<>% as.character() %>% substring(6)
Mapvalues = data.frame(c("Central memory CD8 T cell",  "T-cells CD8 central.memory"),
                       c("Classical monocytes",        "Monocytes classical"),
                       c("Effector memory CD8 T cells","T-cells CD8 effector.memory"),
                       c("Exhausted B cells",          "B-cells exhausted"),
                       c("Follicular helper T cells",  "T-cells follicular.helper"),
                       c("Intermediate monocytes",     "Monocytes intermediate"),
                       c("Low-density basophils",      "Basophils"),
                       c("Low-density neutrophils",    "Neutrophils"),
                       c("MAIT cells",                 "T-cells MAIT"),
                       c("Myeloid dendritic cells",    "Dendritic-cells myeloid"),
                       c("Naive B cells",              "B-cells naive"),
                       c("Naive CD4 T cells",          "T-cells CD4 naive"),
                       c("Naive CD8 T cells",          "T-cells CD8 naive"),
                       c("Natural killer cells",       "NK-cells"),
                       c("Non classical monocytes",    "Monocytes non.classical"),
                       c("Non-switched memory B cells","B-cells non.switched.memory"),
                       c("Non-Vd2 gd T cells",         "T-cells non.vd2.gd"),
                       c("PBMCs",                      "PBMCs"),
                       c("Plasmablasts",               "Plasmablasts"),
                       c("Plasmacytoid dendritic cells","Dendritic-cells plasmacytoid"),
                       c("Progenitor cells",           "Progenitor cells"),
                       c("Switched memory B cells",    "B-cells class-switched.memory"),
                       c("T regulatory cells",         "T-cells CD4 Tregs"),
                       c("Terminal effector CD4 T cells","T-cells CD4 effector.terminal"),
                       c("Terminal effector CD8 T cells","T-cells CD8 effector.terminal"),
                       c("Th1 cells",                   "T-cells Th1"),
                       c("Th1/Th17 cells",              "T-cells Th1.Th17"),
                       c("Th17 cells",                  "T-cells Th17"),
                       c("Th2 cells",                   "T-cells Th2"),
                       c("Vd2 gd T cells",              "T-cells vd2.gd")
                       ) %>% t %>% as.data.frame(stringsAsFactors = F)
colnames(Mapvalues) = c("from","to")
GSE107011_types$main_types = plyr::mapvalues(GSE107011_types$original_names,
                                               from  = Mapvalues$from, 
                                               to = Mapvalues$to)
GSE107011_types = GSE107011_types[,colnames(blue_encode_types)]
GSE107011_types %>% kable() %>% kable_styling()
blue_encode_types %>% kable() %>% kable_styling()

Annotations <- readxl::read_excel("doc/Annotations-with-abbreviations.xlsx")
# 3. merge Blueprint_encode and GSE107011 =====================
blueprint_encode_rm[1:3,1:3];GSE107011_data[1:3,1:3]


blue_encode_GSE107011 <- merge(blueprint_encode_rm, log1p(GSE107011_data),
                                   by="row.names",all=FALSE)
blue_encode_GSE107011 %<>% RemoveDup()

testMMM(blue_encode_GSE107011)
blue_encode_GSE107011[1:3,1:3]
table(Annotations$colnames %in% colnames(blue_encode_GSE107011))

jpeg(paste0(path,"boxplot_blue_encode_GSE107011.jpeg"), units="in", width=10, height=7,res=600)
par(mfrow=c(1,1))
boxplot(blue_encode_GSE107011) #slow
title(main = "boxplot for Blueprint + Encode + GSE107011")
dev.off()

# Create Singler Reference
ref = CreateSinglerReference(name = 'Blue_encode_GSE107011',
                             expr = as.matrix(blue_encode_GSE107011), # the expression matrix
                             types = c(blue_encode_types$sub_types,
                                       GSE107011_types$sub_types),
                             main_types = c(blue_encode_types$main_types,
                                            GSE107011_types$main_types))

saveRDS(ref,file='data/ref_blue_encode_GSE107011_20200609.rds')
