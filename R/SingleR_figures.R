library(SingleR)
library(Seurat)
library(SeuratDisk)
library(magrittr)
library(pheatmap)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

source("R/util.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================
# Step 1: Spearman coefficient
(load(file = "data/MCL_41_harmony_20191231.Rda"))
singler = readRDS("output/singlerT_BladderCancer_20200722.rds")
(load(file="output/singlerT_MCL_41_20200225.Rda"))


##############################
# create singleR data frame
###############################
singlerDF = data.frame("singler1sub" = singler$singler[[1]]$SingleR.single$labels,
                       "sub_types" = singler$singler[[1]]$SingleR.single.main$labels,
                       row.names = rownames(singler$singler[[1]]$SingleR.single$labels),
                       stringsAsFactors = F)
singlerDF$main_types = gsub("\\ .*","", singlerDF$sub_types)
singlerDF = singlerDF[colnames(object),]
table(rownames(singlerDF) == colnames(object))

##############################


##############################
# check the spearman correlation
###############################
#Or by all cell types (showing the top 50 cell types):
jpeg(paste0(path,"DrawHeatmap_sub1.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n = 50,normalize = F))
dev.off()
jpeg(paste0(path,"DrawHeatmap_sub1_N.jpeg"), units="in", width=10, height=7,
     res=600)
print(SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single,top.n = 50,normalize = T))
dev.off()

#Finally, we can also view the labeling as a table compared to the original identities:

table(singlerDF$singler1sub, singlerDF$orig.ident) %>% kable %>%
        kable_styling()
##############################
# adjust cell label
##############################
# combine cell types
singlerDF$singler1sub = gsub("MCL:.*","MCL",singlerDF$singler1sub)
singlerDF$singler1sub = gsub("B_cells:PB","B_cells:Plasma_cells",singlerDF$singler1sub)
singlerDF$cell.types = gsub("B_cells:.*","B_cells",singlerDF$singler1sub)
singlerDF$cell.types = gsub("MEP|CLP|HSC|CMP|GMP|MPP","HSC/progenitors",singlerDF$cell.types)
singlerDF$cell.types = gsub("T_cells:CD4\\+_.*","T_cells:CD4+",singlerDF$cell.types)
singlerDF$cell.types = gsub("T_cells:CD8\\+_.*","T_cells:CD8+",singlerDF$cell.types)
singlerDF$cell.types = gsub("T_cells:Tregs","T_cells:CD4+",singlerDF$cell.types)
singlerDF$cell.types = gsub("DC|Macrophages|Macrophages:M1","Myeloid cells",singlerDF$cell.types)
singlerDF$cell.types = gsub("Erythrocytes","Myeloid cells",singlerDF$cell.types)
singlerDF$cell.types = gsub("Eosinophils|Megakaryocytes|Monocytes","Myeloid cells",singlerDF$cell.types)
singlerDF$cell.types = gsub("Adipocytes|Fibroblasts|mv_Endothelial_cells","Nonhematopoietic cells",singlerDF$cell.types)
table(singlerDF$cell.types, singlerDF$orig.ident) %>% kable() %>% kable_styling()

# reduce false positive results (B cells are labeled as MCL in normal samples)
# and false negative results (MCL cells are labeled as B cells in MCL samples)
# singler1sub false negative results  =========
CCND1 = FetchData(object,"CCND1")
singlerDF$CCND1 = CCND1$CCND1
singlerDF[(singlerDF$CCND1 >0 & singlerDF$cell.types %in% "B_cells"),"cell.types"] = "MCL"
# cell.types false positive results  ========
table(singlerDF$cell.types, object@meta.data$orig.ident) %>% kable %>% kable_styling()
normal_cells <- object$sample %in% c("BH","DJ","MD","NZ") %>% rownames(singlerDF)[.]
singlerDF[normal_cells,"cell.types"] %<>% gsub("MCL","B_cells",.)
# singler1sub false positive results  =========
table(singlerDF$singler1sub, object$orig.ident) %>% kable %>% kable_styling()
singlerDF[normal_cells,"singler1sub"] %<>% gsub("MCL:.*$","B_cells:Memory",.)
table(singlerDF$cell.types, object$orig.ident) %>% kable %>% kable_styling()
table(singlerDF$singler1sub, object$orig.ident)%>% kable %>% kable_styling()
table(singlerDF$cell.types %>% sort)%>% kable %>% kable_styling()
table(singlerDF$cell.types) %>% kable() %>% kable_styling()
##############################
# process color scheme
##############################
#singler_colors <- readxl::read_excel("doc/singler.colors.xlsx")
#singler_colors1 = as.vector(singler_colors$singler.color1[!is.na(singler_colors$singler.color1)])
#singler_colors1[duplicated(singler_colors1)]
#singler_colors2 = as.vector(singler_colors$singler.color2[!is.na(singler_colors$singler.color2)])
singler_colors2 = c("#E6AB02","#6A3D9A", "#2055da", "#FB9A99", "#A65628", "#B3B3B3", "#B3DE69", "#F0027F")
object <- AddMetaData(object = object,metadata = singlerDF["cell.types"])
object <- AddMetaColor(object = object, label= "cell.types", colors = singler_colors2)
Idents(object) <- "cell.types"

lapply(c(UMAPPlot.1,TSNEPlot.1), function(fun)
  fun(object = object, label = F, group.by = "cell.types",
           cols = ExtractMetaColor(object),no.legend = F,
           pt.size = 0.1,label.size = 3, do.print = T,do.return = F,
           title = "Cell type labeling by Blueprint + Encode + MCL"))

save(object,file="data/MCL_41_harmony_20200225.Rda")

##############################
# draw tsne plot
##############################
object <- subset(object,idents = c("HSC/progenitors","Nonhematopoietic cells"), invert = TRUE)
Idents(object) = "cell.types"
object %<>% sortIdent()
table(Idents(object))


cell_Freq <- table(Idents(object)) %>% as.data.frame
cell_Freq = cell_Freq[order(cell_Freq$Var1),]
cell_Freq$col = ExtractMetaColor(object)
cell_Freq = cell_Freq[order(cell_Freq$Freq,decreasing = T),]
cell_Freq$Var1 %<>% factor(levels = as.character(cell_Freq$Var1))
colnames(cell_Freq)[1:2] = c("Cell_Type", "Cell_Number")

jpeg(paste0(path,"cell_type_numbers.jpeg"), units="in", width=10, height=7,res=600)
ggbarplot(cell_Freq, "Cell_Type", "Cell_Number",
          fill = "Cell_Type", color = "Cell_Type",xlab = "",
          palette = cell_Freq$col,x.text.angle = 90,
          title = "Numbers of major cell types in total 43 samples")+NoLegend()+
  theme(plot.title = element_text(hjust = 0.5,size=15))
dev.off()
