library(Seurat)
library(DT)
library(ggplot2)
library(GGally)
library(kableExtra)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
# meta_data
TCGA_deconvolution_clinical <- read.csv("jupyter/TCGA_deconvolution_clinical.csv")

treatments_list = lapply(TCGA_deconvolution_clinical$treatments, function(x) {
                            l = eval(parse(text=x))
                            df <- data.frame(matrix(unlist(l), nrow=length(l)*2,  ncol = 1, 
                                                    byrow = FALSE))
                            df$treatments_id = paste0(rep(names(l),each=2),rep(c("_1","_2", time = length(l))))
                            rownames(df) = 
                            df$sample_id = sub("_treatment.*","",df["submitter_id","X1"])
                            df$treatments_id = names(l)
                            return(df)
                            }
                         )
treatments_df = bind_rows(treatments_list)

gather(treatments_df)

sapply(grep("day",colnames(TCGA_deconvolution_clinical),value =T), function(x){
    temp  = TCGA_deconvolution_clinical[,x]
    temp = temp[!is.na(temp)]
    head(as.numeric(temp))
})

head(TCGA_deconvolution_clinical[,c("sample",
                                    "days_to_last_follow_up",
                                    "days_to_death",
                                    "subtype_Combined.days.to.last.followup.or.death")],10) %>%
    kable %>% kable_styling()
seurat = readRDS("data/20201010_TCGA-BLCA&LAML.rds")
UMAPPlot.1(seurat, group.by = "orig.ident",do.print = T,label = T,label.repel = T,
           no.legend = T, title = "TCGA-BLCA and TCGA-LAML")

PrepareShiny(object,samples = samples, Rshiny_path = Rshiny_path,split.by = "orig.ident",
            reduction = "umap")
all_genes <- c("",rownames(object))
saveRDS(all_genes,file= paste0(Rshiny_path,"data/Rownames.rds"))


BLCA.prop_main = readRDS("data/20201118_BLCA_LAML.prop.rds")
MuSiC_cellIdent = apply(BLCA.prop_main$Est.prop.weighted,1, which.max)
NNLS_cellIdent = apply(BLCA.prop_main$Est.prop.allgene,1, which.max)
cellIdent = data.frame("MuSiC" = colnames(BLCA.prop_main$Est.prop.weighted)[MuSiC_cellIdent],
                       "NNLS"  = colnames(BLCA.prop_main$Est.prop.allgene)[NNLS_cellIdent],
                       row.names = rownames(BLCA.prop_main$Est.prop.allgene))
seurat@meta.data %<>% cbind(cellIdent)

UMAPPlot.1(seurat, group.by = "orig.ident",do.print = T,label = T,label.repel = T,
           width=7, height=7,
           no.legend = T, title = "TCGA-BLCA and TCGA-LAML")

UMAPPlot.1(seurat, group.by = "MuSiC",cols = Singler.colors,label = T, label.repel = T,
           width=7, height=7,
           no.legend = T, do.print = T,
           title = "Deconvolution results by MuSiC")

UMAPPlot.1(seurat, group.by = "NNLS",cols = Singler.colors,label = T, label.repel = T,
           width=7, height=7,
           no.legend = T, do.print = T,
           title = "Deconvolution results by NNLS")
# ============================
BLCA_seurat = readRDS("data/20200920_TCGA-BLCA.rds")
BLCA.prop_main = readRDS("data/20201118_TCGA-BLCA.prop_main.rds")
grep("day",colnames(BLCA_seurat@meta.data),value = T)
Remove = is.na(BLCA_seurat@meta.data[,"subtype_Combined.days.to.last.followup.or.death"]) |
    BLCA_seurat@meta.data[,"subtype_Combined.days.to.last.followup.or.death"] %in% c("[Not Available]","[Discrepancy]") |
    BLCA_seurat@meta.data[,"subtype_Combined.days.to.last.followup.or.death"] < 0
BLCA_seurat = BLCA_seurat[,!Remove]
prop_main = BLCA.prop_main$Est.prop.weighted[!Remove,]
cellGroups <- list("lymphocytes" =c("Monocytes","NK cells","T cells CD4+","Macrophages","CLP","MEP","B cells","T cells CD8+", "GMP","Myocytes","Plasma cells","HSC","DC","Megakaryocytes","Macrophages M1",
                                    "Eosinophils","CMP","Neutrophils","Erythrocytes","Macrophages M2"),
                   "Macrophages_total" = c("Macrophages","Macrophages M1","Macrophages M2"),
                   "HSC_Progenitor" = c("CLP","MEP","HSC","CMP","MPP","GMP")
) 
for(i in seq_along(cellGroups)){
    BLCA_seurat@meta.data[,names(cellGroups)[i]] = rowSums(prop_main[,cellGroups[[i]]])
    
}
meta_data = as.data.frame(cbind(BLCA_seurat@meta.data[,c(names(cellGroups),#"subtype_Neoantigen.load",
                                                        "subtype_Combined.days.to.last.followup.or.death"
                                                        # "days_to_death"
                                                         )],
                                prop_main))
meta_data = as.data.frame(sapply(meta_data,as.numeric))

meta_data[,"scale_log_days_to_death"] = scale(log(meta_data$subtype_Combined.days.to.last.followup.or.death+1))
meta_data = meta_data[,-grep("^subtype_Combined.days.to.last.followup.or.death$",colnames(meta_data))]
glm.fit <- glm(formula = scale_log_days_to_death ~ . -subtype_Combined.days.to.last.followup.or.death,
               data = meta_data, family = gaussian)
summary(glm.fit)

glm.fit1 <- glm(formula = scale_log_days_to_death ~ `Macrophages.M1` + Macrophages.M2 + Macrophages +
                    HSC ,
               data = meta_data, family = gaussian)
summary(glm.fit1)


jpeg(paste0(path,"glm_Macrophages.M1.jpeg"), units="in", width=10, height=7,res=600)
ggplot(meta_data,aes(x=subtype_Combined.days.to.last.followup.or.death,
                     y=`Macrophages.M1`,color=HSC))+geom_point()+
    stat_cor(method = "spearman",label.x = 3000, label.y = 0.025)+
    ggtitle("Days.until.death vs Macrophages M1 %") + TitleCenter()
dev.off()


cor_res <- Hmisc::rcorr(as.matrix(meta_data), type="spearman")
cor_res$P = -log10(cor_res$P)
jpeg(paste0(path,"ggcorr_Macrophages.M1.jpeg"), units="in", width=10, height=7,res=600)

GGally::ggcorr(data =meta_data,
               cor_matrix = cor_res$r, method = c("pairwise", "spearman"),
               digits = 1,
               #nbreaks = 6,
               hjust = 1,
               size = 2,
               angle = -45,
               label = TRUE,
               label_size = 1,
               color = "grey50")+
    ggtitle("Correlation of each attribute in deconvolution data")+
    theme(plot.title = element_text(hjust = 0.5))
dev.off()

GGally::ggcorr(data =meta_data,
               cor_matrix = cor_res$P, method = c("pairwise", "spearman"),
               digits = 1,
               #nbreaks = 6,
               hjust = 1,
               size = 2,
               angle = -45,
               label = TRUE,
               label_size = 1,
               color = "grey50")+
    ggtitle("Correlation of each attribute in deconvolution data")+
    theme(plot.title = element_text(hjust = 0.5))
