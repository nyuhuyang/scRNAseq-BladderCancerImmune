library(Seurat)
library(DT)
library(ggplot2)
library(GGally)
library(kableExtra)
library(kable)

source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat3_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
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
meta_data1 = meta_data[,-grep("^subtype_Combined.days.to.last.followup.or.death$",colnames(meta_data))]
glm.fit <- glm(formula = scale_log_days_to_death ~ .,
               data = meta_data1, family = gaussian)
summary(glm.fit)

glm.fit1 <- glm(formula = scale_log_days_to_death ~ `Macrophages M1` + `Macrophages M2` + Macrophages +
                    HSC ,
               data = meta_data1, family = gaussian)
summary(glm.fit1)


jpeg(paste0(path,"glm_Macrophages.M1.jpeg"), units="in", width=10, height=7,res=600)
ggplot(meta_data,aes(x=subtype_Combined.days.to.last.followup.or.death,
                     y=`Macrophages M1`,color=HSC))+geom_point()+
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


# calculate survive curve=========


meta.data = as.data.frame(cbind(BLCA_seurat@meta.data[,c("days_to_death",
                                                         "subtype_Combined.days.to.last.followup.or.death",
                                                         names(cellGroups))],
                                prop_main))

calculate_survive_curve <- function(meta_data, Cell_types, day, save.path){
    library(survminer)
    library(survival)
    meta.data = meta_data[,c(Cell_types,"vital_status",day)]
    #Create a Survival Object
    meta.data$status = plyr::mapvalues(meta.data$vital_status,
                                             from = c("Alive","Dead","Not Reported"),
                                             to = c(0,1, NA))
    meta.data = meta.data[!is.na(meta.data$status),]
    meta.data$status %<>% as.integer()
    meta.data[,"Overall survival (year)"] = meta.data[,day]/365
    
    if(!dir.exists(save.path)) dir.create(save.path, recursive = T)
    for(i in seq_along(Cell_types)) {
        cell_type = Cell_types[i]
        cut_function = ifelse(median(meta.data[,cell_type]) ==0, "posi_neg","high_low")
        cut_label= switch (cut_function,
                           "posi_neg" = c("Positive","Negative"),
                           "high_low" = c("High","Low"))
        if(all(meta.data[,cell_type] ==0)) next
        meta.data[,cut_function] = switch (cut_function,
                                           "posi_neg" = plyr::mapvalues(meta.data[,cell_type] >0,
                                                                        from = c(TRUE,FALSE),
                                                                        to = c("Positive","Negative")),
                                           "high_low" = plyr::mapvalues(as.integer(cut(meta.data[,cell_type],2)),
                                                                        from = c(1,2),
                                                                        to = c("Low","High"))
        )
        fit <- switch (cut_function,
                       "posi_neg" = survfit(Surv(`Overall survival (year)`, status) ~ posi_neg, data = meta.data),
                       "high_low" = survfit(Surv(`Overall survival (year)`, status) ~ high_low, data = meta.data)
        )
        # if too little data for one side, skip
        lvl <- sort(as.vector(table(meta.data[,cut_function])))
        if(lvl[1] < 3 | lvl[2] < 3) next
        g <- ggsurvplot(
            fit,
            censor.shape="|", censor.size = 4,
            data = meta.data,
            size = 1,                 # change line size
            palette =
                c("#E7B800", "#2E9FDF"),# custom color palettes
            conf.int = TRUE,          # Add confidence interval
            pval = TRUE,              # Add p-value
            legend.title = cell_type,
            risk.table = TRUE,        # Add risk table
            risk.table.col = "strata",# Risk table color by groups
            legend.labs = cut_label,    # Change legend labels
            risk.table.height = 0.25, # Useful to change when you have multiple groups
            ggtheme = theme_bw(),      # Change ggplot2 theme
            xlab = "Overall survival (year)"
        )
        jpeg(paste0(save.path,"ggsurvplot_",cell_type,".jpeg"), units="in", width=5, height=7,res=600)
        print(g)
        dev.off()
        svMisc::progress(i, length(Cell_types))
    }
}


#meta_data[,"Survival rate"] = rank(meta_data$subtype_Combined.days.to.last.followup.or.death)
# range(meta_data[,"Survival rate"]) = [?, 411]
#meta_data[,"Survival rate"] = (meta_data[,"Survival rate"] - min(meta_data[,"Survival rate"]))/max(meta_data[,"Survival rate"])




