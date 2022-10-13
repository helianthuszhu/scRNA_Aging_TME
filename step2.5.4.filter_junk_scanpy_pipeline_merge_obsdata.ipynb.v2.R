setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')

figdir <- file.path("output2_step2.5.3.filter_junk_scanpy_pipeline/stat_tissue_distribution/v2/")
#
#
suppressPackageStartupMessages({
    library("reticulate")
    library("ggplot2")
    #library("scater")
    library("Seurat")
})
sc <- import("scanpy")
#
adataall= sc$read('output2_step2.5.3.filter_junk_scanpy_pipeline/data/adata.harmony.overclustered.filtered.CelltypeDefined.withrawcounts.h5ad')
#
stat= as.data.frame(adataall$obs)
head(stat)
stat$age_num=as.numeric(paste(stat$age))
stat$age_group=cut(
  stat$age_num,
  breaks = c(20, 49, 60,100),
  labels = c("Young", "intermediated", "Aged")
)

######
######
######ROE for mainlineage 
library(dplyr)
options(stringsAsFactors=FALSE)
library(reticulate)
library(tidyr)
library(ggpubr)

ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}
#
##############
summary <- table(stat[,c('age_group','Defined_Cell_Type')])
roe <- as.data.frame(ROIE(summary))
roe$Group <- rownames(roe)
roe_long=roe %>% pivot_longer(cols= colnames(roe[,!(colnames(roe) %in% c('Group'))]),
                                                                     names_to= "Cell_type",
                                                                values_to = "ROE")
roe_long$Group= factor(roe_long$Group, levels=c('Aged','intermediated','Young'))
head(roe_long)


library(Polychrome)
set.seed(723451)
fifth <- createPalette(15, c("#00ffff", "#ff00ff", "#ffff00"), M=1000)
age_color= c('#EE2EE8','#47EFFA','#F0E716')

p1=ggdotchart(roe_long, x = "Cell_type", y = "ROE",
           color = "Group",                                # Color by groups
           palette = age_color, # Custom color palette
           #palette =FLATUI_CLASS,
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 8,                                 # Large dot size
           group = "Group",
           label = round(roe_long$ROE,2),                        # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+ geom_hline(yintercept = 1, linetype = 2, color = "black")+#theme(legend.position='none')
ylab("Ro/e")+ggtitle("")+theme(axis.text.x = element_text(angle = 90,vjust = 1))

pdf(paste0(figdir, '/ROE.ageGroup_vs_Defined_Cell_Subtype_function2.pdf'),width=8,height = 7)
p1
dev.off()
#########
####################################
####combine cell subtype annotation
####
adata_Mal= sc$read('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Epithelial/malignant/data/adata_malignantcell_degs_calculated.h5ad')
#
adata_TNK=sc$read('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm.h5ad')

adata_Myeloid=sc$read('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Myeloidcell_v2/adata_Myeloid_v2_degs_signature_louvain_recluster_renamed_signature.4s_added_subtype_Defined_bad_cluster_rm.h5ad')

adata_Bcell=sc$read('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Bcell_v2/adata_Bcell_v2_degs_subtype_Defined_bad_cluster_rm.h5ad')

adata_Endo=sc$read('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Endothelial_v2/adata_Endothelial_v2_degs_subtype_Defined_aged_clustered.h5ad')

adata_Fib=sc$read('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Fibroblasts_v2/adata_Fibroblasts_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_agegrouped.h5ad')
#
#
meta_Mal=as.data.frame(adata_Mal$obs)
meta_Mal$Defined_Cell_Subtype_function= 'C1-Malignant'
meta_Mal=meta_Mal[, c('Index','age','Defined_Cell_Subtype_function')]

meta_TNK=as.data.frame(adata_TNK$obs)[, c('Index','age','Defined_Cell_Subtype_function')]
meta_Myeloid=as.data.frame(adata_Myeloid$obs)[, c('Index','age','Defined_Cell_Subtype_function')]
meta_Bcell=as.data.frame(adata_Bcell$obs)[, c('Index','age','Defined_Cell_Subtype_function')]
meta_Endo=as.data.frame(adata_Endo$obs)[, c('Index','age','Defined_Cell_Subtype_function')]
meta_Fib=as.data.frame(adata_Fib$obs)[, c('Index','age','Defined_Cell_Subtype_function')]
#
meta_Mast= as.data.frame(adataall$obs)
meta_Mast=meta_Mast[grep('Mast', meta_Mast$Defined_Cell_Type), c('Index','age')]
meta_Mast$Defined_Cell_Subtype_function= 'C1-Mast'
dim(meta_Mast)
table(stat$Defined_Cell_Type)
#
#
metadata= rbind(meta_Mal, meta_TNK, meta_Myeloid,meta_Bcell,  meta_Endo, meta_Fib, meta_Mast)
head(metadata)
metadata$mainlineage= sapply(stringr::str_split(metadata$Defined_Cell_Subtype_function, "-"), `[`, 2)

#
library(DataCombine)
Replaces <- data.frame(from = c('Malignant','CD4','CD8','CD4/CD8','NK','cDC',
                                'Monocyte','Mac','Cycling','cDC2','cDCs','cDC1','pDC','Bcell','ECs','CAF','myCAF','Pericyte','Mast'), 
                       to = c('Malignant','CD4','CD8','mixed','NK','DCs',
                                'Monocyte','Mac','mixed','DCs','DCs','DCs','DCs','Bcell','ECs','Fibs','Fibs','Fibs','Mast'))
metadata$mainlineage2<- FindReplace(data = metadata, Var = "mainlineage", replaceData = Replaces,
                                     from = "from", to = "to",vector=TRUE)

table(metadata$mainlineage2)
#table(metadata$Defined_Cell_Subtype_function)
#
metadata$age_num=as.numeric(paste(metadata$age))
metadata$age_group=cut(
  metadata$age_num,
  breaks = c(20, 49, 60,100),
  labels = c("Young", "intermediated", "Aged")
)
#
#
write.csv(metadata, paste0(figdir, '/metadata_mal_mast_othercells.csv'))
save(metadata, file= paste0(figdir, '/metadata_mal_mast_othercells.RData'))
###
###
###
load('output2_step2.5.filter_junk/TotalTissue.combined.not.doubletFiltered.RData')
sceall=TotalTissue.combined
#
sceall_sel= subset(sceall, subset = Index %in% metadata$Index)
colnames(sceall_sel@meta.data)
#
mm= metadata[rownames(sceall_sel@meta.data), c('Defined_Cell_Subtype_function',	
                                               'mainlineage',	'mainlineage2',	'age_num',	'age_group')]

sceall_sel <- AddMetaData(
  object = sceall_sel,
  metadata = mm
)


saveRDS(sceall_sel, file=paste0(figdir, '/sceall_sel.for.cellchat.rds') )


#
#
#
###########
###########
#
color_CLASS = c('#0000FF', # 0
                 '#FF0000', # 1   
                 '#FF8C00', # 2
                 '#FFD700', # 3  
                 '#808000', # 4
                 '#9ACD32', # 5   
                 '#44F544', # 6
                 '#20B2AA', # 7  
                 '#2F4F4F', # 8
                 '#00BFFF', # 9   
                 '#147ADF', # 10
                 '#191970', # 11 
                 '#6200A4', # 12
                 '#8B008B', # 13  
                 '#FF00FF', # 14
                 '#FF1493', # 15 
                 '#8B4513', # 16
                 '#D2691E', # 17 
                 '#B0C4DE', # 18
                 '#696969', # 19
                 '#FF586E', # 20
                 '#FAEBD7', # 21
                 '#E2C792', # 22
                 '#FFF7A4', # 23
                 '#83715A', # 24
                 '#BC8F8F', # 25
                 '#4C718E', # 26
                 '#519395', # 27
                 '#7FFFD4', # 28  
                 '#0F830F', # 29
                 '#D4FFF0', # 30
                 '#B72222', # 31
                 '#FDBABA', # 32
                 '#1B1010', # 33
                 '#CD57FF', # 34
                 '#F0FFFF')
###
library(ggplot2)
aaaa= read.csv('output2_step2.5.3.filter_junk_scanpy_pipeline/stat_tissue_distribution/stacked_plot_ordered/adata.obs.for.stacked.plot.csv', header=T, row.names=1)
head(aaaa)
table(aaaa$Sample,aaaa$age_group)

aaaa$age_group= factor(aaaa$age_group, levels= c('Aged','Intermediated', 'Young'))
aaaa$Sample= factor(aaaa$Sample, levels= rev(c('E1','Y1'
,'Y2'
,'Y3'
,'Y4'
,'E2'
,'E3'
,'E4'
,'LUNG_T06'
,'LUNG_T08'
,'LUNG_T19'
,'LUNG_T20'
,'LUNG_T30'
,'LUNG_T09'
,'LUNG_T18'
,'LUNG_T25'
,'LUNG_T28'
,'LUNG_T31'
,'LUNG_T34')))

pdf('output2_step2.5.3.filter_junk_scanpy_pipeline/stat_tissue_distribution/stacked_plot_ordered/adata.obs.for.stacked.plot_percentage.pdf',width=5,height = 4)

ggplot(aaaa, aes(x=Sample, fill=Defined_Cell_Type)) + 
  geom_bar(position = "fill")+scale_fill_manual(values = color_CLASS)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  coord_flip()+
  ylab('percentage')
dev.off()

pdf('output2_step2.5.3.filter_junk_scanpy_pipeline/stat_tissue_distribution/stacked_plot_ordered/adata.obs.for.stacked.plot_abs.pdf',width=5,height = 4)

ggplot(aaaa, aes(x=Sample, fill=Defined_Cell_Type)) + 
  geom_bar()+scale_fill_manual(values = color_CLASS)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  coord_flip()+
  ylab('percentage')
dev.off()
