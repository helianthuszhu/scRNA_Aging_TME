############
######time
######2021-05-14
############method h5ad adata to seurat/singlecellexperiment object
#method 1
#
#run /data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/anndata2ri_h5ad2singlecellexperiment.ipynb
#method 2
#conda activate /home/chenh/.cache/basilisk/1.2.1/zellkonverter-1.0.3/anndata_env
#library(zellkonverter)
#mysce=readH5AD('adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd.h5ad')
#
######the two methods can do when next time with raw counts
######
######
###
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
###get raw counts object
###
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')
datadir='output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/'
figdir='output2_step2.5.3.filter_junk_scanpy_pipeline/figure_Tcell_v2/'

suppressPackageStartupMessages({
    library("reticulate")
    library("ggplot2")
    #library("scater")
    library("Seurat")
})

sc <- import("scanpy")


adata <- sc$read_h5ad(paste0(datadir, 'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd.h5ad'))
adata

###load raw counts seurat object
#
load('output2_step2.5.filter_junk/TotalTissue.combined.not.doubletFiltered.RData')
scetotal= TotalTissue.combined
#
#change raw counts and meta data
######
meta= adata$obs
head(meta)

all(rownames(meta) %in% colnames(scetotal))

sce=subset(scetotal[, rownames(meta)])
####
sce@meta.data= meta
####
# Add embedding
embedding <- adata$obsm["X_umap"]
rownames(embedding) <- adata$obs_names$to_list()
colnames(embedding) <- c("umap_1", "umap_2")
sce[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")
#
saveRDS(sce, file= paste0(datadir, 'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd_with_rawcounts.rds'))
#############
######signature score compare between three age groups
######
library(Seurat)
sce= readRDS(paste0(datadir, 'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd_with_rawcounts.rds'))
###
library(ggplot2)
library(ggpubr)
library(Polychrome)
set.seed(723451)
fifth <- createPalette(15, c("#00ffff", "#ff00ff", "#ffff00"), M=1000)
agecol= as.vector(fifth)
age_color= c('#F0E716','#47EFFA','#EE2EE8')
###
statdata= sce@meta.data
statdata$Tcellsubtype= sapply(stringr::str_split(statdata$Defined_Cell_Subtype_function, "[-]"), `[`, 2)
table(statdata$Tcellsubtype)

my_comparisons=list(c('Young', 'Aged'), c('Intermediated', 'Aged'), c('Young', 'Intermediated'))

pdf(paste0(figdir, '/boxplot_cytoscore.ageGroup_vs_Defined_Cell_Subtype_function.pdf'),width=25,height = 5)
ggplot(statdata)+
  aes(x = age_group, y = cytotoxic_score, fill = age_group) + # add color to boxes with fill
  #geom_boxplot(varwidth = FALSE) + # vary boxes width according to n obs.
  geom_boxplot(aes(fill=age_group),outlier.colour = "black",outlier.size = 1)+
  #geom_jitter(alpha = 0.8, width = 0,size=1) + # adds random noise and limit its width
  facet_wrap(~Defined_Cell_Subtype_function, nrow=1) + # divide into 2 panels
  theme(legend.position = "none") +# remove legend
  scale_fill_manual(values =  age_color)+
  stat_compare_means(aes(group = age_group),label.y = 5, label = "p.signif")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  #ggtitle(paste0('Age','\n','No.',dim(clindata)[1]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())

dev.off()
#
pdf(paste0(figdir, '/boxplot_cytoscore.ageGroup_vs_Defined_Cell_Subtype_function_NKcells.pdf'),width=4,height = 5)
ggplot(statdata[grep('NK', statdata$Defined_Cell_Subtype_function),])+
  aes(x = age_group, y = cytotoxic_score, fill = age_group) + # add color to boxes with fill
  #geom_boxplot(varwidth = FALSE) + # vary boxes width according to n obs.
  geom_boxplot(aes(fill=age_group),outlier.colour = "black",outlier.size = 1)+
  #geom_jitter(alpha = 0.8, width = 0,size=1) + # adds random noise and limit its width
  facet_wrap(~Defined_Cell_Subtype_function, nrow=1) + # divide into 2 panels
  theme(legend.position = "none") +# remove legend
  scale_fill_manual(values =  age_color)+
  stat_compare_means(aes(group = age_group),label.y = 5, label = "p.signif")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  #ggtitle(paste0('Age','\n','No.',dim(clindata)[1]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())

dev.off()




pdf(paste0(figdir, '/boxplot_cyto_exh_score.ageGroup_vs_Tcellsubtype.pdf'),width=6,height = 5)
ggplot(statdata)+
  aes(x = age_group, y = cytotoxic_score, fill = age_group) + # add color to boxes with fill
  #geom_boxplot(varwidth = FALSE) + # vary boxes width according to n obs.
  geom_boxplot(aes(fill=age_group),outlier.colour = "black",outlier.size = 0.5)+
  #geom_jitter(alpha = 0.8, width = 0,size=1) + # adds random noise and limit its width
  facet_wrap(~Tcellsubtype, nrow=1) + # divide into 2 panels
  theme(legend.position = "none") +# remove legend
  scale_fill_manual(values =  age_color)+
  stat_compare_means(aes(group = age_group),label.y = 5, label = "p.signif")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  #ggtitle(paste0('Age','\n','No.',dim(clindata)[1]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())

ggplot(statdata)+
  aes(x = age_group, y = exhausted_score, fill = age_group) + # add color to boxes with fill
  #geom_boxplot(varwidth = FALSE) + # vary boxes width according to n obs.
  geom_boxplot(aes(fill=age_group),outlier.colour = "black",outlier.size = 0.5)+
  #geom_jitter(alpha = 0.8, width = 0,size=1) + # adds random noise and limit its width
  facet_wrap(~Tcellsubtype, nrow=1) + # divide into 2 panels
  theme(legend.position = "none") +# remove legend
  scale_fill_manual(values =  age_color)+
  stat_compare_means(aes(group = age_group),label.y = 5, label = "p.signif")+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  #ggtitle(paste0('Age','\n','No.',dim(clindata)[1]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())


ggplot(statdata)+
  aes(x = age, y = cytotoxic_score, fill = age) + # add color to boxes with fill
  #geom_boxplot(varwidth = FALSE) + # vary boxes width according to n obs.
  geom_boxplot(aes(fill=age_group),outlier.colour = "black",outlier.size = 0.5)+
  #geom_jitter(alpha = 0.8, width = 0,size=1) + # adds random noise and limit its width
  facet_wrap(~Tcellsubtype, nrow=1) + # divide into 2 panels
  theme(legend.position = "none") +# remove legend
  scale_fill_manual(values =  age_color)+
  stat_compare_means(aes(group = age),label.y = 5, label = "p.signif")+
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  #ggtitle(paste0('Age','\n','No.',dim(clindata)[1]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())

ggplot(statdata)+
  aes(x = age, y = exhausted_score, fill = age) + # add color to boxes with fill
  #geom_boxplot(varwidth = FALSE) + # vary boxes width according to n obs.
  geom_boxplot(aes(fill=age_group),outlier.colour = "black",outlier.size = 0.5)+
  #geom_jitter(alpha = 0.8, width = 0,size=1) + # adds random noise and limit its width
  facet_wrap(~Tcellsubtype, nrow=1) + # divide into 2 panels
  theme(legend.position = "none") +# remove legend
  scale_fill_manual(values =  age_color)+
  stat_compare_means(aes(group = age),label.y = 5, label = "p.signif")+
  #stat_compare_means(comparisons = my_comparisons,label = "p.signif")+
  #ggtitle(paste0('Age','\n','No.',dim(clindata)[1]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())


dev.off()

#################################################################################################################
#################################################################################################################
#################################################################################################################
######ROE 
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
stat= sce@meta.data
summary <- table(stat[,c('age_group','Defined_Cell_Subtype_function')])
roe <- as.data.frame(ROIE(summary))
roe$Group <- rownames(roe)
roe_long=roe %>% pivot_longer(cols= colnames(roe[,!(colnames(roe) %in% c('Group'))]),
                                                                     names_to= "Cell_type",
                                                                values_to = "ROE")
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

######only show in Aged group
td= subset(roe_long, Group=='Aged')

p1=ggdotchart(td, x = "Cell_type", y = "ROE",
           color = "Group",                                # Color by groups
           palette = age_color, # Custom color palette
           #palette =FLATUI_CLASS,
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 8,                                 # Large dot size
           group = "Group",
           label = round(td$ROE,2),                        # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+ geom_hline(yintercept = 1, linetype = 2, color = "black")+#theme(legend.position='none')
ylab("Ro/e (in Aged group)")+ggtitle("")+theme(axis.text.x = element_text(angle = 90,vjust = 1))

pdf(paste0(figdir, '/ROE.ageGroup_Aged_vs_Defined_Cell_Subtype_function2.pdf'),width=8,height = 7)
p1
dev.off()
#################################################################################################################
#################################################################################################################
#################################################################################################################
#########calculate the percentage of each cell subtype if each sample
#########
stat1=sce@meta.data
#stat1$pos=paste(meta.sub$defined.subtype_detailed,meta.sub$subtype_Clusters,sep="_")
#stat1$pos=paste0('cluster',stat1$louvain)
#stat1$gg=stat1$Sample
#!!!!!stacked malignancy across samples
x=as.data.frame.matrix(table(stat1$Defined_Cell_Subtype_function,stat1$Sample))
#colnames(x)=paste0("pt", colnames(x))
# Transform this data in %
#data_percentage <- apply(x, 2, function(x){x*100/sum(x,na.rm=T)})
data_percentage <- apply(x, 2, function(x){x/sum(x,na.rm=T)})
#data_percentage = as.data.frame(data_percentage )
data_percentage=as.data.frame(t(data_percentage))
#
ageinfo= stat1[!(duplicated(stat1$Sample)), c('Sample', 'intage')]
rownames(ageinfo)= ageinfo$Sample
#
prostat= data.frame(age= ageinfo[rownames(data_percentage), 'intage'],data_percentage )
#
prostat

#save(prostat, file=paste0(datadir, 'prostat.RData'))
#
library(tidyr)
prostat_long <- prostat %>% pivot_longer(cols=colnames(prostat[-1]),
                            names_to= "Defined_Cell_Subtype_function",
                            values_to = "priportion")
prostat_long$Defined_Cell_Subtype_function=gsub('[.]','-',prostat_long$Defined_Cell_Subtype_function)
head(prostat_long)
#
prostat_long$Defined_Cell_Subtype_function= factor(prostat_long$Defined_Cell_Subtype_function, 
     levels=c('C1-CD4-CCR7','C2-CD4-HSPA1A', 'C3-CD4-CXCR6','C4-CD4-CXCL13','C5-CD4-FOXP3', 'C6-CD4-RORC',
     'C7-CD8-SLC4A10','C8-CD8-GZMK', 'C9-CD8-CX3CR1','C10-CD8-ZNF683','C11-CD8-LAYN', 'C12-CD8-TRDV1',
     'C13-NK-FCGR3A','C14-NK-XCL1','C15-CD4-CD8-ISG15','C16-CD4-CD8-MKI67'))
###split based on cell type
library(ggpubr)
pdf(paste0(figdir, '/Cor_age_vs_cellPercentage2.pdf'),width=8,height = 8)
ggscatter(prostat_long, x = "age", y = "priportion", size = 0.5,
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "spearman",
          xlab = "age", ylab = "percentage (%) per sample")+
  ggtitle("age_vs_cellPercentage_spearman")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")+facet_wrap("Defined_Cell_Subtype_function")

ggscatter(prostat_long, x = "age", y = "priportion", size = 0.5,
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "pearson",
          xlab = "age", ylab = "percentage (%) per sample")+
  ggtitle("age_vs_cellPercentage_pearson")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")+facet_wrap("Defined_Cell_Subtype_function")
dev.off()


tindex= unique(prostat_long$Defined_Cell_Subtype_function)
plist= list()
corlist= list()
for (i in 1: length(tindex)){
   pp= ggscatter(subset(prostat_long,  Defined_Cell_Subtype_function %in%  tindex[i]), x = "age", y = "priportion", size = 0.5,
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "pearson",
          xlab = "age", ylab = "percentage (%) per sample")+
  ggtitle(paste0("Pearson cor",'\n', 'age vs ',tindex[i]))+
    #geom_point(aes(color="#dd1c77"))
    theme(plot.title = element_text(size = 8, face = "bold"))+
  geom_point(fill="#dd1c77",color="#dd1c77")
    plist[[i]]=pp
    ##
    
  ttmp= subset(prostat_long,  Defined_Cell_Subtype_function %in%  tindex[i])
  rescor= cor.test(ttmp$age,ttmp$priportion, method="pearson")
  rescor_dataframe = data.frame(cor= rescor$estimate,  Pvalue= rescor$p.value)
  rownames(rescor_dataframe)=tindex[i]
  corlist[[i]]=rescor_dataframe  
}

library(patchwork)
pdf(paste0(figdir, '/Cor_age_vs_cellPercentage2.indi.pdf'),width=10,height = 10)

plist[[1]]+plist[[2]]+plist[[3]]+plist[[4]]+plist[[5]]+plist[[6]]+plist[[7]]+plist[[8]]+plist[[9]]+plist[[10]]+plist[[11]]+plist[[12]]+
plist[[13]]+plist[[14]]+plist[[15]]+plist[[16]]+
  plot_layout(ncol = 4)

dev.off()

#
cordataframe= do.call(rbind, corlist)
cordataframe=as.data.frame(cordataframe)
cordataframe$sig= ifelse(cordataframe$Pvalue < 0.05, 'significant', 'non significant')
cordataframe$name= rownames(cordataframe)
write.csv(cordataframe, paste0(figdir, '/Cor_age_vs_cellPercentage2.indi.csv'))
#
#https://r-coder.com/dot-plot-r/
#https://r-graphics.org/recipe-bar-graph-dot-plot

tcellcor= read.csv(paste0(figdir, '/Cor_age_vs_cellPercentage2.indi.csv'), header = T)
head(tcellcor)
tcellcor$cor2=abs(tcellcor$cor)
tcellcor$name= factor(tcellcor$name, levels = tcellcor[order(tcellcor$cor, decreasing = T),'name'])
#######
minvalue= round(min(tcellcor$cor), digits = 1)
maxvalue= round(max(tcellcor$cor), digits = 1)
library(ggplot2)
library(ggpubr)
pdf(paste0(figdir, '/Cor_needle_age_vs_cellPercentage.pdf'),width = 6,height = 5)
ggplot(tcellcor, aes(x = cor, y = name,color=Pvalue)) +
    geom_segment(aes(yend = name), xend = 0, colour = "grey50") +
    #geom_point(size = Pvalue, aes(colour = Pvalue)) +
    geom_point(aes( size = cor2), alpha = 1) +
    #scale_colour_brewer(palette = "Set1", limits = c("significant", "non significant")) +
    theme_classic2()+
    scale_colour_gradientn(limits=c(5-05,0.05), colours = c("darkred", "orange", "yellow", "white"),
                           breaks = c( 0.01, 0.03,0.05),
                           labels = c(0.01, 0.03,0.05))+
    #scale_colour_gradient2(limits=c(5-05,0.05), low="darkred",  mid ="orange", high ="Black")
    xlab('Pearson coef (cell proportion with age)')+ylab('Cluster')+
    labs(color="P value", size= 'Coef')+
    geom_vline(xintercept=0, linetype="dashed", color = "black")+
    scale_x_continuous(breaks=seq(minvalue, maxvalue, 0.2), limits = c(minvalue, maxvalue))
dev.off()
#################################################################################################################
## DirichletReg
stat1=sce@meta.data
#stat1$pos=paste(meta.sub$defined.subtype_detailed,meta.sub$subtype_Clusters,sep="_")
#stat1$pos=paste0('cluster',stat1$louvain)
#stat1$gg=stat1$Sample
#!!!!!stacked malignancy across samples
x=as.data.frame.matrix(table(stat1$Defined_Cell_Subtype_function,stat1$Sample))
#colnames(x)=paste0("pt", colnames(x))
# Transform this data in %
#data_percentage <- apply(x, 2, function(x){x*100/sum(x,na.rm=T)})
data_percentage <- apply(x, 2, function(x){x/sum(x,na.rm=T)})
#data_percentage = as.data.frame(data_percentage )
data_percentage=as.data.frame(t(data_percentage))
#
ageinfo= stat1[!(duplicated(stat1$Sample)), c('Sample', 'intage')]
rownames(ageinfo)= ageinfo$Sample
head(ageinfo)
#
prostat= data.frame(Age= as.numeric(ageinfo[rownames(data_percentage), 'intage']),data_percentage )
#
prostat


ttd= prostat
######
library("DirichletReg")
AL <- DR_data(ttd[, -c(1)])
#plot(AL, cex = 0.5, a2d = list(colored = FALSE, c.grid = FALSE))
lake1 <- DirichReg(AL ~ Age, ttd)
lake1
coef(lake1)

coda= data.frame(intercept=lake1$coefficients)
head(coda)

lake2 <- update(lake1, . ~ . + I(Age^2) | . + I(Age^2) | . + I(Age^2)|. + I(Age^2)|. + I(Age^2) | . + I(Age^2) | . + I(Age^2)|. + I(Age^2)|. + I(Age^2) | . + I(Age^2) | . + I(Age^2)|. + I(Age^2)|. + I(Age^2) | . + I(Age^2) | . + I(Age^2)|. + I(Age^2)
)
anova(lake1, lake2)
summary(lake2)
#
mycol= c('#0000FF', # 0
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
                 '#F0FFFF' # 35 
)

figdir= paste0(figdir, '/DirichletReg/')
if(!dir.exists(figdir)) dir.create(figdir,recursive=TRUE)

write.csv(coda, paste0(figdir, '/DirichletReg_TNK.coef.csv'))

pdf(paste0(figdir, '/DirichletReg_TNK.pdf'),width=7,height = 9)
par(mar = c(4, 4, 4, 4) + 0.1)
plot(rep(ttd$Age, dim(AL)[2]), as.numeric(AL), pch = 21, bg = rep(mycol[1:dim(AL)[2]], each = dim(AL)[1]), xlab = "age (yr)", ylab = "Proportion",
       ylim = c(0, 1), main = "Tcell Composition in Leader et al Lung Cancer by DirichletReg")

Xnew <- data.frame(Age = seq(min(ttd$Age), max(ttd$Age),
                  length.out = 100))
 for (i in 1:dim(AL)[2]) lines(cbind(Xnew, predict(lake1, Xnew)[, i]), col = mycol[1:dim(AL)[2]][i], lwd = 2)
legend('topright', legend = colnames(AL), lwd = 2, 
       col =  mycol[1:dim(AL)[2]], pt.bg =  mycol[1:dim(AL)[2]], pch = 21,bty = "n")
par(new = TRUE)
plot(cbind(Xnew, predict(lake1, Xnew, F, F, T)), lty = "24", type = "l", ylim = c(0,
                max(predict(lake1, Xnew, F, F, T))), axes = F, ann = F, lwd = 2)
axis(4)
mtext(expression(paste("Precision (", phi, ")", sep = "")), 4, line = 3)
legend("top", legend = c(expression(hat(mu[c] == hat(alpha)[c]/hat(alpha)[0])),
                     expression(hat(phi) == hat(alpha)[0])), lty = c(1, 2), lwd = c(3, 2), bty = "n")
dev.off()









#################################################################################################################
#################################################################################################################
####cumulative curve
#https://stackoverflow.com/questions/17033513/how-to-specify-color-of-lines-and-points-in-ecdf-ggplot2
####
library(Seurat)
sce = readRDS(paste0(datadir, 'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd_with_rawcounts.rds'))
####
library(ggplot2)
library(tidyr)
tm=sce@meta.data[,c('age_group','Defined_Cell_Subtype_function','naive_score','cytotoxic_score','exhausted_score','NK_score')]
head(tm)
####################1 combined plot
prostat_long <- as.data.frame(tm) %>% pivot_longer(cols=colnames(tm[-c(1,2)]),
                            names_to= "fun_state",
                            values_to = "score")
head(prostat_long)
#
pdf(paste0(figdir, '/CumulativeFractions.pdf'),width=20,height = 6)

ggplot(prostat_long, aes(score, colour = Defined_Cell_Subtype_function)) + 
             stat_ecdf(geom = "step")+
             labs(title="Empirical Cumulative Density Function",
     y = "Cumulative fraction", x="Score")+
            theme_classic()+
     scale_color_manual(values =color_CLASS)+
     #scale_size_manual(values = rep(5,length(unique(sce@meta.data$Defined_Cell_Subtype_function))))
     geom_line( lwd=1.5, stat="ecdf")+
     #facet_wrap("fun_state")
     facet_grid(.~ fun_state)

ggplot(prostat_long, aes(score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title="Empirical Cumulative Density Function",
     y = "Cumulative fraction", x="Score")+
            theme_classic()+
     scale_color_manual(values =color_CLASS)+
     #scale_size_manual(values = rep(5,length(unique(sce@meta.data$Defined_Cell_Subtype_function))))
     geom_line( lwd=1.5, stat="ecdf")+
     #facet_wrap("fun_state")
     facet_grid(.~ fun_state)

dev.off()
####################2
age_color= c('#F0E716','#EE2EE8','#47EFFA')
#
p1=ggplot(tm, aes(naive_score, colour = Defined_Cell_Subtype_function)) + 
             stat_ecdf(geom = "step")+
             labs(title="Empirical Cumulative Density Function",
     y = "Cumulative fraction", x="Naive_score")+
            theme_classic()+
     scale_color_manual(values =color_CLASS)+
     geom_line( lwd=1.5, stat="ecdf")
#
p2=ggplot(tm, aes(cytotoxic_score, colour = Defined_Cell_Subtype_function)) + 
             stat_ecdf(geom = "step")+
             labs(title="Empirical Cumulative Density Function",
     y = "Cumulative fraction", x="cytotoxic_score")+
            theme_classic()+
     scale_color_manual(values =color_CLASS)+
     geom_line( lwd=1.5, stat="ecdf")
#
p3=ggplot(tm, aes(exhausted_score, colour = Defined_Cell_Subtype_function)) + 
             stat_ecdf(geom = "step")+
             labs(title="Empirical Cumulative Density Function",
     y = "Cumulative fraction", x="exhausted_score")+
            theme_classic()+
     scale_color_manual(values =color_CLASS)+
     geom_line( lwd=1.5, stat="ecdf")
#
p4=ggplot(tm, aes(NK_score, colour = Defined_Cell_Subtype_function)) + 
             stat_ecdf(geom = "step")+
             labs(title="Empirical Cumulative Density Function",
     y = "Cumulative fraction", x="NK_score")+
            theme_classic()+
     scale_color_manual(values =color_CLASS)+
     geom_line( lwd=1.5, stat="ecdf")
#

p5=ggplot(tm, aes(naive_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title="Empirical Cumulative Density Function",
     y = "Cumulative fraction", x="Naive_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")
#
p6=ggplot(tm, aes(cytotoxic_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title="Empirical Cumulative Density Function",
     y = "Cumulative fraction", x="cytotoxic_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")
#
p7=ggplot(tm, aes(exhausted_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title="Empirical Cumulative Density Function",
     y = "Cumulative fraction", x="exhausted_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")
#
p8=ggplot(tm, aes(NK_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title="Empirical Cumulative Density Function",
     y = "Cumulative fraction", x="NK_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")
#
pdf(paste0(figdir, '/CumulativeFractions2.pdf'),width=28,height = 5)

gridExtra::grid.arrange(p1,p2,p3,p4,
             nrow = 1)   
dev.off()

pdf(paste0(figdir, '/CumulativeFractions_age.pdf'),width=24,height = 4)

gridExtra::grid.arrange(p5,p6,p7,p8,
             nrow = 1)   
dev.off()
#########
#2021-05-15
#
#subset CD8 T cells
tm_CD8= tm[grepl('-CD8-', tm$Defined_Cell_Subtype_function),]
head(tm_CD8)
#
gn= 'CD8_Tcells'
#
p5=ggplot(tm_CD8, aes(naive_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title=paste0(gn),
     y = "Cumulative fraction", x="Naive_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")

p6=ggplot(tm_CD8, aes(cytotoxic_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title=paste0(gn),
     y = "Cumulative fraction", x="cytotoxic_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")
#
p7=ggplot(tm_CD8, aes(exhausted_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title=paste0(gn),
     y = "Cumulative fraction", x="exhausted_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")
#
p8=ggplot(tm_CD8, aes(NK_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title=paste0(gn),
     y = "Cumulative fraction", x="NK_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")
#
pdf(paste0(figdir, '/CumulativeFractions_age_CD8subset.pdf'),width=24,height = 4)

gridExtra::grid.arrange(p5,p6,p7,p8,
             nrow = 1)   
dev.off()

#
#subset CD4 T cells
tm_CD4= tm[grepl('-CD4-', tm$Defined_Cell_Subtype_function),]
head(tm_CD4)
#
gn= 'CD4_Tcells'
#
p5=ggplot(tm_CD4, aes(naive_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title=paste0(gn),
     y = "Cumulative fraction", x="Naive_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")

p6=ggplot(tm_CD4, aes(cytotoxic_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title=paste0(gn),
     y = "Cumulative fraction", x="cytotoxic_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")
#
p7=ggplot(tm_CD4, aes(exhausted_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title=paste0(gn),
     y = "Cumulative fraction", x="exhausted_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")
#
p8=ggplot(tm_CD4, aes(NK_score, colour = age_group)) + 
             stat_ecdf(geom = "step")+
             labs(title=paste0(gn),
     y = "Cumulative fraction", x="NK_score")+
            theme_classic()+
     scale_color_manual(values =age_color)+
     geom_line( lwd=1.5, stat="ecdf")
#
pdf(paste0(figdir, '/CumulativeFractions_age_CD4subset.pdf'),width=24,height = 4)

gridExtra::grid.arrange(p5,p6,p7,p8,
             nrow = 1)   
dev.off()

#################################################################################################################
#################################################################################################################
#################################################################################################################
#####density plot test naive CD8 T cells
#
library("Nebulosa")
#
pdf(paste0(figdir, '/_density_Tcells_naive.pdf'),height=10,width=21)
p5 <- plot_density(sce,reduction = "umap", c('CCR7','TCF7','SELL','LEF1'), joint = TRUE)
p5 + plot_layout(ncol = 3)
dev.off()

#################################################################################################################
#################################################################################################################
#################################################################################################################
######
#########
#2021-05-16
#########
#########aged degs between three group for each cell type
#########
datadir='output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/aged_degs_seurat/'
if(!dir.exists(datadir)) dir.create(datadir,recursive=TRUE)
figdir='output2_step2.5.3.filter_junk_scanpy_pipeline/figure_Tcell_v2/aged_degs_seurat/'
if(!dir.exists(figdir)) dir.create(figdir,recursive=TRUE)


library(Seurat)
sce = readRDS(paste0('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/', 'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd_with_rawcounts.rds'))
sce
#
sce <- AddMetaData(
  object = sce,
  metadata = gsub('/', '-', sce$Defined_Cell_Subtype_function),
  col.name = 'Defined_Cell_Subtype_function'
)
Idents(sce)='age_group'
levels(sce)
####
scelist= SplitObject(sce, split.by = "Defined_Cell_Subtype_function")


cmarker_list <- lapply(X = scelist, FUN = function(x) {
    x <- FindAllMarkers(x)
    ####filter based \logFC\ > 0.25 & adj pvalue < 0.05
    x <- subset(x, abs(avg_log2FC) > 0.25 & p_val_adj < 0.05 )
    
})

cName= names(cmarker_list)
for (i in 1:length(cName)) {    
    aa= as.data.frame(cmarker_list[cName[i]])
    
    write.csv(aa,  paste0(datadir, cName[i],'_significant_degs_aged_cellSubtype.csv'))
}
save(cmarker_list, file=paste0(datadir, 'significant_degs_aged_cellSubtype_seurat.RData'))
##merge degs in to a dataframe

#########################################
#########################################
#up genes

cmarker_list_test <- lapply(X = scelist, FUN = function(x) {
    x <- FindAllMarkers(x)
    ####filter based \logFC\ > 0.25 & adj pvalue < 0.05
    x <- subset(x, avg_log2FC > 0.25 & p_val_adj < 0.05 )
    
    x <- data.frame(x, group= 'up')
    
    x <- x[, 6:8]
    
    #library(tidyr)
    x <- spread(x, key = cluster, value = group)
    #x= as.data.frame(x)
})

#####
cmlist= list()
for(i in 1:length(cmarker_list_test)) {
    tmd= cmarker_list_test[i]
    gn= names(cmarker_list_test[i])
    aa= as.data.frame(tmd, check.names=F)
    colnames(aa)[1]='gname'
    cmlist[[i]]= aa
    
}

df_merge  <- Reduce(function(...) merge(..., by="gname", all=TRUE), 
                    cmlist)

rownames(df_merge)=df_merge$gname
df_merge=df_merge[,-1]
head(df_merge)

all_groupName= unique(paste(sce$Defined_Cell_Subtype_function, sce$age_group, sep='.'))
length(all_groupName)

missedName= all_groupName[!all_groupName %in% colnames(df_merge)]

df_miss= data.frame(matrix(NA, nrow = nrow(df_merge), ncol =  length(missedName)),row.names= rownames(df_merge))
colnames(df_miss)= missedName
head(df_miss)

df_merge= cbind(df_merge, df_miss)
head(df_merge)
write.csv(df_merge, file= paste0(datadir, 'significant_up_degs_aged_cellSubtype_seurat.merged.csv'))
########################
df_merge=read.csv(paste0(datadir, 'significant_up_degs_aged_cellSubtype_seurat.merged.csv'), header=T, row.names=1, check.names=F)
####summarize the overlaps
#cn= names(cmarker_list_test)
cn= unique(sce$Defined_Cell_Subtype_function)
cnlist=list()
for(i in 1: length(cn)){
    tmd=df_merge[, grep(cn[i], colnames(df_merge))]
    head(tmd)
    nc1= sjmisc::row_count(tmd, count=NA, var=cn[i], append=FALSE)
    nc1[,paste(cn[i], 'state', sep='_')]= ifelse(nc1[,1] == ncol(tmd), 'not_sig', 'sig')
    #head(nc1)
    cnlist[[i]]= nc1
}
ctype_count= do.call(cbind, cnlist)
ctype_count= ctype_count[, grepl('state', colnames(ctype_count))]
head(ctype_count)
ctype_count_group= sjmisc::row_count(ctype_count, count='sig', var='ctype_count_group', append=FALSE)
#
up_count_group = sjmisc::row_count(df_merge, count='up', var='up_count', append=FALSE)
#down_count_group = sjmisc::row_count(df_merge, count='down', var='down_count', append=FALSE)
#
#stat_count= data.frame(ctype_count_group=ctype_count_group,up_count_group=up_count_group, down_count_group= down_count_group )
stat_count= data.frame(ctype_count_group=ctype_count_group,up_count_group=up_count_group)



stat_count$group='type3'
stat_count$group=ifelse(stat_count$ctype_count_group >=2,'type1', ifelse(stat_count$ctype_count_group ==1 , 'type2', 'type3'))

write.csv(stat_count, file= paste0(datadir, 'significant_up_degs_aged_cellSubtype_seurat.merged_stat_count.csv'))

head(stat_count)
dim(stat_count)
######
######
########
library(RColorBrewer)
library(pheatmap)
library(viridis)
#hmcols=colorRampPalette(rev(brewer.pal(10, "RdBu")))(256)
hmcols=c(colorRampPalette(c("white", rev(plasma(323, begin = 0.15))[1]))(10), rev(plasma(323, begin = 0.18)))
hmcols<-colorRampPalette(c("white","white","red"))(100)
#######
annotation_row= data.frame(tmp= colnames(df_merge))
rownames(annotation_row)= annotation_row$tmp
annotation_row$cell_type= sapply(stringr::str_split(annotation_row$tmp, "[.]"), `[`, 1)
annotation_row$age= sapply(stringr::str_split(annotation_row$tmp, "[.]"), `[`, 2)
annotation_row$cluster= sapply(stringr::str_split(annotation_row$tmp, "[-]"), `[`, 1) 
annotation_row$cluster= sapply(stringr::str_split(annotation_row$cluster, "[C]"), `[`, 2)
#annotation_row$cluster= gsub('-','',  annotation_row$cluster) 
annotation_row$cluster=as.numeric(paste(annotation_row$cluster))                        
annotation_row=annotation_row[order(annotation_row$cluster, decreasing=F),c('cell_type','age')] 
                           
head(annotation_row)


ann_colors  = list(age =c('Young'='#47EFFA','Intermediated'='#EE2EE8','Aged'='#F0E716'),
cell_type= c('C1-CD4-CCR7'=color_CLASS[1],
             'C2-CD4-HSPA1A'=color_CLASS[2],
             'C3-CD4-CXCR6'=color_CLASS[3],
             'C4-CD4-CXCL13'=color_CLASS[4],
             'C5-CD4-FOXP3'=color_CLASS[5],
             'C6-CD4-RORC'=color_CLASS[6],
             'C7-CD8-SLC4A10'=color_CLASS[7],
             'C8-CD8-GZMK'=color_CLASS[8],
             'C9-CD8-CX3CR1'=color_CLASS[9],
             'C10-CD8-ZNF683'=color_CLASS[10],
             'C11-CD8-LAYN'=color_CLASS[11],
             'C12-CD8-TRDV1'=color_CLASS[12],
             'C13-NK-FCGR3A'=color_CLASS[13],
             'C14-NK-XCL1'=color_CLASS[14],
             'C15-CD4-CD8-ISG15'=color_CLASS[15],
             'C16-CD4-CD8-MKI67'=color_CLASS[16])
             )
#
tt= stat_count[order(stat_count$ctype_count_group, decreasing=T),]

tt_sel= subset(tt, group=='type1')

heatmap_matrix=df_merge[rownames(tt_sel),rownames(annotation_row)]
heatmap_matrix[is.na(heatmap_matrix)] <- 0                    
heatmap_matrix[heatmap_matrix=="up"]<-1

#heatmap_matrix[] <- lapply(heatmap_matrix, function(x) as.numeric(as.character(x)))

#annotation_col= tt_sel

ph1 <- pheatmap::pheatmap(as.matrix(t(heatmap_matrix)), 
	           main= paste0('No.of.degs:', dim(heatmap_matrix)[1], '(shared by at least two cell type)'),
               #useRaster = T,
                treeheight_col=0,
               cluster_cols=T, 
               cluster_rows=F, 
               show_rownames=F, 
               show_colnames=F, 
               fontsize=5,
                          legend=F,
                          annotation_legend=F,
               #clustering_distance_rows=row_dist,
               #clustering_method = 'ward.D2',
               #cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               #breaks=bks,
               border_color = NA,
               color=hmcols,
               cellwidth=0.4,
               cellheight=5,
               annotation_colors = ann_colors,
               annotation_row= annotation_row
               #annotation_col=annotation_col
                        )
pdf(paste0(figdir, '/degs_seurat_summary11.pdf'),width=10)
ph1
dev.off()
                                               
#
                       
                           
tt_sel= subset(tt, group=='type2')

heatmap_matrix=df_merge[rownames(tt_sel),rownames(annotation_row)]
heatmap_matrix[is.na(heatmap_matrix)] <- 0                    
heatmap_matrix[heatmap_matrix=="up"]<-1#

#heatmap_matrix[] <- lapply(heatmap_matrix, function(x) as.numeric(as.character(x)))
                           
                         
ph2 <- pheatmap::pheatmap(as.matrix(t(heatmap_matrix)), 
	           main= paste0('No.of.degs:', dim(heatmap_matrix)[1],'(unique in one group in each cell type)'),
               #useRaster = T,
               treeheight_col=0,
               cluster_cols=T, 
               cluster_rows=F, 
               show_rownames=T, 
               show_colnames=F, 
               fontsize=5,
               #clustering_distance_rows=row_dist,
               #clustering_method = 'ward.D2',
               #cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               #breaks=bks,
               border_color = NA,
               color=hmcols,
               cellwidth=0.4,
               cellheight=5,
               annotation_colors = ann_colors,
               annotation_row= annotation_row
               #annotation_col=annotation_col
                        )
pdf(paste0(figdir, '/degs_seurat_summary22.pdf'),width=10)
ph2
dev.off()

                           
plot_list=list()
plot_list[[1]] = ph1[[4]]                            
plot_list[[2]] = ph2[[4]]                          
#g<-do.call(gridExtra::grid.arrange,plot_list)
                           
g <- gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs= plot_list,ncol=2))
                           
ggsave(paste0(figdir, '/degs_seurat_summary_up.pdf'),g,width =15,height = 10)                    
#######

                 
                           
                        
                      
#########################################
#down genes

cmarker_list_test <- lapply(X = scelist, FUN = function(x) {
    x <- FindAllMarkers(x)
    ####filter based \logFC\ > 0.25 & adj pvalue < 0.05
    x <- subset(x, avg_log2FC < -0.25 & p_val_adj < 0.05 )
    
    x <- data.frame(x, group= 'down')
    
    x <- x[, 6:8]
    
    #library(tidyr)
    x <- spread(x, key = cluster, value = group)
    #x= as.data.frame(x)
})

#####
cmlist= list()
for(i in 1:length(cmarker_list_test)) {
    tmd= cmarker_list_test[i]
    gn= names(cmarker_list_test[i])
    aa= as.data.frame(tmd, check.names=F)
    colnames(aa)[1]='gname'
    cmlist[[i]]= aa
    
}

df_merge  <- Reduce(function(...) merge(..., by="gname", all=TRUE), 
                    cmlist)
rownames(df_merge)=df_merge$gname
df_merge=df_merge[,-1]
head(df_merge)

all_groupName= unique(paste(sce$Defined_Cell_Subtype_function, sce$age_group, sep='.'))
length(all_groupName)

missedName= all_groupName[!all_groupName %in% colnames(df_merge)]

df_miss= data.frame(matrix(NA, nrow = nrow(df_merge), ncol =  length(missedName)),row.names= rownames(df_merge))
colnames(df_miss)= missedName
head(df_miss)

df_merge= cbind(df_merge, df_miss)
head(df_merge)
write.csv(df_merge, file= paste0(datadir, 'significant_down_degs_aged_cellSubtype_seurat.merged.csv'))
########################
df_merge=read.csv(paste0(datadir, 'significant_down_degs_aged_cellSubtype_seurat.merged.csv'), header=T, row.names=1, check.names=F)

####summarize the overlaps
#cn= names(cmarker_list_test)
cn= unique(sce$Defined_Cell_Subtype_function)
cnlist=list()
for(i in 1: length(cn)){
    tmd=df_merge[, grep(cn[i], colnames(df_merge))]
    head(tmd)
    nc1= sjmisc::row_count(tmd, count=NA, var=cn[i], append=FALSE)
    nc1[,paste(cn[i], 'state', sep='_')]= ifelse(nc1[,1] == ncol(tmd), 'not_sig', 'sig')
    #head(nc1)
    cnlist[[i]]= nc1
}
ctype_count= do.call(cbind, cnlist)
ctype_count= ctype_count[, grepl('state', colnames(ctype_count))]
head(ctype_count)
ctype_count_group= sjmisc::row_count(ctype_count, count='sig', var='ctype_count_group', append=FALSE)
#
down_count_group = sjmisc::row_count(df_merge, count='down', var='down_count', append=FALSE)
#down_count_group = sjmisc::row_count(df_merge, count='down', var='down_count', append=FALSE)
#
#stat_count= data.frame(ctype_count_group=ctype_count_group,up_count_group=up_count_group, down_count_group= down_count_group )
stat_count= data.frame(ctype_count_group=ctype_count_group,down_count_group=down_count_group)


stat_count$group='type3'
stat_count$group=ifelse(stat_count$ctype_count_group >=2,'type1', ifelse(stat_count$ctype_count_group ==1 , 'type2', 'type3'))

write.csv(stat_count, file= paste0(datadir, 'significant_down_degs_aged_cellSubtype_seurat.merged_stat_count.csv'))

head(stat_count)
dim(stat_count)
######
######
########
library(RColorBrewer)
library(pheatmap)
library(viridis)
#hmcols=colorRampPalette(rev(brewer.pal(10, "RdBu")))(256)
hmcols=c(colorRampPalette(c("white", rev(viridis(323, begin = 0.15))[1]))(10), rev(viridis(323, begin = 0.18)))
hmcols<-colorRampPalette(c("white","white",'blue'))(100)
#######
annotation_row= data.frame(tmp= colnames(df_merge))
rownames(annotation_row)= annotation_row$tmp
annotation_row$cell_type= sapply(stringr::str_split(annotation_row$tmp, "[.]"), `[`, 1)
annotation_row$age= sapply(stringr::str_split(annotation_row$tmp, "[.]"), `[`, 2)
annotation_row$cluster= sapply(stringr::str_split(annotation_row$tmp, "[-]"), `[`, 1) 
annotation_row$cluster= sapply(stringr::str_split(annotation_row$cluster, "[C]"), `[`, 2)
#annotation_row$cluster= gsub('-','',  annotation_row$cluster) 
annotation_row$cluster=as.numeric(paste(annotation_row$cluster))                        
annotation_row=annotation_row[order(annotation_row$cluster, decreasing=F),c('cell_type','age')] 
                           
head(annotation_row)
#######


ann_colors  = list(age =c('Young'='#47EFFA','Intermediated'='#EE2EE8','Aged'='#F0E716'),
cell_type= c('C1-CD4-CCR7'=color_CLASS[1],
             'C2-CD4-HSPA1A'=color_CLASS[2],
             'C3-CD4-CXCR6'=color_CLASS[3],
             'C4-CD4-CXCL13'=color_CLASS[4],
             'C5-CD4-FOXP3'=color_CLASS[5],
             'C6-CD4-RORC'=color_CLASS[6],
             'C7-CD8-SLC4A10'=color_CLASS[7],
             'C8-CD8-GZMK'=color_CLASS[8],
             'C9-CD8-CX3CR1'=color_CLASS[9],
             'C10-CD8-ZNF683'=color_CLASS[10],
             'C11-CD8-LAYN'=color_CLASS[11],
             'C12-CD8-TRDV1'=color_CLASS[12],
             'C13-NK-FCGR3A'=color_CLASS[13],
             'C14-NK-XCL1'=color_CLASS[14],
             'C15-CD4-CD8-ISG15'=color_CLASS[15],
             'C16-CD4-CD8-MKI67'=color_CLASS[16])
             )
#
tt= stat_count[order(stat_count$ctype_count_group, decreasing=T),]

tt_sel= subset(tt, group=='type1')

heatmap_matrix=df_merge[rownames(tt_sel),rownames(annotation_row)]
heatmap_matrix[is.na(heatmap_matrix)] <- 0                    
heatmap_matrix[heatmap_matrix=="down"]<-1#

                           
#heatmap_matrix[] <- lapply(heatmap_matrix, function(x) as.numeric(as.character(x)))
                           
#annotation_col= tt_sel

   
ph1 <- pheatmap::pheatmap(as.matrix(t(heatmap_matrix)), 
	           main= paste0('No.of.degs:', dim(heatmap_matrix)[1],'(shared by at least two cell type)'),
               #useRaster = T,
                treeheight_col=0,
               cluster_cols=T, 
               cluster_rows=F, 
               show_rownames=F, 
               show_colnames=F, 
               fontsize=5,
                          legend=F,
                          annotation_legend=F,
               #clustering_distance_rows=row_dist,
               #clustering_method = 'ward.D2',
               #cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               #breaks=bks,
               border_color = NA,
               color=hmcols,
               cellwidth=0.4,
               cellheight=5,
               annotation_colors = ann_colors,
               annotation_row= annotation_row
               #annotation_col=annotation_col
                        )
pdf(paste0(figdir, '/degs_seurat_summary11_down.pdf'),width=10)
ph1
dev.off()
                                               
#
                       
                           
tt_sel= subset(tt, group=='type2')

heatmap_matrix=df_merge[rownames(tt_sel),rownames(annotation_row)]
heatmap_matrix[is.na(heatmap_matrix)] <- 0                    
heatmap_matrix[heatmap_matrix=="down"]<-1#

#heatmap_matrix[] <- lapply(heatmap_matrix, function(x) as.numeric(as.character(x)))
                           
ph2 <- pheatmap::pheatmap(as.matrix(t(heatmap_matrix)), 
	           main= paste0('No.of.degs:', dim(heatmap_matrix)[1],'(unique in one group in each cell type)'),
               #useRaster = T,
               treeheight_col=0,
               cluster_cols=T, 
               cluster_rows=F, 
               show_rownames=T, 
               show_colnames=F, 
               fontsize=5,
               #clustering_distance_rows=row_dist,
               #clustering_method = 'ward.D2',
               #cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               #breaks=bks,
               border_color = NA,
               color=hmcols,
               cellwidth=0.4,
               cellheight=5,
               annotation_colors = ann_colors,
               annotation_row= annotation_row
               #annotation_col=annotation_col
                        )
pdf(paste0(figdir, '/degs_seurat_summary22_down.pdf'),width=10)
ph2
dev.off()

                           
plot_list=list()
plot_list[[1]] = ph1[[4]]                            
plot_list[[2]] = ph2[[4]]                          
#g<-do.call(gridExtra::grid.arrange,plot_list)
                           
g <- gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs= plot_list,ncol=2))
                           
ggsave(paste0(figdir, '/degs_seurat_summary_down.pdf'),g,width =15,height = 10)                    
#######                           
##############################

########################################################################################################################################################################
#########degs only between young and old groups for each cell type
datadir='output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/aged_degs_seurat/'
if(!dir.exists(datadir)) dir.create(datadir,recursive=TRUE)
figdir='output2_step2.5.3.filter_junk_scanpy_pipeline/figure_Tcell_v2/aged_degs_seurat/'
if(!dir.exists(figdir)) dir.create(figdir,recursive=TRUE)

library(Seurat)
sce = readRDS(paste0('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/', 'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd_with_rawcounts.rds'))
sce
#
sce <- AddMetaData(
  object = sce,
  metadata = gsub('/', '-', sce$Defined_Cell_Subtype_function),
  col.name = 'Defined_Cell_Subtype_function'
)
Idents(sce)='age_group'
levels(sce)
####exclude Intermediated group
sce= subset(sce, idents=c('Young','Aged' ))
sce
####
scelist= SplitObject(sce, split.by = "Defined_Cell_Subtype_function")


cmarker_list <- lapply(X = scelist, FUN = function(x) {
    x <- FindAllMarkers(x)
    ####filter based \logFC\ > 0.25 & adj pvalue < 0.05
    x <- subset(x, abs(avg_log2FC) > 0.25 & p_val_adj < 0.05 )
    
})

cName= names(cmarker_list)
for (i in 1:length(cName)) {    
    aa= as.data.frame(cmarker_list[cName[i]])
    
    write.csv(aa,  paste0(datadir, cName[i],'_significant_degs_aged_cellSubtype_only_Young_vs_Old.csv'))
}
save(cmarker_list, file=paste0(datadir, 'significant_degs_aged_cellSubtype_seurat_only_Young_vs_Old.RData'))

###
###
###################################################################################################################################################use logFC to draw annotated with aging genes
#2021-05-17
l#oad(paste0(datadir, 'significant_degs_aged_cellSubtype_seurat.RData'))
load(paste0(datadir, 'significant_degs_aged_cellSubtype_seurat_only_Young_vs_Old.RData'))
###########

cmarker_list_test <- lapply(X = cmarker_list, FUN = function(x) {
    x=subset(x, cluster == 'Aged')
    x=x[,c('avg_log2FC','gene')]
    #library(tidyr)
    #x <- spread(x, key = cluster, value = group)
    x= as.data.frame(x)
})
#
cmlist= list()
for(i in 1:length(cmarker_list_test)) {
    tmd= cmarker_list_test[i]
    gn= names(cmarker_list_test[i])
    aa= as.data.frame(tmd, check.names=F)
    head(aa)
    colnames(aa)[2]='gname'
    rownames(aa)= aa$gname
    cmlist[[i]]= aa
    
}
#
df_merge  <- Reduce(function(...) merge(..., by="gname", all=TRUE), 
                    cmlist)

rownames(df_merge)=df_merge$gname
df_merge=df_merge[,-1]
head(df_merge)
###
#############
########
library(RColorBrewer)
library(pheatmap)
library(viridis)
#hmcols=colorRampPalette(rev(brewer.pal(10, "RdBu")))(256)
hmcols=c(colorRampPalette(c("white", rev(viridis(323, begin = 0.15))[1]))(10), rev(viridis(323, begin = 0.18)))

#######

#
annotation_row= data.frame(tmp= colnames(df_merge))
rownames(annotation_row)= annotation_row$tmp
annotation_row$cell_type= sapply(stringr::str_split(annotation_row$tmp, "[.]"), `[`, 1)
annotation_row$age= sapply(stringr::str_split(annotation_row$tmp, "[.]"), `[`, 2)
annotation_row$cluster= sapply(stringr::str_split(annotation_row$tmp, "[-]"), `[`, 1) 
annotation_row$cluster= sapply(stringr::str_split(annotation_row$cluster, "[C]"), `[`, 2)
#annotation_row$cluster= gsub('-','',  annotation_row$cluster) 
annotation_row$cluster=as.numeric(paste(annotation_row$cluster))                        
annotation_row=annotation_row[order(annotation_row$cluster, decreasing=F),c('cell_type','age')] 
                           
head(annotation_row)
#
gage= read.csv('/data/Zhuxq/young_LC_analysis/gseaDB/genage_human.csv', header=T)
head(gage)
#
ann_colors  = list(#age =c('Young'='#47EFFA','Intermediated'='#EE2EE8','Aged'='#F0E716'),
                   state= c('non_Aging_gene'= '#F0E716', 'Aging_gene'= '#EE2EE8'),
cell_type= c('C1-CD4-CCR7'=color_CLASS[1],'C2-CD4-HSPA1A'=color_CLASS[2],'C3-CD4-CXCR6'=color_CLASS[3],
             'C4-CD4-CXCL13'=color_CLASS[4],'C5-CD4-FOXP3'=color_CLASS[5],'C6-CD4-RORC'=color_CLASS[6],
             'C7-CD8-SLC4A10'=color_CLASS[7],'C8-CD8-GZMK'=color_CLASS[8],'C9-CD8-CX3CR1'=color_CLASS[9],
             'C10-CD8-ZNF683'=color_CLASS[10],'C11-CD8-LAYN'=color_CLASS[11],'C12-CD8-TRDV1'=color_CLASS[12],
             'C13-NK-FCGR3A'=color_CLASS[13],'C14-NK-XCL1'=color_CLASS[14],'C15-CD4-CD8-ISG15'=color_CLASS[15],
             'C16-CD4-CD8-MKI67'=color_CLASS[16])
             )
#
gann= data.frame(state= ifelse(rownames(df_merge) %in% gage$symbol, 'Aging_gene','non_Aging_gene'))
rownames(gann)=rownames(df_merge)
head(gann)
##########
########
library(RColorBrewer)
library(pheatmap)
library(viridis)
hmcols=colorRampPalette(rev(brewer.pal(10, "RdBu")))(256)
#hmcols=c(colorRampPalette(c("white", rev(plasma(323, begin = 0.15))[1]))(10), rev(plasma(323, begin = 0.18)))

#heatmap_matrix[] <- lapply(heatmap_matrix, function(x) as.numeric(as.character(x)))

annotation_col= gann
heatmap_matrix= df_merge[, rownames(annotation_row)]
heatmap_matrix[is.na(heatmap_matrix)] <- 0

############
#
#
gene_name=rownames(subset(annotation_col, state=='Aging_gene'))

library(grid)


paletteLength <- 100
#c1=viridis(24)[1]
#c2=viridis(24)[10]
#mycolor= colorRampPalette(c(c1,"white",c2))(100)
mycolor<-colorRampPalette(c("blue","white","red"))(100)

mycolor
myBreaks <- unique(c(seq(min(heatmap_matrix), 0, length.out=ceiling(paletteLength/2) + 1), 
                     seq(max(heatmap_matrix)/paletteLength, max(heatmap_matrix),
                         length.out=floor(paletteLength/2))))

#
source('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/useful_R_Fun.pheatmap.add.flag.R')


heatmap_matrix_lim=heatmap_matrix
heatmap_matrix_lim[heatmap_matrix_lim< -2] <- -2
heatmap_matrix_lim[heatmap_matrix_lim>2] <- 2

myBreaks <- unique(c(seq(min(heatmap_matrix_lim), 0, length.out=ceiling(paletteLength/2) + 1), 
                     seq(max(heatmap_matrix_lim)/paletteLength, max(heatmap_matrix_lim),
                         length.out=floor(paletteLength/2))))


ph1 <- pheatmap::pheatmap(as.matrix((heatmap_matrix_lim)), 
	           main= paste0('No.of.degs:', dim(heatmap_matrix)[1],' in Aged_vs_Young group', '\n','No.of.degs overlapped with GenAge database (n=307):',length(gene_name)),
               #useRaster = T,
               treewidth_col=20,
               cluster_cols=F, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               fontsize=5,
                          legend=T,
                          annotation_legend=T,
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               breaks=myBreaks,
               border_color = NA,
               color=mycolor,
               cellwidth=8,
               cellheight=0.4,
               na_col = "#DDDDDD",
               annotation_colors = ann_colors,
               annotation_row= annotation_col,
               annotation_col=annotation_row
                        )

pdf(paste0(figdir, '/degs_seurat_summary_sig_in_Aged_group.logFC_vertical_2.pdf'),height=10)

add.flag(ph1,
         kept.labels = gene_name,
         repel.degree = 0.2)
dev.off()

#
####count the sig down and up genes in each cell type
###
count_group = data.frame(up=colSums(heatmap_matrix>0), down= colSums(heatmap_matrix<0))
#rownames(count_group)= rownames(heatmap_matrix)
count_group= as.data.frame(t(count_group))
count_group= data.frame(group= rownames(count_group), count_group)
head(count_group)
#
library(tidyr)
dlong <- count_group %>% pivot_longer(cols=colnames(count_group)[-1],
                                                                     names_to= "cell_type",
                                                                values_to = "count")
dlong$cell_type= gsub('.avg_log2FC','',dlong$cell_type)

dlong$cell_type= gsub('[.]','-',dlong$cell_type)
dlong$cell_type=factor(dlong$cell_type, levels=rev(c('C1-CD4-CCR7','C2-CD4-HSPA1A','C3-CD4-CXCR6',
             'C4-CD4-CXCL13','C5-CD4-FOXP3','C6-CD4-RORC',
             'C7-CD8-SLC4A10','C8-CD8-GZMK','C9-CD8-CX3CR1',
             'C10-CD8-ZNF683','C11-CD8-LAYN','C12-CD8-TRDV1',
             'C13-NK-FCGR3A','C14-NK-XCL1','C15-CD4-CD8-ISG15',
             'C16-CD4-CD8-MKI67'
             )))
dlong 

library(ggplot2)

pdf(paste0(figdir, '/proportion.sig_in_Aged_vs_Young_updown_number_genes.pdf'),width=8,height = 4)

ggplot(dlong, aes(x=cell_type,y= count, fill=group)) + 
geom_bar(position="stack", stat="identity")+scale_fill_manual(values = c('blue','red'))+
 theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())+
        coord_flip()+
ylab('No. DEGs')
dev.off()
##################
#################################################################################################################
#################################################################################################################
#############################################################################################################################################
#2021-05-18
#start to use monocle smooth spline method to get degs with aging
#need to use raw counts, find genes in each cell type with aging
#
#load seurat object data and find markers in each age sample 
###
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')
#
library(Seurat)
sce = readRDS(paste0('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/', 'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd_with_rawcounts.rds'))
sce
#
datadir='output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/aged_degs_smooth/'
if(!dir.exists(datadir)) dir.create(datadir,recursive=TRUE)
figdir='output2_step2.5.3.filter_junk_scanpy_pipeline/figure_Tcell_v2/aged_degs_smooth/'
if(!dir.exists(figdir)) dir.create(figdir,recursive=TRUE)
#
######
sce <- AddMetaData(
  object = sce,
  metadata = gsub('/', '-', sce$Defined_Cell_Subtype_function),
  col.name = 'Defined_Cell_Subtype_function'
)
####
#####*******delete ribosome and MT genes
RPS.genes <- grep(pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA", x = rownames(sce), value = TRUE)
MT.genes=  grep(pattern = "^MT-", x = rownames(sce), value = TRUE)
congene= rownames(sce)[!(rownames(sce) %in% c(RPS.genes, MT.genes))]
length(congene)
#
sce=subset(sce, features= congene)
Idents(sce)='age'
levels(sce)
sce
#######
scelist= SplitObject(sce, split.by = "Defined_Cell_Subtype_function")
####
table(sce$Defined_Cell_Subtype_function, sce$age)
####marker list in each cell type with aging 
cmarker_list <- lapply(X = scelist, FUN = function(x) {
    x <- FindAllMarkers(x,only.pos = F, min.pct = 0.25, logfc.threshold = 0.2)
    ####filter based \logFC\ > 0.25 & adj pvalue < 0.05
    x <- subset(x, abs(avg_log2FC) > 0.25 & p_val_adj < 0.05 )
    
})
#######save function
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


cname= names(cmarker_list)

for(i in 1:length(cname)){
    #sig gene name
    agedeglist= unique(cmarker_list[[i]]$gene)
    length(agedeglist)
   #get the matrix
   sce_sel= scelist[[i]]
   pt.matrix= as.matrix(sce_sel@assays$RNA@counts[agedeglist,order(sce_sel$intage, decreasing = F)])
    aa= pt.matrix
   dim(pt.matrix)
    ##https://github.com/crickbabs/ZebrafishDevelopingHindbrainAtlas
    pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
    pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)})) 
    rownames(pt.matrix)=agedeglist
    colnames(pt.matrix)= colnames(aa)
    #
    #####order as age increasing
    heatmap_matrix=pt.matrix
norder= sce_sel@meta.data[order(sce_sel@meta.data$intage, decreasing = F), ]
head(norder)
heatmap_matrix= as.data.frame(t(heatmap_matrix))[rownames(norder), ]
heatmap_matrix=t(heatmap_matrix)
dim(heatmap_matrix)
####clustering
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1
####color
bks <- seq(-3.1,3.1, by = 0.1)
hmcols <- colorRamps::blue2green2red(length(bks) - 1)
####
num_clusters <- 6

    
ph <- pheatmap::pheatmap(heatmap_matrix, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=F, 
               show_colnames=F, 
               clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols)

annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))	

write.csv(heatmap_matrix, file=paste0(datadir,cname[i],'_heatmap_matrix.csv'))
write.csv(annotation_row, file=paste0(datadir,cname[i],'_annotation_row_clusters.csv'))


annotation_col= data.frame(age=norder$intage )	
rownames(annotation_col)= rownames(norder)	
head(annotation_col)	
#####
ann_colors = list(
    #age = c("#edf8b1", "#081d58")
	age = rev(viridis::viridis(100))
)

ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= cname[i],
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=F, 
               show_colnames=F, 
               clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
  annotation_colors = ann_colors,
  annotation_row= annotation_row, annotation_col=annotation_col)

    
save_pheatmap_pdf(ph, paste0(figdir, cname[i], '_age_pheatmap_total_clusters.pdf'))

save(heatmap_matrix,annotation_row, annotation_col,ann_colors,cmarker_list,cname,
     file= paste0(datadir,cname[i],'_data_for_subheatmap_draw.RData') )
#####################
}
##################
##################chcek the output and confirm the clsuter to extract
datadir='output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/aged_degs_smooth/'
#if(!dir.exists(datadir)) dir.create(datadir,recursive=TRUE)
figdir='output2_step2.5.3.filter_junk_scanpy_pipeline/figure_Tcell_v2/aged_degs_smooth/'
#if(!dir.exists(figdir)) dir.create(figdir,recursive=TRUE)

cname
up_cluster= c(2,2,1,1,2,2,1,1,1,3,1,1,3,1,1,1)
down_cluster=c(3,4,5,5,3,6,2,3,4,2,4,3,5,3,5,3)
plot_list=list()
for(i in 1:length(cname)){
    load(paste0(datadir,cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster==up_cluster[i] | Cluster==down_cluster[i]))
    lenupg= length(rownames(subset(annotation_row, Cluster==up_cluster[i] )))
    lendowng= length(rownames(subset(annotation_row, Cluster==down_cluster[i] )))           
    heatmap_matrix= heatmap_matrix[gsel,]
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'No.genes.up:',up_cluster[i],':',lenupg,';','No.genes.down:',down_cluster[i],':',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=F, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
}
                  
g <- gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs= plot_list,ncol=4))
library(ggplot2)    
#save_pheatmap_pdf(g, paste0(figdir, '/degs_age_dependent_smooth_merged.pdf'),width =16,height = 16)
ggsave(paste0(figdir, '/degs_age_dependent_smooth_merged2.pdf'),g,width =16,height = 20) 
####################################################
###########age smooth gene enrichment
###########
##hallMarker and c2cp
library(clusterProfiler)
c5 <- read.gmt('/data/Zhuxq/young_LC_analysis/gseaDB/h.all.v7.4.entrez.gmt')
c2cp <- read.gmt('/data/Zhuxq/young_LC_analysis/gseaDB/c2.cp.v7.4.entrez.gmt')

cname
up_cluster= c(2,2,1,1,2,2,1,1,1,3,1,1,3,1,1,1)
down_cluster=c(3,4,5,5,3,6,2,3,4,2,4,3,5,3,5,3)

enrich_list=list()
for(i in 1:length(cname)){
    load(paste0(datadir,cname[i],'_data_for_subheatmap_draw.RData'))
    gseup= rownames(subset(annotation_row, Cluster==up_cluster[i] ))
    gsedown= rownames(subset(annotation_row, Cluster==down_cluster[i]))
    #lenupg= length(rownames(subset(annotation_row, Cluster==up_cluster[i] )))
    #lendowng= length(rownames(subset(annotation_row, Cluster==down_cluster[i] ))) 
    
    ######enrichment for up genes
    genes <- gseup
    genes = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    
    ec5 <- enricher(genes$ENTREZID, TERM2GENE=c5,minGSSize = 3,pvalueCutoff = 0.05)
    #head(ec5)
    ec2cp <- enricher(genes$ENTREZID, TERM2GENE=c2cp,minGSSize = 3,pvalueCutoff = 0.05)
    #head(ec2cp) 
    ######
    up_enrich= rbind(as.data.frame(ec5), as.data.frame(ec2cp))
    #
    if (dim(up_enrich)[1]==0){
        aa=  data.frame(matrix(NA, nrow = 1, ncol = dim(up_enrich)[2]) )
        colnames(aa)= colnames(up_enrich)
        up_enrich=aa
    } 
    up_enrich$Group= 'age_Young_up'
    head(up_enrich)
    ####
    ####
    ######enrichment for up genes
    genes <- gsedown
    genes = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    
    ec5 <- enricher(genes$ENTREZID, TERM2GENE=c5,minGSSize = 3,pvalueCutoff = 0.05)
    #head(ec5)
    ec2cp <- enricher(genes$ENTREZID, TERM2GENE=c2cp,minGSSize = 3,pvalueCutoff = 0.05)
    #head(ec2cp) 
    ######
    down_enrich= rbind(as.data.frame(ec5), as.data.frame(ec2cp))
    
    if (dim(down_enrich)[1]==0){
        aa=  data.frame(matrix(NA, nrow = 1, ncol = dim(down_enrich)[2]))
        colnames(aa)= colnames(down_enrich)
        down_enrich= aa
    } 
    #
    down_enrich$Group= 'age_Old_up'
    head(down_enrich)
    ######
    enrich_out= rbind(up_enrich, down_enrich)
    enrich_out$celltype= cname[i]
    
    enrich_list[[i]]= enrich_out
}
#
enrich_merged= do.call(rbind, enrich_list)
save(enrich_merged, file= paste0(datadir,'enrichment_results_merged_for_each_celltype.RData'))
write.csv(enrich_merged, file= paste0(datadir,'enrichment_results_merged_for_each_celltype.csv'))
#####
#####if non enrichment was found fro some group, use metascape maybe
#####

#####
###########################################################################################################################################################
###########################################################################################################################################################
#2021-05-19
#################aging dependent TFs
#################
################load data
#library(Seurat)
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output4_step4.0.scenic.allcells.py/output_Tcell')

figdir <- file.path("figure_Tcell/")
if(!dir.exists(figdir)) dir.create(figdir,recursive=TRUE)

datadir <- file.path("data_Tcell/")
if(!dir.exists(datadir)) dir.create(datadir,recursive=TRUE)

#####
#load seurat data to get cell type information
#
scdir='/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/'
library(Seurat)
sce=readRDS(paste0(scdir,'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd_with_rawcounts.rds'))
sce
#
cname=  gsub('/','-',unique(sce$Defined_Cell_Subtype_function))
cname= c( "C1-CD4-CCR7" , "C2-CD4-HSPA1A" ,"C3-CD4-CXCR6", "C4-CD4-CXCL13" ,"C5-CD4-FOXP3",
           "C6-CD4-RORC" , "C7-CD8-SLC4A10" ,"C8-CD8-GZMK" ,"C9-CD8-CX3CR1" ,"C10-CD8-ZNF683"  ,
		   "C11-CD8-LAYN" , "C12-CD8-TRDV1",    "C13-NK-FCGR3A", "C14-NK-XCL1" ,"C15-CD4-CD8-ISG15" , "C16-CD4-CD8-MKI67")
#
source('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/useful_R_Fun.pheatmap.add.flag.R')
#
gage= read.csv('/data/Zhuxq/young_LC_analysis/gseaDB/genage_human.csv', header=T)
head(gage)
#
plot_list=list()
#
for (i in 1:length(cname)){
    #read the rss data
    auc_mtx= read.csv(paste0(cname[i],'_auc_sel.csv'),header=T, row.names=1, check.names=F)
    rss= read.csv(paste0(cname[i],'_rss.csv'),header=T, row.names=1, check.names=F)
    adata.obs=read.csv(paste0(cname[i],'_adata.obs.csv'),header=T, row.names=1, check.names=F)
    #set the order with age increasing
    adata.obs$age_num=as.numeric(paste(adata.obs$age))
    adata.obs=adata.obs[order(adata.obs$age_num, decreasing=F), ]

    #####select age specific regulons
    colnames(rss)= gsub('[(+)]', '', colnames(rss))
    rss_t= as.data.frame(t(rss))
    colnames(rss_t)= paste('age', colnames(rss_t), sep='_')
    rss_t= data.frame(gname= rownames(rss_t), rss_t)
    rss_t[1:3,1:4]
    #
    library(tidyr)
    library(dplyr)
    datalong <- rss_t %>% pivot_longer(cols=colnames(rss_t)[-1],names_to= "group",values_to = "auc")
    head(datalong)

    datalong2 <- datalong %>%   
    arrange(desc(auc)) %>% 
    group_by(group) %>%
    slice(1:20)
    datalong2
    #
    tflist= unique(datalong2$gname)
    length(tflist)
    #####prepare for heatmap
    all(rownames(adata.obs) %in% rownames(auc_mtx))
    colnames(auc_mtx)=gsub('[(+)]', '', colnames(auc_mtx))
    auc_mtx= as.data.frame(t(auc_mtx[rownames(adata.obs),]))

    auc_mtx[1:4,1:4]
    pt.matrix = auc_mtx[tflist, ]
    dim(pt.matrix)
    ##
    aa= pt.matrix
    ##https://github.com/crickbabs/ZebrafishDevelopingHindbrainAtlas
    pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
    pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)})) 
    rownames(pt.matrix)=tflist
    colnames(pt.matrix)= colnames(aa)
    #
    #####order as age increasing
    heatmap_matrix=pt.matrix
    ####clustering
    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
    row_dist[is.na(row_dist)] <- 1
    ####color
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    ####
    num_clusters <- 6

    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=F, 
               show_colnames=F, 
               clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols)

    annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))	

    write.csv(heatmap_matrix, file=paste0(datadir,cname[i],'_heatmap_matrix.csv'))
    write.csv(annotation_row, file=paste0(datadir,cname[i],'_annotation_row_clusters.csv'))

    norder= adata.obs
    annotation_col= data.frame(age=norder$age_num )	
    rownames(annotation_col)= rownames(norder)	
    head(annotation_col)	
    #####
    ann_colors = list(
    #age = c("#edf8b1", "#081d58")
	age = rev(viridis::viridis(100))
    )

    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= cname[i],
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=num_clusters,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
               annotation_colors = ann_colors,
               annotation_row= annotation_row, annotation_col=annotation_col)
    #
    plot_list[[i]] = ph[[4]]
    #
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    library(grid)
    
    #save_pheatmap_pdf(ph, paste0(figdir, cname[i], '_age_pheatmap_total_clusters.pdf'))
    pdf(paste0(figdir, cname[i], '_age_pheatmap_total_clusters.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()
    
    
    save(heatmap_matrix,annotation_row, annotation_col,ann_colors,datalong2,tflist,cname,gene_name,
     file= paste0(datadir,cname[i],'_data_for_subheatmap_draw.RData') )
    #
}
#
g <- gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs= plot_list,ncol=4))
library(ggplot2)    
#save_pheatmap_pdf(g, paste0(figdir, '/degs_age_dependent_smooth_merged.pdf'),width =16,height = 16)
ggsave(paste0(figdir, '/TFs_age_dependent_smooth_initial_clusters_merged.pdf'),g,width =16,height = 20) 

###############
#2021-05-20
##################chcek the output and confirm the clsuter to extract
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output4_step4.0.scenic.allcells.py/output_Tcell')

figdir <- file.path("figure_Tcell/reduced/")
if(!dir.exists(figdir)) dir.create(figdir,recursive=TRUE)

datadir <- file.path("data_Tcell/reduced/")
if(!dir.exists(datadir)) dir.create(datadir,recursive=TRUE)


cname= c( "C1-CD4-CCR7" , "C2-CD4-HSPA1A" ,"C3-CD4-CXCR6", "C4-CD4-CXCL13" ,"C5-CD4-FOXP3",
           "C6-CD4-RORC" , "C7-CD8-SLC4A10" ,"C8-CD8-GZMK" ,"C9-CD8-CX3CR1" ,"C10-CD8-ZNF683"  ,
		   "C11-CD8-LAYN" , "C12-CD8-TRDV1",    "C13-NK-FCGR3A", "C14-NK-XCL1" ,"C15-CD4-CD8-ISG15" , "C16-CD4-CD8-MKI67")
gage= read.csv('/data/Zhuxq/young_LC_analysis/gseaDB/genage_human.csv', header=T)
head(gage)
library(grid)
source('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/useful_R_Fun.pheatmap.add.flag.R')

plot_list=list()
tfover=list()
tflist_young= list()
tflist_old= list()
######C1-CD4-CCR7
up_cluster= c(2,6)
down_cluster=c(1,4)
i=1
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]

    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down
    
    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()
#    



######C2-CD4-HSPA1A
up_cluster= c(3,5)
down_cluster=c(2,1)
i=2
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down

    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()


######C3-CD4-CXCR6
up_cluster= c(4,3)
down_cluster=c(2,6)
i=3
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()


######C4-CD4-CXCL13
up_cluster= c(4,2)
down_cluster=c(1,3,5)
i=4
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()

###

######C5-CD4-FOXP3
up_cluster= c(3,1)
down_cluster=c(5,2)
i=5
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()

######C6-CD4-RORC
up_cluster= c(3,2)
down_cluster=c(1,6)
i=6
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()

######C7-CD8-SLC4A10
up_cluster= c(4,5)
down_cluster=c(3,6)
i=7
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()

##########C8-CD8-GZMK
up_cluster= c(1,2)
down_cluster=c(4,6)
i=8
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()


#########C9-CD8-CX3CR1
up_cluster= c(1,3)
down_cluster=c(2,6)
i=9
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()

#########C10-CD8-ZNF683
up_cluster= c(2,3)
down_cluster=c(1)
i=10
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()



#########C11-CD8-LAYN
up_cluster= c(2,1)
down_cluster=c(4,5,6)
i=11
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()



######C12-CD8-TRDV1
up_cluster= c(5,1)
down_cluster=c(3,4,6)
i=12
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()

####C13-NK-FCGR3A
up_cluster= c(6,3)
down_cluster=c(5,2)
i=13
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()


####C14-NK-XCL1
up_cluster= c(1,2,3)
down_cluster=c(5,4)
i=14
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()


####C15-CD4-CD8-ISG15
up_cluster= c(1,3)
down_cluster=c(2,4)
i=15
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()


#######
up_cluster= c(1,5)
down_cluster=c(3,4)
i=16
    load(paste0('data_Tcell/',cname[i],'_data_for_subheatmap_draw.RData'))
    gsel= rownames(subset(annotation_row, Cluster %in% c(up_cluster,down_cluster )))
    lenupg= length(rownames(subset(annotation_row, Cluster %in% up_cluster )))
    lendowng= length(rownames(subset(annotation_row, Cluster %in% down_cluster )))   
   
    heatmap_matrix= heatmap_matrix[gsel,]
    gene_name= intersect(rownames(heatmap_matrix), gage$symbol)
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    
    ph <- pheatmap::pheatmap(heatmap_matrix, 
               main= paste0(cname[i],'\n', 'TFs.up:',lenupg,';','TFs.down:',lendowng),
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               #clustering_distance_rows=row_dist,
               clustering_method = 'ward.D2',
               #cutree_rows=2,
               silent=TRUE,
               filename=NA,
               breaks=bks,
               border_color = NA,
               color=hmcols,
             annotation_colors = ann_colors,
             annotation_row= annotation_row, 
             annotation_col=annotation_col)
    plot_list[[i]] = ph[[4]]  
    
    tfover[[i]]= gene_name

    gene_up= intersect(rownames(subset(annotation_row, Cluster %in% up_cluster )), gage$symbol)
    gene_down= intersect(rownames(subset(annotation_row, Cluster %in% down_cluster )), gage$symbol)
   
    tflist_young[[i]]= gene_up
    tflist_old[[i]]= gene_down


    pdf(paste0(figdir, cname[i], '_age_pheatmap_clusters_reduced.pdf'))
    add.flag(ph,kept.labels = gene_name,repel.degree = 0.5)
    dev.off()
#####
names(tfover)= cname
names(tflist_young)=cname
names(tflist_old)=cname
save(tfover, tflist_young, tflist_old, file= paste0(datadir, '/TFs_lists_overlappedwithagedatabase_upanddownlist.RData'))

g <- gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs= plot_list,ncol=4))
library(ggplot2)    
#save_pheatmap_pdf(g, paste0(figdir, '/degs_age_dependent_smooth_merged.pdf'),width =16,height = 16)
ggsave(paste0(figdir, '/TFs_age_dependent_smooth_merged_reduced.pdf'),g,width =20,height = 20) 

#
######################
######################
#2021-06-29
######gene pathways enrichment analysis
######
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')
datadir='output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/'
figdir='output2_step2.5.3.filter_junk_scanpy_pipeline/figure_Tcell_v2/'


load(paste0(datadir, 'aged_degs_smooth/enrichment_results_merged_for_each_celltype.RData' ))
head(enrich_merged)
enrich_merged$celltype= factor(enrich_merged$celltype, levels=c( "C1-CD4-CCR7" , "C2-CD4-HSPA1A" ,"C3-CD4-CXCR6", "C4-CD4-CXCL13" ,"C5-CD4-FOXP3",
           "C6-CD4-RORC" , "C7-CD8-SLC4A10" ,"C8-CD8-GZMK" ,"C9-CD8-CX3CR1" ,"C10-CD8-ZNF683"  ,
		   "C11-CD8-LAYN" , "C12-CD8-TRDV1",    "C13-NK-FCGR3A", "C14-NK-XCL1" ,"C15-CD4-CD8-ISG15" , "C16-CD4-CD8-MKI67"))
#enrich_merged$logp=-log10(enrich_merged$pvalue)

enrich_merged_sel= subset(enrich_merged, Group=='age_Young_up')
freqdata= as.data.frame.matrix(table(enrich_merged_sel$ID, enrich_merged_sel$celltype))
freqdata$count= rowSums(freqdata)
freqdata=freqdata[order(freqdata$count, decreasing=T),]
write.csv(freqdata, paste0(datadir, 'aged_degs_smooth/enrichment_results_merged_for_each_celltype.frequency.stat_age_Young_up.csv' ))
##
enrich_merged_sel= subset(enrich_merged, Group=='age_Old_up')
freqdata= as.data.frame.matrix(table(enrich_merged_sel$ID, enrich_merged_sel$celltype))
freqdata$count= rowSums(freqdata)
freqdata=freqdata[order(freqdata$count, decreasing=T),]
write.csv(freqdata, paste0(datadir, 'aged_degs_smooth/enrichment_results_merged_for_each_celltype.frequency.stat_age_Old_up.csv' ))

#####
#####show selected pathways using bubble plot
#####
library(patchwork)
library(ggplot2)
library(ggpubr)
selpaths= read.csv(paste0(datadir, 'aged_degs_smooth/selected_pathways.csv' ), header=T)

selpaths$ID=gsub('KEGG_','',selpaths$paths)
selpaths$ID=gsub('REACTOME_','',selpaths$ID)
selpaths$ID=gsub('BIOCARTA_','',selpaths$ID)
selpaths$ID=gsub('HALLMARK_','',selpaths$ID)
selpaths$ID=gsub('PID_','',selpaths$ID)
selpaths$ID=gsub('WP_','',selpaths$ID)
selpaths$ID=tolower(selpaths$ID)

head(selpaths)
#####
age_Y_path= enrich_merged[enrich_merged$ID %in% subset(selpaths, group=='age_Young_up')$paths, ]
age_Y_path=subset(age_Y_path, Group=='age_Young_up')

age_Y_path$ID=gsub('KEGG_','',age_Y_path$ID)
age_Y_path$ID=gsub('REACTOME_','',age_Y_path$ID)
age_Y_path$ID=gsub('BIOCARTA_','',age_Y_path$ID)
age_Y_path$ID=gsub('HALLMARK_','',age_Y_path$ID)
age_Y_path$ID=gsub('PID_','',age_Y_path$ID)
age_Y_path$ID=gsub('WP_','',age_Y_path$ID)
age_Y_path$ID=tolower(age_Y_path$ID)
age_Y_path$ID= factor(age_Y_path$ID, levels=rev(subset(selpaths, group=='age_Young_up')$ID))
age_Y_path


#####
age_O_path= enrich_merged[enrich_merged$ID %in% subset(selpaths, group=='age_Old_up')$paths, ]
age_O_path=subset(age_O_path, Group=='age_Old_up')

age_O_path$ID=gsub('KEGG_','',age_O_path$ID)
age_O_path$ID=gsub('REACTOME_','',age_O_path$ID)
age_O_path$ID=gsub('BIOCARTA_','',age_O_path$ID)
age_O_path$ID=gsub('HALLMARK_','',age_O_path$ID)
age_O_path$ID=gsub('PID_','',age_O_path$ID)
age_O_path$ID=gsub('WP_','',age_O_path$ID)
age_O_path$ID=tolower(age_O_path$ID)
age_O_path$ID= factor(age_O_path$ID, levels=rev(subset(selpaths, group=='age_Old_up')$ID))
age_O_path

#
g1=ggplot(age_Y_path, aes(x = celltype, y = ID,color=p.adjust)) +
    geom_point(aes( size = Count), alpha = 1) +
    theme_classic2()+
    scale_colour_gradientn(limits=c(5-05,0.05), #colours = c("darkred", "orange", "yellow", "white"),
                           #colours = c('#EE2EE8','#47EFFA','#F0E716', "white"),
                           colours = c('#47EFFA', "white"),
                           breaks = c( 0.01, 0.03,0.05),
                           labels = c(0.01, 0.03,0.05))+
    xlab('cell subtype')+ylab('Age-dependent Young enriched')+
    labs(color="adjusted P value", size= 'Count')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())

g2=ggplot(age_O_path, aes(x = celltype, y = ID,color=p.adjust)) +
    geom_point(aes( size = Count), alpha = 1) +
    theme_classic2()+
    scale_colour_gradientn(limits=c(5-05,0.05), #colours = c("darkred", "orange", "yellow", "white"),
                           #colours = c('#F0E716','#47EFFA','#EE2EE8', "white"),
                           colours = c('#EE2EE8', "white"),
                           breaks = c( 0.01, 0.03,0.05),
                           labels = c(0.01, 0.03,0.05))+
    xlab('cell subtype')+ylab('Age-dependent Old enriched')+
    labs(color="adjusted P value", size= 'Count')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank())



#####
pdf(paste0(figdir, 'aged_degs_smooth/degs_enriched_pathways.pdf'),width = 18,height = 5)
g1+g2+plot_layout(ncol = 2)
dev.off()

###############################################################
###############################################################
#2021-10-26
#write out age-dependent gene list


setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')
#
datadir='output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/aged_degs_smooth/'
###
###jus load for cname list
#load(paste0(datadir,'C1-CD4-CCR7','_data_for_subheatmap_draw.RData'))
cname= c('C3-CD4-CXCR6','C11-CD8-LAYN','C8-CD8-GZMK','C15-CD4-CD8-ISG15','C6-CD4-RORC','C10-CD8-ZNF683','C1-CD4-CCR7','C9-CD8-CX3CR1','C4-CD4-CXCL13','C2-CD4-HSPA1A','C12-CD8-TRDV1','C13-NK-FCGR3A','C5-CD4-FOXP3','C7-CD8-SLC4A10','C14-NK-XCL1','C16-CD4-CD8-MKI67')
up_cluster= c(2,2,1,1,2,2,1,1,1,3,1,1,3,1,1,1)
down_cluster=c(3,4,5,5,3,6,2,3,4,2,4,3,5,3,5,3)

gene_list=list()
for(i in 1:length(cname)){
    load(paste0(datadir,cname[i],'_data_for_subheatmap_draw.RData'))
    age_dep_down= data.frame(geneName= rownames(subset(annotation_row, Cluster==up_cluster[i] )),
                            group= 'age_dependent_down')
    
    age_dep_up= data.frame(geneName= rownames(subset(annotation_row, Cluster==down_cluster[i])),
                          group= 'age_dependent_up')
    aa= rbind(age_dep_down, age_dep_up)
    aa$celltype= paste0(cname[i])
    #
    gene_list[[i]]= aa
}
#
gene_list_merged= do.call(rbind, gene_list)
write.csv(gene_list_merged, file= paste0(datadir,'Tcell_age_dependent_DEGs_list_merged.csv'))
#####


#######20220410
#######---------------------------------------------------------------
#######test for DirichletReg for T cells
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')
datadir='output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/'
figdir='output2_step2.5.3.filter_junk_scanpy_pipeline/figure_Tcell_v2/'
#
library(Seurat)
sce= readRDS(paste0(datadir, 'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd_with_rawcounts.rds'))
#
stat1=sce@meta.data
#
x=as.data.frame.matrix(table(stat1$Defined_Cell_Subtype_function,stat1$Sample))
data_percentage <- apply(x, 2, function(x){x/sum(x,na.rm=T)})
data_percentage=as.data.frame(t(data_percentage))

ageinfo= stat1[!(duplicated(stat1$Sample)), c('Sample', 'intage')]
rownames(ageinfo)= ageinfo$Sample
ageinfo
#
prostat= data.frame('age'= ageinfo[rownames(data_percentage), 'intage'],data_percentage )
prostat
######
library("DirichletReg")
AL <- DR_data(prostat[, -c(1)])
#plot(AL, cex = 0.5, a2d = list(colored = FALSE, c.grid = FALSE))
lake1 <- DirichReg(AL ~ age, prostat)
lake1
coef(lake1)


lake2 <- update(lake1, . ~ . + I(age^2) | . + I(age^2) | . + I(age^2)|. + I(age^2)|
. + I(age^2) | . + I(age^2) | . + I(age^2)|. + I(age^2)|
. + I(age^2) | . + I(age^2) | . + I(age^2)|. + I(age^2)|
. + I(age^2) | . + I(age^2) | . + I(age^2)|. + I(age^2)
)

anova(lake1, lake2)
#
mycol= c('#0000FF', # 0
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
                 '#F0FFFF' # 35 
)

pdf(paste0(figdir, '/DirichletReg_Aged_vs_Defined_Cell_Subtype_function_rawmodel.pdf'),width=7,height = 9)
par(mar = c(4, 4, 4, 4) + 0.1)
plot(rep(prostat$age, 16), as.numeric(AL), pch = 21, bg = rep(mycol[1:dim(AL)[2]], each = dim(AL)[1]), xlab = "age (yr)", ylab = "Proportion",
       ylim = c(0, 0.4), main = "Tcell Composition in aging Lung Cancer by DirichletReg")

Xnew <- data.frame(age = seq(min(prostat$age), max(prostat$age),
                  length.out = 100))
 for (i in 1:dim(AL)[2]) lines(cbind(Xnew, predict(lake1, Xnew)[, i]), col = mycol[1:dim(AL)[2]][i], lwd = 2)
legend('topright', legend = colnames(AL), lwd = 2, 
       col =  mycol[1:dim(AL)[2]], pt.bg =  mycol[1:dim(AL)[2]], pch = 21,bty = "n")
par(new = TRUE)
plot(cbind(Xnew, predict(lake1, Xnew, F, F, T)), lty = "24", type = "l", ylim = c(0,
                max(predict(lake1, Xnew, F, F, T))), axes = F, ann = F, lwd = 2)
axis(4)
mtext(expression(paste("Precision (", phi, ")", sep = "")), 4, line = 3)
legend("top", legend = c(expression(hat(mu[c] == hat(alpha)[c]/hat(alpha)[0])),
                     expression(hat(phi) == hat(alpha)[0])), lty = c(1, 2), lwd = c(3, 2), bty = "n")
dev.off()
#
#draw coef
library("ggpubr")
d2= as.data.frame(t(as.data.frame(coef(lake1))))
d2$celltype= rownames(d2)
d2= d2[order(d2$age, decreasing=T),]

#d2$celltype= factor(d2$celltype, levels= c(
#'C1.CD4.CCR7','C2.CD4.HSPA1A','C3.CD4.CXCR6','C4.CD4.CXCL13','C5.CD4.FOXP3','C6.CD4.RORC','C7.CD8.SLC4A10','C8.CD8.GZMK', #'C9.CD8.CX3CR1','C10.CD8.ZNF683','C11.CD8.LAYN','C12.CD8.TRDV1','C13.NK.FCGR3A','C14.NK.XCL1','C15.CD4.CD8.ISG15','C16.CD4.CD8.MKI67')
#)

# Get p-values
u = summary(lake1)
pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
v = names(pvals)
pvals = matrix(pvals, ncol=length(u$varnames))
#rownames(pvals) = gsub('Health', '', v[1:nrow(pvals)])
rownames(pvals)='pvalue'
colnames(pvals) = u$varnames
pvals= as.data.frame(t(pvals))
pvals
###
d2= cbind(d2[rownames(pvals),], pvals )

d2$adjpvalue= p.adjust(as.vector(d2$pvalue), method='BH')
d2
####



pdf(paste0(figdir, '/DirichletReg_Aged_vs_Defined_Cell_Subtype_function_rawmodel_coef.pdf'),width=7,height = 5)

ggbarplot(d2, x = "celltype", y = "age",
          fill = "celltype",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          #color = "#525252", 
          palette = mycol,            # jco journal color palett. see ?ggpar
          #sort.val = "asc",          # Sort the value in dscending order
          sort.by.groups = TRUE,     # Don't sort inside each group
          x.text.angle = 90,           # Rotate vertically x axis texts
          ggtheme = theme_minimal()
)+coord_flip()+
#scale_y_continuous(breaks=c(0,1,3,6, 9,12),position = "right")+
  #scale_y_continuous(position = "right")+
  #font("x.text", size = 5, vjust = 0.5)+
  theme(axis.text.x = element_text(size=8))+
  ylab('coef of Dirichlet Regression')+ggtitle("T cell Composition")+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),legend.position="right",
        # Hide panel borders and remove grid lines
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank()
  )+scale_fill_manual(values=c(
  'C1.CD4.CCR7'=mycol[1],
      'C2.CD4.HSPA1A'=mycol[2],
      'C3.CD4.CXCR6'=mycol[3],
      'C4.CD4.CXCL13'=mycol[4],
      'C5.CD4.FOXP3'=mycol[5],
      'C6.CD4.RORC'=mycol[6],'C7.CD8.SLC4A10'=mycol[7],'C8.CD8.GZMK'=mycol[8], 'C9.CD8.CX3CR1'=mycol[9],'C10.CD8.ZNF683'=mycol[10],'C11.CD8.LAYN'=mycol[11],
      'C12.CD8.TRDV1'=mycol[12],'C13.NK.FCGR3A'=mycol[13],'C14.NK.XCL1'=mycol[14],
      'C15.CD4.CD8.ISG15'=mycol[15],'C16.CD4.CD8.MKI67'=mycol[16]))

dev.off()
####------------------------------
tcellcor=d2
colnames(tcellcor)[2]='cor'
minvalue= round(min(d2$cor), digits = 1)
maxvalue= round(max(d2$cor), digits = 1)
tcellcor$cor2=abs(tcellcor$cor)
tcellcor= tcellcor[order(tcellcor$cor, decreasing=T),]
tcellcor$celltype=  factor(tcellcor$celltype, levels= unique(tcellcor$celltype))

library(ggplot2)
library(ggpubr)
pdf(paste0(figdir, '/Cor_needle_age_vs_cellPercentage_DirichletReg.pdf'),width = 6,height = 5)
ggplot(tcellcor, aes(x = cor, y = celltype,color=adjpvalue)) +
    geom_segment(aes(yend = celltype), xend = 0, colour = "grey50") +
    #geom_point(size = Pvalue, aes(colour = Pvalue)) +
    geom_point(aes( size = cor2), alpha = 1) +
    #scale_colour_brewer(palette = "Set1", limits = c("significant", "non significant")) +
    theme_classic2()+
    scale_colour_gradientn(limits=c(5-05,0.05), colours = c("darkred", "orange", "yellow", "white"),
                           breaks = c( 0.01, 0.03,0.05),
                           labels = c(0.01, 0.03,0.05))+
    #scale_colour_gradient2(limits=c(5-05,0.05), low="darkred",  mid ="orange", high ="Black")
    xlab('coef of Dirichlet Regression')+ylab('Cluster')+
    labs(color="adj P value", size= 'Coef')+
    geom_vline(xintercept=0, linetype="dashed", color = "black")#+
    #scale_x_continuous(breaks=seq(minvalue, maxvalue, 0.01), limits = c(minvalue, maxvalue))
dev.off()


save(tcellcor, file = paste0(datadir, '/coef_pvalue_DirichletReg.analysis.csv'))

###########
#########20220414
#########-----------------------------------------------------------------
#### check the co-expression of ICB markers
datadir='output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/coexp_ExhaustM/'
if(!dir.exists(datadir)) dir.create(datadir,recursive=TRUE)
figdir='output2_step2.5.3.filter_junk_scanpy_pipeline/figure_Tcell_v2/coexp_ExhaustM/'
if(!dir.exists(figdir)) dir.create(figdir,recursive=TRUE)

library(Seurat)
sce = readRDS(paste0('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Tcell_v2/', 'adata_Tcell_v2_degs_signature_louvain_recluster_renamed_subtype_Defined_bad_cluster_rm_NaiveCytoExhNK_scoreAdd_with_rawcounts.rds'))
sce
#
sce <- AddMetaData(
  object = sce,
  metadata = gsub('/', '-', sce$Defined_Cell_Subtype_function),
  col.name = 'Defined_Cell_Subtype_function'
)
#
library(scater)
sce.sc=as.SingleCellExperiment(sce)
cpm(sce.sc) <- calculateCPM(sce.sc)
assay(sce.sc,"logcpm")= log(as.matrix(assay(sce.sc,"cpm"))+1)
#
expm_sel= as.data.frame(assay(sce.sc, 'logcpm'))[c('CD274','PDCD1','TIGIT','CTLA4'), ]
expm_sel[1:4,1:4]
library(dplyr)
stat_exp= cbind(colData(sce.sc), t(expm_sel)) %>% as.data.frame
head(stat_exp)
#
stat_exp$age_group= factor(stat_exp$age_group, levels= c('Young' ,'Intermediated' ,'Aged'))

####-----
stat_exp$Defined_Cell_Subtype_function= gsub('[-]', '_',  stat_exp$Defined_Cell_Subtype_function)
stat_exp$bioGroup= stat_exp$age_group

library(ggpubr)
#library(paletteer) 
#colgroup=as.vector(paletteer_d("unikn::pal_signal"))
#colgroup= rev(colgroup)

colgroup= c('#F0E716','#47EFFA','#EE2EE8')

my_color2=c("#7F3C8DFF", '#11A579FF',  '#3969ACFF',  '#F2B701FF',  '#E73F74FF',
           '#80BA5AFF', '#E68310FF', '#008695FF', '#CF1C90FF', '#F97B72FF', '#4B4B8FFF',  
    '#33A02CFF',  '#FB9A99FF',  '#E31A1CFF','#A5AA99FF','#B2DF8AFF', 
           '#FDBF6FFF',  '#FF7F00FF',  '#CAB2D6FF',  '#6A3D9AFF',  '#FFFF99FF',  '#B15928FF', '#A6CEE3FF',  '#1F78B4FF'

)


my_color2=c('#0000FF', # 0
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


cd4list= c('C3_CD4_CXCR6','C6_CD4_RORC',
           'C1_CD4_CCR7','C4_CD4_CXCL13','C2_CD4_HSPA1A','C5_CD4_FOXP3'
           )
cd8list= c('C11_CD8_LAYN','C8_CD8_GZMK','C10_CD8_ZNF683',
           'C9_CD8_CX3CR1','C12_CD8_TRDV1',
           'C7_CD8_SLC4A10')

#tmp= subset(stat_exp,  Defined_Cell_Subtype_function== 'C11_CD8_LAYN')
#tmp$TIGIT_bin= ifelse(tmp$TIGIT> 0, 'pos','neg')
#tmp$PDCD1_bin= ifelse(tmp$PDCD1> 0, 'pos','neg')

tmp= stat_exp[stat_exp$Defined_Cell_Subtype_function %in% cd4list, ]

cc= dim(subset(tmp, TIGIT >0  & CTLA4 >0 & bioGroup == 'Aged'))[1]/(table(tmp$bioGroup)['Aged']) * 100 
cc= as.numeric(cc) %>% round(digits = 1)
cc

bb= dim(subset(tmp, TIGIT >0  & CTLA4 >0 & bioGroup == 'Intermediated'))[1]/(table(tmp$bioGroup)['Intermediated']) * 100 
bb= as.numeric(bb) %>% round(digits = 1)
bb

aa= dim(subset(tmp, TIGIT >0  & CTLA4 >0 & bioGroup == 'Young'))[1]/(table(tmp$bioGroup)['Young']) * 100 
aa= as.numeric(aa) %>% round(digits = 1)
aa

g1=ggscatter(tmp, x = "TIGIT", y = "CTLA4",
          #add = "reg.line",                                 # Add regression line
          #conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black",fill = "lightgray"),
          color = "Defined_Cell_Subtype_function", palette =  my_color2[1:6],size= 2,
)+facet_grid(. ~ bioGroup)+
geom_point(color='lightgray', shape=21, size = 3, stroke=0.1, alpha= 0.3)+

annotate("text", x=3, y=5, label=paste0('Young: ',aa,'%','\n',
                                       'Intermediated: ',bb,'%','\n',
                                       'Aged: ',cc,'%','\n'),size=4)+

theme(legend.position = "right",
        #legend.justification = c(0,1),
        #legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        aspect.ratio = 1,
        plot.margin = margin(t=5,r=20,b=5,l=20,unit="pt"),
        axis.title.y = element_text(size=20,colour="black"),
        axis.title.x = element_text(size=20,colour="black"),
        #axis.title.x = element_blank(),
        axis.text = element_text(size=20,colour="black",angle=0),
        axis.text.x = element_text(size=20,colour="black",angle=0,hjust=0.5,vjust=0.5),
        #axis.text.x = element_blank(),
        plot.title = element_text(size=15,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(.1, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.5)
        )

#------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------

tmp= subset(stat_exp,  Defined_Cell_Subtype_function== 'C4_CD4_CXCL13')
#tmp$TIGIT_bin= ifelse(tmp$TIGIT> 0, 'pos','neg')
#tmp$PDCD1_bin= ifelse(tmp$PDCD1> 0, 'pos','neg')

tmp= stat_exp[stat_exp$Defined_Cell_Subtype_function %in% cd8list, ]


cc= dim(subset(tmp, TIGIT >0  & CTLA4 >0 & bioGroup == 'Aged'))[1]/(table(tmp$bioGroup)['Aged']) * 100 
cc= as.numeric(cc) %>% round(digits = 1)
cc

bb= dim(subset(tmp, TIGIT >0  & CTLA4 >0 & bioGroup == 'Intermediated'))[1]/(table(tmp$bioGroup)['Intermediated']) * 100 
bb= as.numeric(bb) %>% round(digits = 1)
bb

aa= dim(subset(tmp, TIGIT >0  & CTLA4 >0 & bioGroup == 'Young'))[1]/(table(tmp$bioGroup)['Young']) * 100 
aa= as.numeric(aa) %>% round(digits = 1)
aa

g3=ggscatter(tmp, x = "TIGIT", y = "CTLA4",
          #add = "reg.line",                                 # Add regression line
          #conf.int = TRUE,                                  # Add confidence interval
          add.params = list(color = "black",fill = "lightgray"),
          color = "Defined_Cell_Subtype_function", palette =  my_color2[7:12],size= 2,
)+facet_grid(. ~ bioGroup)+
geom_point(color='lightgray', shape=21, size = 3, stroke=0.1, alpha= 0.3)+

annotate("text", x=3, y=5, label=paste0('Young: ',aa,'%','\n',
                                       'Intermediated: ',bb,'%','\n',
                                       'Aged: ',cc,'%','\n'),size=4)+
theme(legend.position = "right",
        #legend.justification = c(0,1),
        #legend.title = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        aspect.ratio = 1,
        plot.margin = margin(t=5,r=20,b=5,l=20,unit="pt"),
        axis.title.y = element_text(size=20,colour="black"),
        axis.title.x = element_text(size=20,colour="black"),
        #axis.title.x = element_blank(),
        axis.text = element_text(size=20,colour="black",angle=0),
        axis.text.x = element_text(size=20,colour="black",angle=0,hjust=0.5,vjust=0.5),
        #axis.text.x = element_blank(),
        plot.title = element_text(size=15,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(.1, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.5)
        )
############



pdf(paste0(figdir, '/scater_coexp_ICB_markers_Tcell.pdf'),width = 10,height = 16)
#cowplot::plot_grid(plotlist[[1]],plotlist[[2]],plotlist[[3]],plotlist[[4]],plotlist[[5]],plotlist[[6]],ncol=3)
library(gridExtra)
grid.arrange(g1,g3,ncol=1
    )
dev.off()

save(stat_exp, file= paste0(datadir, '/scater_coexp_ICB_markers_Tcell.csv'))




####END



#1.umaps-clusters, qc, age, doublet score,batch effect
#2.umaps-annotation,seprated by age group
#3.top selective genes
#4.functional gene sets including inhibitory stimulatory genes among clusters
#5.functional gene sets cyto exh genes expression between clusters seprated by age group
#5.proportion of clusters in terms of age,samples(maybe not necessary), correlation with age variable
#6.ROE
7.trajectories of selective clusters, with age ???
8.degs, pathways, metabolism
9.scenic, TFs enriched
#10.cumulative fraction curve along with naive, cytotoxic score, exh score between age group
#11.check the coexpression of naive markers in CD8/CD4 clusters, denploty 
#
#



#1  degs in each cell type between three groups (use heatmap to show, up, down, not sig, three groups for each cell type)
#2  ven diagram of degs among cell types using degs that are differentiated expressed between young and old two groups and show some specific genes (on ven or dot plot)

#4  number of degs (up/down separated) using histgram or Rose diagrams (http://rstudio-pubs-static.s3.amazonaws.com/179757_4fd1fbdcaa3f4f28868bfada3b9855d4.html)   
#5  how degs overlap with well known aging related genes

3  degs during aging using smooth.spline function, heatmap

6  pathways enrichment using aging dependent degs for each cell type (use ggplot dotplot when have many groups)
7 cell cycle phase differences