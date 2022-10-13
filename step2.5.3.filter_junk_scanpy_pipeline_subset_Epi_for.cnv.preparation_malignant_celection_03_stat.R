
####
############################heatmap of priportion of Malignant cell clusters among Samples 
############################
############################
############################
############################
############################
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')
outDir <- file.path("output2_step2.5.3.filter_junk_scanpy_pipeline/data_Epithelial/malignant")
outDir2 <- file.path("output2_step2.5.3.filter_junk_scanpy_pipeline/data_Epithelial/malignant/figure")
#
stat1=read.csv(paste0(outDir, '/data/adata_malignantcell.csv'), header=T, row.names=1)
#stat1$pos=paste(meta.sub$defined.subtype_detailed,meta.sub$subtype_Clusters,sep="_")
stat1$pos=paste0('cluster',stat1$louvain)
stat1$gg=stat1$Sample
#!!!!!stacked malignancy across samples
x=as.data.frame.matrix(table(stat1$pos,stat1$gg))
#colnames(x)=paste0("pt", colnames(x))
# Transform this data in %
#data_percentage <- apply(x, 2, function(x){x*100/sum(x,na.rm=T)})
data_percentage <- apply(x, 2, function(x){x/sum(x,na.rm=T)})
data_percentage = as.data.frame(data_percentage )
###
###
ann_colors = list(
    #age = c("#edf8b1", "#081d58")
	age = rev(viridis::viridis(100))
)

cc= stat1[!duplicated(stat1$Sample), ]
cc$age_num= as.numeric(cc$age)
annotation_col= data.frame(age=cc$age_num )	
rownames(annotation_col)= cc$Sample
head(annotation_col)	

bks <- seq(0,1, by = 0.1)
hmcols <- colorRamps::blue2green2red(length(bks) - 1)

library(RColorBrewer)
library(pheatmap)
hmcols=colorRampPalette(rev(brewer.pal(10, "RdBu")))(256)


hmcols=c(colorRampPalette(c("white", rev(plasma(323, begin = 0.15))[1]))(10), rev(plasma(323, begin = 0.18)))
#hmcols=c(colorRampPalette(c(rev(plasma(323, begin = 0.15))[1]))(10), rev(plasma(323, begin = 0.18)))
#
ph <- pheatmap::pheatmap(data_percentage, 
               #useRaster = T,
               cluster_cols=T, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=T, 
               #clustering_distance_rows=row_dist,
               #clustering_method = 'ward.D2',
               #cutree_rows=num_clusters,
               #silent=TRUE,
               #filename=NA,
               #breaks=bks,
               border_color = 'black',
               color=hmcols,
			   annotation_colors = ann_colors,
			   #annotation_row= annotation_row, 
			   annotation_col=annotation_col
			   )

pdf(paste0(outDir2, '/Epi.Mal.heatmap.proportion_clustersVSSample_v2.pdf'))
ph
dev.off()
############################
############################
############################
############################
############################calculate cell cycle scores
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')
#outDir <- file.path("output2_step2.5.3.filter_junk_scanpy_pipeline/data_Epithelial/malignant")
outDir2 <- file.path("output2_step2.5.3.filter_junk_scanpy_pipeline/data_Epithelial/malignant/figure")
#
library(Seurat)
load('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Epithelial/malignant/sceMalignant_cells.RData')
#
#calculate cell score and cell phase assignment
###
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
#all.sce.epi.mal <- CellCycleScoring(all.sce.epi.mal, s.features = s.genes, g2m.features = g2m.genes,set.ident = TRUE)
#head(all.sce.epi.mal@meta.data)
#assign cell cycle status based on 2 MADs 
#https://github.com/Michorlab/tnbc_scrnaseq/blob/master/code/analysis1.R
#Unravelling subclonal heterogeneity andaggressive disease states in TNBC throughsingle-cell RNA-seq
#
# function to average the expression over a vector of gene indices, for all cells
avg_expr_genes <- function(mat_to_fit, indices){
  scores_cells <- apply(mat_to_fit, 2, function(x){mean(x[indices])})
  return(scores_cells)
}
#
## cell cycle assignment as done in melanoma in Tirosch et al 2016
#melanoma_cellcycle <- read.table(here("data", "melanoma_cellcycle.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#melanoma_g1s <- melanoma_cellcycle$G1S
#melanoma_g1s <- match_clean_vector_genes(mat_norm, melanoma_g1s)
#scores_g1s <- avg_expr_genes(mat_norm, melanoma_g1s$index)
#melanoma_g2m <- melanoma_cellcycle$G2M
#melanoma_g2m <- match_clean_vector_genes(mat_norm, melanoma_g2m)
#scores_g2m <- avg_expr_genes(mat_norm, melanoma_g2m$index)
mat_norm=sce[["RNA"]]@data
scores_g1s <- avg_expr_genes(mat_norm, s.genes[s.genes %in% rownames(sce)])
scores_g2m <- avg_expr_genes(mat_norm, g2m.genes[g2m.genes  %in% rownames(sce)])

cycling_mel <- rep(NA, length(scores_g1s))
for (i in 1:length(cycling_mel)) {
  if (scores_g1s[i] >= (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] < (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "cycling"
  if (scores_g1s[i] < (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] >= (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "cycling"
  if (scores_g1s[i] < (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] < (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "non-cycling"
  if (scores_g1s[i] >= (median(scores_g1s) + 2 * mad(scores_g1s)) && scores_g2m[i] >= (median(scores_g2m) + 2 * mad(scores_g2m)))
    cycling_mel[i] <- "cycling"
}
#
MAD_dericed.cellcycle=data.frame(scores_g1s=scores_g1s, scores_g2m=scores_g2m, cycling_mel=cycling_mel)
head(MAD_dericed.cellcycle)
#######
sce=AddMetaData(sce,metadata = MAD_dericed.cellcycle)

head(sce@meta.data)
#
meta=sce@meta.data
meta$age_num=as.numeric(paste(meta$age))
meta$age_group=cut(
  meta$age_num,
  breaks = c(20, 49, 60,100),
  labels = c("Young(20-40)", "intermediated(41-60)", "Aged(61-80)")
)
sce <- AddMetaData(
  object = sce,
  metadata = meta[,c('age_group', 'age_num')]
)
#
#
color_CLASS = c(
                 '#0000FF', # 0
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
library(ggplot2)
pdf(paste0(outDir2, '/proportion.ageGroup_vs_cellcycling.pdf'),width=3,height = 5)
ggplot(sce@meta.data, aes(x=age_group, fill=cycling_mel)) + 
geom_bar(position = "fill")+scale_fill_manual(values = color_CLASS)+
 theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())+
#coord_flip()+
ylab('percentage')
dev.off()

############################################################################################
############################################################################################
#2021-06-08
############################calculate biodiversity score
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')
#outDir <- file.path("output2_step2.5.3.filter_junk_scanpy_pipeline/data_Epithelial/malignant")
outputdir <- file.path("output2_step2.5.3.filter_junk_scanpy_pipeline/data_Epithelial/malignant/biodiversity")
if(!dir.exists(outputdir)) dir.create(outputdir,recursive=TRUE)
#

library(Seurat)
load('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Epithelial/malignant/sceMalignant_cells.RData')
#
sce <- RunPCA(object = sce, features = VariableFeatures(object = sce),npcs=100)

embeds1 = Embeddings(sce[["pca"]])
tmpmeta1=cbind(embeds1, sce@meta.data)
tmpmeta1$ptid= tmpmeta1$Sample
tmpmeta1$age_num=as.numeric(paste(tmpmeta1$age))
tmpmeta1$age_group=cut(
  tmpmeta1$age_num,
  breaks = c(20, 49, 60,100),
  labels = c("Young(20-40)", "intermediated(41-60)", "Aged(61-80)")
)
#

#########
##########calculate biodiversity score
biodata=tmpmeta1
id=names(table(biodata$ptid))
#1.detect extreme value for each tumor
#based on the first three PCs
datalist=list()

for (k in 1:length(id)) {
  sub.biodata=subset(biodata, ptid==id[k])
  
  rangePC1.low=colMeans(sub.biodata[,1:3])[1]-3*apply(sub.biodata[,1:3],2,sd)[1]
  rangePC1.up=colMeans(sub.biodata[,1:3])[1]+3*apply(sub.biodata[,1:3],2,sd)[1]
  #
  rangePC2.low=colMeans(sub.biodata[,1:3])[2]-3*apply(sub.biodata[,1:3],2,sd)[2]
  rangePC2.up=colMeans(sub.biodata[,1:3])[2]+3*apply(sub.biodata[,1:3],2,sd)[2]
  #
  rangePC3.low=colMeans(sub.biodata[,1:3])[3]-3*apply(sub.biodata[,1:3],2,sd)[3]
  rangePC3.up=colMeans(sub.biodata[,1:3])[3]+3*apply(sub.biodata[,1:3],2,sd)[3]
  #judge if a cell with extreme value
  dim(sub.biodata)
  subdata_sel=subset(sub.biodata, (PC_1 > rangePC1.low & PC_1 < rangePC1.up ) |  
                       (PC_2 > rangePC2.low & PC_2 < rangePC2.up ) | 
                       (PC_3 > rangePC3.low & PC_3 < rangePC3.up )
  )
  datalist[[k]]=subdata_sel
}
biodata_filtered=do.call(rbind, datalist)
#2.calculate score by chosing different PCs numbers
head(biodata_filtered)
biopcs=biodata_filtered[,1:60]
PCsmean=colMeans(biopcs)
biopcs2=sweep(as.matrix(biopcs), 2, PCsmean, FUN = "-")
biopcs3=biopcs2^2
#####
bioscore=data.frame(ptid=biodata_filtered$ptid,age=as.numeric(paste(biodata_filtered$age)),age_group= biodata_filtered$age_group,stage= biodata_filtered$Stages,
biodiversity.20s=rowMeans(biopcs3[,1:20])
,biodiversity.25s=rowMeans(biopcs3[,1:25])
,biodiversity.30s=rowMeans(biopcs3[,1:30])
,biodiversity.35s=rowMeans(biopcs3[,1:35])
,biodiversity.40s=rowMeans(biopcs3[,1:40])
,biodiversity.45s=rowMeans(biopcs3[,1:45])
,biodiversity.50s=rowMeans(biopcs3[,1:50])
,biodiversity.55s=rowMeans(biopcs3[,1:55])
,biodiversity.60s=rowMeans(biopcs3[,1:60])
)
head(bioscore)
#save(bioscore, file="data/scRNA/scPanCompare/6.seurat.epiMal.filtered/noMTadjusted/bioscore.blueprint.RData")
####
####
p<-ggplot(bioscore, aes(x=age_group, y=biodiversity.30s, fill=age_group)) +
  geom_boxplot()

pdf(paste0(outputdir, '/boxplot.bioscore.age_group.pdf'))
p
dev.off()
##############
library(ggpubr)
pdf(paste0(outputdir, '/Cor_age_vs_bioscore.pdf'),width=8,height = 8)
ggscatter(bioscore, x = "age", y = "biodiversity.30s", size = 0.5,
          add = "reg.line", conf.int = F, add.params = list(color = "#253494"),
          cor.coef = T, cor.method = "pearson",
          xlab = "age", ylab = "biodiversity")+
  ggtitle("age_vs_cellPercentage_pearson")+#geom_point(aes(color="#dd1c77"))
  geom_point(fill="#dd1c77",color="#dd1c77")#+facet_wrap("Defined_Cell_Subtype_function")
dev.off()

#####
#####
library(dplyr)
library(tidyr)
bioscore.patient= bioscore[, -c(2,3,4)] %>% group_by(ptid) %>% summarise_all(mean)
bioscore.patient=as.data.frame(bioscore.patient)
rownames(bioscore.patient)=bioscore.patient$ptid
head(bioscore.patient)
####add age information
ageinfo= bioscore[!duplicated(bioscore$ptid), c(1:4)]
head(ageinfo)
####
bioscore.patient=merge(ageinfo, bioscore.patient, by='ptid', all=F)
rownames(bioscore.patient)=bioscore.patient$ptid
bioscore.patient
#####


####
bioscore.patient.draw=bioscore.patient %>% pivot_longer(cols=colnames(bioscore.patient)[-c(1,2,3,4)],
                                                                          names_to= "No.PCs",
                                                                          values_to = "Bdiversity.score")
bioscore.patient.draw=as.data.frame(bioscore.patient.draw)
bioscore.patient.draw
#bioscore.patient.draw$ptid=factor(bioscore.patient.draw$ptid,levels = c("pt33","pt38","pt31","pt32","pt37","pt35"))
#########rescale
#
library(scales)
bioscore.patient=bioscore.patient[,-c(1:4)]
bioscore.patient.rescaled=apply(bioscore.patient,2,rescale, to = c(0,1))
bioscore.patient.rescaled=as.data.frame(bioscore.patient.rescaled)
bioscore.patient.rescaled$ptid=rownames(bioscore.patient.rescaled)
bioscore.patient.rescaled$age=ageinfo$age
bioscore.patient.rescaled$age_group=ageinfo$age_group
bioscore.patient.rescaled$stage=ageinfo$stage
head(bioscore.patient.rescaled)

library(tidyr)
bioscore.patient.rescaled.draw=bioscore.patient.rescaled %>% pivot_longer(cols=colnames(bioscore.patient.rescaled)[1:9],
                                                                              names_to= "No.PCs",
                                                                              values_to = "Scaled.diversity.score")
bioscore.patient.rescaled.draw=as.data.frame(bioscore.patient.rescaled.draw)
#bioscore.patient.rescaled.draw$ptid=factor(bioscore.patient.rescaled.draw$ptid,levels = c("pt33","pt38","pt31","pt32","pt37","pt35"))
#
#save(bioscore.patient.draw,bioscore.patient.rescaled.draw,
 #    file="data/scRNA/scPanCompare/6.seurat.epiMal.filtered/noMTadjusted/bioscore.draw.across.pt.blueprint.RData")
######
library(ggplot2)
library(gridExtra)
pdf(paste0(outputdir, "/biscore.vs.PCs.combination.pdf"),width = 12,height = 5)
grid.arrange(
ggplot(data=bioscore.patient.draw, aes(x=age, y=Bdiversity.score, group=No.PCs, color=No.PCs)) +
  geom_line() + geom_point()+
  #scale_color_brewer(palette="Paired")+
  theme_minimal(),
#
ggplot(data=bioscore.patient.rescaled.draw, aes(x=age, y=Scaled.diversity.score, group=No.PCs, color=No.PCs)) +
  geom_line() + geom_point()+
  #scale_color_brewer(palette="Paired")+
  theme_minimal()
,ncol=2
)
dev.off()
#####
bioscore.patient.draw.group=subset(bioscore.patient.draw, No.PCs=="biodiversity.30s")
bioscore.patient.draw.group

bioscore.patient.rescaled.draw_p30= subset(bioscore.patient.rescaled.draw,  No.PCs== 'biodiversity.30s')
bioscore.patient.rescaled.draw_p30
#
pdf(paste0(outputdir, "/bioscore.median.group.pdf"),width = 5,height = 6)
ggboxplot(bioscore.patient.draw.group, x = "age_group", y = "Bdiversity.score",
          color = "age_group", #palette = c("Div-High"="#d01c8b","Div-Low"="#4dac26"),
          add = "jitter")+stat_compare_means(aes(group = age_group))

dev.off()

library(ggpubr)
######
pdf(paste0(outputdir, "/bioscore.median.stage_group.pdf"),width = 5,height = 6)
ggboxplot(bioscore.patient.draw.group, x = "stage", y = "Bdiversity.score",
          color = "stage", #palette = c("Div-High"="#d01c8b","Div-Low"="#4dac26"),
          add = "jitter")+stat_compare_means(aes(group = age_group))
ggboxplot(bioscore.patient.rescaled.draw_p30, x = "stage", y = "Scaled.diversity.score",
          color = "stage", #palette = c("Div-High"="#d01c8b","Div-Low"="#4dac26"),
          add = "jitter")+stat_compare_means(aes(group = age_group))

dev.off()
#
#
save(bioscore.patient.draw.group,bioscore.patient.rescaled.draw_p30 , file= paste0(outputdir, "/bioscore.patient_level.RData" ))