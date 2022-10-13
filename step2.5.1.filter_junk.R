setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')

outDir <- file.path("output2_step2.5.filter_junk")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

source('/data/Zhuxq/young_LC_analysis/my_colors.R')
#Load required packages

library(Seurat)
library(dplyr)
library(magrittr)
library(data.table)
library(Matrix)
library(devtools)
library(RcppArmadillo)
library(Rcpp)
library(scales)
library(pheatmap)
library(gplots)
library(ggplot2)
library(cowplot)
library(tibble)
library(harmony)
#
#
load('output2_step2.1.qc_scaledata/TotalTissue.combined.tmp.RData')
load('output2_step2.2.unimodel.prediction/my.combined.unimodel.celltype.prediced_metadata.RData')
load('output2_step2.3.harmony.overcluster/TotalTissue.combined.harmony_overclustered.RData')
scrublet_1= read.csv('output2_step2.4.doublet.scrublet.python/mergedadatalistScrubletOutput_default.csv', header=T, row.names=1)
scrublet_2= read.csv('output2_step2.4.doublet.scrublet.python/mergedadatalistScrubletOutput_neighbors30.csv', header=T, row.names=1)
#add data into metadata
#1
aa= meta.merge.unimodel
aa=aa[rownames(TotalTissue.combined@meta.data),]
TotalTissue.combined<- AddMetaData(object = TotalTissue.combined, metadata = aa$Cell_subtype, col.name = 'Cell_subtype.unimodel')
TotalTissue.combined<- AddMetaData(object = TotalTissue.combined, metadata = aa$Cell_type.refined, col.name = 'Cell_type.unimodel')
#2
aa=TotalTissue.combined.harmony_overclustered@meta.data
aa=aa[rownames(TotalTissue.combined@meta.data),]
TotalTissue.combined<- AddMetaData(object = TotalTissue.combined, metadata = aa$UMAP_Clusters.harmony, col.name = 'UMAP_Clusters.harmony_overclustered')
#3
aa=scrublet_1
aa=aa[rownames(TotalTissue.combined@meta.data),]
TotalTissue.combined<- AddMetaData(object = TotalTissue.combined, metadata = aa$scrubletScore, col.name = 'scrubletScore_default')
TotalTissue.combined<- AddMetaData(object = TotalTissue.combined, metadata = aa$scrublet_predict, col.name = 'scrublet_predict_default')
#4
aa=scrublet_2
aa=aa[rownames(TotalTissue.combined@meta.data),]
TotalTissue.combined<- AddMetaData(object = TotalTissue.combined, metadata = aa$scrubletScore, col.name = 'scrubletScore_n30')
TotalTissue.combined<- AddMetaData(object = TotalTissue.combined, metadata = aa$scrublet_predict, col.name = 'scrublet_predict_n30')
#
table(TotalTissue.combined$scrublet_predict_default,  TotalTissue.combined$scrublet_predict_n30)

save(TotalTissue.combined, file=paste0(outDir, '/TotalTissue.combined.not.doubletFiltered.RData'))

#use default one
#filter
#scesel=subset(TotalTissue.combined,  scrubletScore_default <= 0.7)
#scesel
#
aa=TotalTissue.combined@meta.data[, c('scrubletScore_default', 'UMAP_Clusters.harmony_overclustered','nFeature_RNA')]
head(aa)
#calculate mean scrublet score in overclusterd clusters
library(dplyr)
clusterscore=aa %>% group_by(UMAP_Clusters.harmony_overclustered) %>% summarise_all(mean)
clusterscore=clusterscore[order(clusterscore$scrubletScore_default, decreasing =T),]
clusterscore$UMAP_Clusters.harmony_overclustered=factor(clusterscore$UMAP_Clusters.harmony_overclustered, levels=clusterscore$UMAP_Clusters.harmony_overclustered)

p1<-ggplot(data=as.data.frame(clusterscore), aes(x=UMAP_Clusters.harmony_overclustered, y=scrubletScore_default)) + 
geom_bar(stat="identity")+coord_flip()

clusterscore=clusterscore[order(clusterscore$nFeature_RNA, decreasing =F),]
clusterscore$UMAP_Clusters.harmony_overclustered=factor(clusterscore$UMAP_Clusters.harmony_overclustered, levels=clusterscore$UMAP_Clusters.harmony_overclustered)

p2<-ggplot(data=as.data.frame(clusterscore), aes(x=UMAP_Clusters.harmony_overclustered, y=nFeature_RNA)) + 
geom_bar(stat="identity")+coord_flip()
library(gridExtra)
pdf(paste0(outDir, "/Doubletscore.barplot.pdf"), width=8, height=40)
grid.arrange(p1,p2, ncol=2)
dev.off()

#
#no clusters with mean score more than 0.6
#so no filtering further for potentional boublets 
#######
p3=ggplot(TotalTissue.combined@meta.data, aes(x=UMAP_Clusters.harmony_overclustered, fill=scrublet_predict_default)) + 
geom_bar(position = "fill")+coord_flip()+ggtitle('scrublet')+theme_classic()+geom_hline(yintercept=0.05)
p4=ggplot(TotalTissue.combined@meta.data, aes(x=UMAP_Clusters.harmony_overclustered, fill=scrublet_predict_default)) + 
geom_bar()+coord_flip()+ggtitle('scrublet')+theme_classic()+geom_hline(yintercept=100)

#p4=ggplot(scesel@meta.data, aes(x=UMAP_Clusters.harmony_overclustered, fill=scrublet_predict_default)) + geom_bar(position = "fill")+coord_flip()+ggtitle('scrublet')+theme_classic()
#scale_fill_manual(values = FLATUI_CLASS[3:2])+
#p2=ggplot(drawdata)+aes(x = RNA_snn_res.0.5, y = nFeature_RNA, fill = RNA_snn_res.0.5) + # add color to boxes with fillgeom_boxplot(varwidth = F)
library(gridExtra)
pdf(paste0(outDir, "/Doublet.proportion.pdf"), width=12, height=40)
grid.arrange(p3,p4, ncol=2)
dev.off()
#



Epithelialfeatures = c('EPCAM', 'KRT19', 'CDH1', 'KRT18')
TCELLSfeatures = c('CD3D', 'CD3E', 'CD3G', 'TRAC')
NKfeatures = c('NCAM1', 'NKG7', 'GNLY', 'KLRD1', 'KLRF1')
BCELLSfeatures = c('CD79A', 'IGHM', 'IGHG3', 'IGHA2')
MYELOIDfeatures = c('CD68', 'MARCO', 'FCGR3A', 'LYZ','AIF1')
MASTCELLSfeatures = c('KIT', 'MS4A2','GATA2','TPSB2')
Fibroblastsfeatures = c('DCN', 'COL1A1', 'COL1A2', 'THY1' )
ENDOTHELIALfeatures = c('PECAM1', 'CLDN5', 'FLT1', 'RAMP2')
Cyclingfeatures=c('STMN1', 'MKI67', 'TOP2A', 'CDK1')

features <- list("Epithelial" = Epithelialfeatures,  "T cell" = TCELLSfeatures, "NK" = NKfeatures, "B cell" = BCELLSfeatures, "Myeloid" = MYELOIDfeatures,
                  'MAST cell'= MASTCELLSfeatures, 'Fibroblast'= Fibroblastsfeatures, 'Endothelial'= ENDOTHELIALfeatures, 'Cycling'= Cyclingfeatures)
#
pdf(paste0(outDir, '/dotplot.plots.cellmarkers.main.lineage.pdf'),width = 14,height = 80)
DotPlot(object = TotalTissue.combined.harmony_overclustered, features=features, cluster.idents=T) + theme(axis.text.x = element_text(angle = 90))
dev.off()


#
#
#
junkcluster=c('61','83','102','133','255','285','354','358','379','400','403','404',
'405','423','428','441','448','457','460','462','465','466','467','469','472')

dcls1 <- WhichCells(TotalTissue.combined.harmony_overclustered, idents = c(junkcluster))

p5=DimPlot(TotalTissue.combined.harmony_overclustered, label=F, group.by="UMAP_Clusters.harmony",
cells.highlight= list(dcls1), 
cols.highlight = FLATUI_CLASS, cols= "grey")
pdf(paste0(outDir, "/Junk.clusters.highlight.pdf"))
p5
dev.off()
#
#
#now subset the cluster
Idents(TotalTissue.combined)= 'UMAP_Clusters.harmony_overclustered'
levels(TotalTissue.combined)
#
#
TotalTissue.combined_fil=subset(TotalTissue.combined, idents=junkcluster, invert = TRUE)
TotalTissue.combined_fil
save(TotalTissue.combined_fil, file=paste0(outDir, '/TotalTissue.combined_fil.RData'))
#

###

