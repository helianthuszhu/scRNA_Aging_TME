setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')

outDir <- file.path("output2_step2.1.qc_scaledata")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)
	
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
# 2. Data Pre-processing
#Merge all objects
TotalTissue.combined =readRDS('total.tissue.combined.addMTRPS.percentage.rds')

# Changing between meta.data for identities- you can change this by altering what you input into your metadata
Idents(object = TotalTissue.combined) <- 'Sample'
# Check active identity
levels(TotalTissue.combined)
#QC
TotalTissue.combined[["percent.mt.merge"]] <- PercentageFeatureSet(TotalTissue.combined, pattern = "^MT-")
TotalTissue.combined[["percent.rb.merge"]] <- PercentageFeatureSet(TotalTissue.combined, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")

# Visualize QC metrics as a violin plot
p1=VlnPlot(TotalTissue.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rb', "percent.mt.merge",'percent.rb.merge'), ncol = 4)
p2 <- FeatureScatter(TotalTissue.combined, feature1 = "nCount_RNA", feature2 = "percent.mt.merge")
p3 <- FeatureScatter(TotalTissue.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#####
pdf(paste0(outDir, '/QC.plot.pdf'))
p1;p2;p3
dev.off()
#filter cells
#
TotalTissue.combined <- subset(TotalTissue.combined, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt.merge < 20 & nCount_RNA > 2000 & nCount_RNA < 40000)

saveRDS(TotalTissue.combined, file=paste0(outDir, "/TotalTissue.combined.qc.filtered.rds"))

p1=VlnPlot(TotalTissue.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt.merge",'percent.rb.merge'), ncol = 4)
p2 <- FeatureScatter(TotalTissue.combined, feature1 = "nCount_RNA", feature2 = "percent.mt.merge")
p3 <- FeatureScatter(TotalTissue.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#####
pdf(paste0(outDir, '/QC.plot2.pdf'))
p1;p2;p3
dev.off()

#
#
#
#NORMALIZE DATA
TotalTissue.combined <- NormalizeData(object = TotalTissue.combined, normalization.method = "LogNormalize", 
                                      scale.factor = 10000)

#FIND VARIABLE GENES
TotalTissue.combined<- FindVariableFeatures(TotalTissue.combined, selection.method = "vst", nfeatures = 2000)


#Calculate Cell Cycle Score (S-G2M Difference)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
TotalTissue.combined<- CellCycleScoring(TotalTissue.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
TotalTissue.combined$CC.Difference <- TotalTissue.combined$S.Score - TotalTissue.combined$G2M.Score

#THIS STEP MAY TAKE A VERY LONG TIME
#Scale Data
TotalTissue.combined<- ScaleData(object = TotalTissue.combined, vars.to.regress = c("nCount_RNA",'CC.Difference', 'percent.mt.merge'), 
                           features = rownames(TotalTissue.combined))
#
save(TotalTissue.combined,file=paste0(outDir, "/TotalTissue.combined.tmp.RData"))
#
#
#
library(DropletUtils)
idx= unique(TotalTissue.combined$Sample)
da= TotalTissue.combined
Idents(da)= 'Sample'
for (1 in 1:length(idx)) {
	aa=subset(da, idents=idx[1])
	
	outDir <- file.path("output2_step2.1.qc_scaledata")
	
}
output.path=paste0(outDir, '/10xFormat')
write10xCounts(x = TotalTissue.combined@assays$RNA@counts, path = output.path)
write.table(TotalTissue.combined@meta.data, file=paste0(output.path, '/meta_data.tsv'), sep='\t', quote=F)




###
###
###2022-10-12
###*********draw qc plot for revirew
###*********************************
###*********************************
###*********************************
###*********************************
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')

outDir <- file.path("output2_step2.1.qc_scaledata")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)
	
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
# 2. Data Pre-processing
#Merge all objects
TotalTissue.combined =readRDS('total.tissue.combined.addMTRPS.percentage.rds')
#load(paste0(outDir, '/TotalTissue.combined.not.doubletFiltered.RData'))

### set qc para
nFeature_lower <- 500
nFeature_upper <- 5000
nCount_lower <- 2000
nCount_upper <- 40000
pMT_lower <- 0
pMT_upper <- 20
pHB_lower <- 0
pHB_upper <- 5
pRP_lower <- 3
pRP_upper <- 60



seu_obj=TotalTissue.combined

### calculate mitochondrial, hemoglobin and ribosomal gene counts
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^MT-", col.name = "pMT")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^HBA|^HBB", col.name = "pHB")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^RPS|^RPL", col.name = "pRP")


qcparams <- c("nFeature_RNA", "nCount_RNA", "pMT"#, "pHB", "pRP"
             )
gp1= VlnPlot(seu_obj, features ='nFeature_RNA', pt.size = 0, group.by = "orig.ident", ncol = 1,log =F)+
      geom_hline(yintercept = nFeature_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = nFeature_upper, color = "red", linetype = 2) 

gp2= VlnPlot(seu_obj, features ='nCount_RNA', pt.size = 0, group.by = "orig.ident", ncol = 1,log =F)+
      geom_hline(yintercept = nCount_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = nCount_upper, color = "red", linetype = 2) 

gp3= VlnPlot(seu_obj, features ='pMT', pt.size = 0, group.by = "orig.ident", ncol = 1,log =F)+
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) 

#gp2= RidgePlot(object = seu_obj, features = qcparams, group.by = "orig.ident", ncol = 1,log=T)

pdf(paste0(outDir,"/qc_raw_binder_study.pdf" ), width = 4, 
        height = 9)
cowplot::plot_grid(gp1,gp2,gp3,ncol=1)
dev.off()




qcparams <- c("nFeature_RNA", "nCount_RNA", "pMT")
seu_obj$Sample= factor(seu_obj$Sample, levels=c('LUNG_T09','LUNG_T34','LUNG_T18','LUNG_T25','LUNG_T31','LUNG_T28','LUNG_T19','LUNG_T30','LUNG_T06','LUNG_T08','LUNG_T20',
                                                'Y1','Y2','Y3','Y4','E1','E2','E3','E4'))
#sample each
gp1= VlnPlot(seu_obj, features =qcparams[1], pt.size = 0, group.by = "Sample", ncol = 1,log=F)+
geom_hline(yintercept = nFeature_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = nFeature_upper, color = "red", linetype = 2) 
gp2= VlnPlot(seu_obj, features =qcparams[2], pt.size = 0, group.by = "Sample", ncol = 1,log=F)+
geom_hline(yintercept = nCount_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = nCount_upper, color = "red", linetype = 2) 

gp3= VlnPlot(seu_obj, features =qcparams[3], pt.size = 0, group.by = "Sample", ncol = 1,log=F)+
geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2)
#gp2= RidgePlot(object = seu_obj, features = qcparams, group.by = "Sample", ncol = 1,log=T)

pdf(paste0(outDir,"/qc_raw_binder_sample.pdf" ), width = 10, 
        height = 10)
cowplot::plot_grid(gp1,gp2,gp3,ncol=1)
dev.off()






# Data Filtering 

qc_std_plot_helper <- function(x) x + 
  scale_color_viridis() +
  geom_point(size = 0.01, alpha = 0.3)

qc_std_plot <- function(seu_obj) {
  qc_data <- as_tibble(FetchData(seu_obj, c("nCount_RNA", "nFeature_RNA", "pMT", "pHB", "pRP")))
  plot_grid(
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pMT))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pHB))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pRP))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pMT, color = nFeature_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    #qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pHB, color = nFeature_RNA))) + 
     # geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
     # geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
     # geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
     # geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    #qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pRP, color = nFeature_RNA))) + 
     # geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
     # geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2)+
     # geom_hline(yintercept = pRP_lower, color = "red", linetype = 2) +
     # geom_hline(yintercept = pRP_upper, color = "red", linetype = 2),
    
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pMT, color = nCount_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    #qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pHB, color = nCount_RNA))) + 
    #  geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
    #  geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
    #  geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
    #  geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    #qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pRP, color = nCount_RNA))) + 
    #  geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
    #  geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2)+
    #  geom_hline(yintercept = pRP_lower, color = "red", linetype = 2) +
    #  geom_hline(yintercept = pRP_upper, color = "red", linetype = 2),
    
    #qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nCount_RNA))) + 
    #  geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
    #  geom_hline(yintercept = pMT_upper, color = "red", linetype = 2)+
    #  geom_vline(xintercept = pRP_lower, color = "red", linetype = 2) +
    #  geom_vline(xintercept = pRP_upper, color = "red", linetype = 2),
    #qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nFeature_RNA))) + 
    #  geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
    #  geom_hline(yintercept = pMT_upper, color = "red", linetype = 2)+
    #  geom_vline(xintercept = pRP_lower, color = "red", linetype = 2) +
    #  geom_vline(xintercept = pRP_upper, color = "red", linetype = 2),
    
    
    #ggplot(gather(qc_data, key, value), aes(key, value)) +
    #  geom_violin() +
    #  facet_wrap(~key, scales = "free", ncol = 5),
    
    ncol = 5, align = "hv"
  )
}
#

seu_obj_unfiltered <- seu_obj

pdf(paste0(outDir,'/qc_raw_binder_cor.pdf'),width=40,height=30)
qc_std_plot(seu_obj_unfiltered)
dev.off()

png(paste0(outDir,'/qc_raw_binder_cor.png'),width = 40, height = 30, units = "cm",res=150)
qc_std_plot(seu_obj_unfiltered)
dev.off()

png(paste0(outDir,'/qc_raw_binder_cor2.png'),width = 60, height = 10, units = "cm",res=150)
qc_std_plot(seu_obj_unfiltered)
dev.off()

pdf(paste0(outDir,'/qc_raw_binder_cor2.pdf'),width=60,height=10)
qc_std_plot(seu_obj_unfiltered)
dev.off()


qc_data <- as_tibble(FetchData(seu_obj, c("nCount_RNA", "nFeature_RNA", "pMT", "pHB", "pRP")))


g1=qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pMT))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2)

 g2=   qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pHB))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2)

  g3=  qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pRP))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2)
    
  g4=  qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pMT, color = nFeature_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2)

g5= qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pMT, color = nCount_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2)


pdf(paste0(outDir,"/qc_raw_binder_cor_f1.pdf" ), width = 6, 
        height = 5)
cowplot::plot_grid(g1, ncol=1)
dev.off()

pdf(paste0(outDir,"/qc_raw_binder_cor_f2.pdf" ), width = 6, 
        height = 5)
cowplot::plot_grid(g2, ncol=1)
dev.off()

pdf(paste0(outDir,"/qc_raw_binder_cor_f3.pdf" ), width = 6, 
        height = 5)
cowplot::plot_grid(g3, ncol=1)
dev.off()

pdf(paste0(outDir,"/qc_raw_binder_cor_f4.pdf" ), width = 6, 
        height = 5)
cowplot::plot_grid(g4, ncol=1)
dev.off()

pdf(paste0(outDir,"/qc_raw_binder_cor_f5.pdf" ), width = 6, 
        height = 5)
cowplot::plot_grid(g5, ncol=1)
dev.off()


png(paste0(outDir,"/qc_raw_binder_cor_f1.png" ), width = 10, height = 10,units = "cm",res=150)
cowplot::plot_grid(g1, ncol=1)
dev.off()

png(paste0(outDir,"/qc_raw_binder_cor_f2.png" ), width = 10, height = 10,units = "cm",res=150)
cowplot::plot_grid(g2, ncol=1)
dev.off()

png(paste0(outDir,"/qc_raw_binder_cor_f3.png" ), width = 10, height = 10,units = "cm",res=150)
cowplot::plot_grid(g3, ncol=1)
dev.off()

png(paste0(outDir,"/qc_raw_binder_cor_f4.png" ), width = 10, height = 10,units = "cm",res=150)
cowplot::plot_grid(g4, ncol=1)
dev.off()

png(paste0(outDir,"/qc_raw_binder_cor_f5.png" ), width = 10, height = 10,units = "cm",res=150)
cowplot::plot_grid(g5, ncol=1)
dev.off()