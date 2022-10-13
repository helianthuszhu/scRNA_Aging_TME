setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')

outDir <- file.path("output2_step2.3.harmony.overcluster")
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
load(paste0('output2_step2.1.qc_scaledata', "/TotalTissue.combined.tmp.RData"))

#
#Data Visualization

#Run PCA and Determine Dimensions for 90% Variance
TotalTissue.combined <- RunPCA(object = TotalTissue.combined, features = VariableFeatures(object = TotalTissue.combined))
#
#
stdev <- TotalTissue.combined@reductions$pca@stdev
var <- stdev^2

EndVar = 0

for(i in 1:length(var)){
  total <- sum(var)
  numerator <- sum(var[1:i])
  expvar <- numerator/total
  if(EndVar == 0){
    if(expvar > 0.9){
      EndVar <- EndVar + 1
      PCNum <- i
    }
  }
}
#Confirm #PC's determined explain > 90% of variance
sum(var[1:PCNum])/ sum(var)
#
#run harmony
#
#
TotalTissue.combined.harmony <- TotalTissue.combined %>% RunHarmony("Sample", plot_convergence = F,max.iter.harmony = 20)

TotalTissue.combined.harmony <- TotalTissue.combined.harmony %>% 
    RunUMAP(reduction = "harmony", dims = 1:PCNum) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:PCNum) %>% 
	FindClusters(resolution = 50) %>%
    identity()


TotalTissue.combined.harmony[["UMAP_Clusters.harmony"]] <- Idents(object = TotalTissue.combined.harmony)
#
#
TotalTissue.combined.harmony_overclustered=TotalTissue.combined.harmony
save(TotalTissue.combined.harmony_overclustered, file = paste0(outDir, '/TotalTissue.combined.harmony_overclustered.RData') )

levels(TotalTissue.combined.harmony)
####
####
pdf(paste0(outDir, '/umaps.harmony.overclusterd.pdf'), width=50)
#library(gridExtra)
DimPlot(TotalTissue.combined.harmony, group.by = c("UMAP_Clusters.harmony"), label=T)
dev.off()
#
#
#
#