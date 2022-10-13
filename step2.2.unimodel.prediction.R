
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')

outDir <- file.path("output2_step2.1.qc_scaledata")

outDir2 <- file.path("output2_step2.2.unimodel.prediction")
if(!dir.exists(outDir2)) dir.create(outDir2,recursive=TRUE)
#Load required packages

library(Seurat)
#
#
# 2. Data Pre-processing
#Merge all objects
load(paste0(outDir, "/TotalTissue.combined.tmp.RData"))



pancreas.list <- SplitObject(TotalTissue.combined, split.by = "Sample")

refidx=names(pancreas.list)[grep('LUNG', names(pancreas.list))]

reference.list <- pancreas.list[refidx]
pancreas.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

pancreas.integrated <- IntegrateData(anchorset = pancreas.anchors, dims = 1:30)

library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically set during
# IntegrateData
DefaultAssay(pancreas.integrated) <- "integrated"
# Run the standard workflow for visualization and clustering
pancreas.integrated <- ScaleData(pancreas.integrated, verbose = FALSE)
pancreas.integrated <- RunPCA(pancreas.integrated, npcs = 30, verbose = FALSE)
pancreas.integrated <- RunUMAP(pancreas.integrated, reduction = "pca", dims = 1:30, verbose = FALSE)

#

pancreas.integrated <- RunUMAP(pancreas.integrated, dims = 1:30, reduction = "pca", return.model = TRUE)


Idents(TotalTissue.combined)= 'orig.ident'

pancreas.query=subset(TotalTissue.combined,  idents= 'mysce')


pancreas.anchors <- FindTransferAnchors(reference = pancreas.integrated, query = pancreas.query, 
    dims = 1:30)

pancreas.query <- MapQuery(anchorset = pancreas.anchors, reference = pancreas.integrated, query = pancreas.query, 
    refdata = list(Cell_subtype = "Cell_subtype", Cell_type.refined= 'Cell_type.refined' ), reference.reduction = "pca", reduction.model = "umap")

#
p1 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "Cell_type.refined", label = F, label.size = 3, 
    repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "Cell_subtype", label = F) + ggtitle("Reference annotations")	

p2 <- DimPlot(pancreas.integrated, reduction = "umap", group.by = "Cell_subtype", label = F) + ggtitle("Reference annotations")	

p3 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.Cell_type.refined", label = F) + ggtitle("Query transferred labels")
p4 <- DimPlot(pancreas.query, reduction = "ref.umap", group.by = "predicted.Cell_subtype", label = F) + ggtitle("Query transferred labels")


pdf(paste0(outDir2, '/umaps.integrated.ref.Lung.pdf'), width=14, height=6)
p1 + p2 
p3 +p4
dev.off()
#####
#####
#####
ameta=pancreas.integrated@meta.data
bmeta=pancreas.query@meta.data
bmeta$Cell_subtype= bmeta$predicted.Cell_subtype
bmeta$Cell_type.refined= bmeta$predicted.Cell_type.refined
meta.merge=rbind(ameta[,c('Cell_subtype','Cell_type.refined')],
           bmeta[,c('Cell_subtype','Cell_type.refined')]
)
####
TotalTissue.combined<- AddMetaData(object = TotalTissue.combined, 
	metadata = meta.merge$Cell_subtype, col.name = 'Cell_subtype.merge')
TotalTissue.combined<- AddMetaData(object = TotalTissue.combined, 
		metadata = meta.merge$Cell_type.refined, col.name = 'Cell_type.merge')	
	
saveRDS(TotalTissue.combined, file=paste0(outDir2, '/my.combined.unimodel.celltype.prediced.rds'))

save(pancreas.query, pancreas.integrated, file=paste0(outDir2, '/my.combined.unimodel.results.celltype.prediced.RData'))

meta.merge.unimodel=meta.merge
save(meta.merge.unimodel, file=paste0(outDir2, '/my.combined.unimodel.celltype.prediced_metadata.RData'))