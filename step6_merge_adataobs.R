#load the seruat object
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')

library("Seurat")

load('output2_step2.5.filter_junk/TotalTissue.combined.not.doubletFiltered.RData')
scetotal= TotalTissue.combined
#
obsdata= read.csv('output6_merge_adataobs/obs_merged.csv', header=T,row.names=1)
head(obsdata)
#filter
obsdata_sel= obsdata[!grepl('Epithelial', obsdata$Defined_Cell_Subtype_function),]
obsdata_sel= obsdata_sel[!grepl('Doublet', obsdata_sel$Defined_Cell_Subtype_function),]
obsdata_sel= obsdata_sel[!grepl('FewMarker', obsdata_sel$Defined_Cell_Subtype_function),]
#
obsdata_sel$Defined_Cell_Subtype_function = gsub('/', '-', obsdata_sel$Defined_Cell_Subtype_function)
obsdata_sel$Defined_Cell_Subtype_function = gsub(' ', '-', obsdata_sel$Defined_Cell_Subtype_function)

#
head(scetotal@meta.data)
#
sce=subset(scetotal[, rownames(obsdata_sel)])
head(sce@meta.data)

sce <- AddMetaData(
  object = sce,
  metadata = obsdata_sel$Defined_Cell_Subtype_function,
  col.name = 'Defined_Cell_Subtype_function'
)
head(sce@meta.data)
####################
####################findmarker
Idents(sce)='Defined_Cell_Subtype_function'
out_marker <- FindAllMarkers(sce,min.pct = 0.15)
save(out_marker, sce, file= 'output6_merge_adataobs/out_marker_0.15pct.RData')
####filter based \logFC\ > 0.25 & adj pvalue < 0.05
#x <- subset(x, abs(avg_log2FC) > 0.25 & p_val_adj < 0.05 )


###2021 11 13
###write out count matrix 
###
sce=subset(scetotal[, rownames(obsdata)])
head(sce@meta.data)

sce <- AddMetaData(
  object = sce,
  metadata = obsdata$Defined_Cell_Subtype_function,
  col.name = 'Defined_Cell_Subtype_function'
)
head(sce@meta.data)
sce
######select out 
sce_mine= subset(sce[, rownames(subset(sce@meta.data,  orig.ident == 'mysce'))])
sce_mine
######
count.data <- GetAssayData(object = sce_mine[["RNA"]], slot = "counts")
count.data = as.matrix(count.data)
count.data[1:5,1:4]




write.csv(sce_mine@meta.data, file = 'output6_merge_adataobs/cell_annotation_sce_mine.csv')
write.csv(count.data, file = 'output6_merge_adataobs/cell_matrix_sce_mine.csv')