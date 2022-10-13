######merge 8 samples with public data
######
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata')
######
pubdata=readRDS('/data/Zhuxq/lung.GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.rds')
pubcellab=read.table('/data/Zhuxq/lung.GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt', header=T, sep='\t' )
pubclin=read.csv('/data/Zhuxq/lung.GSE131907/GSE131907_Lung_Cancer_Feature_Summary.csv', header=T)
######
#select primary samples
pubcellab=subset(pubcellab, Sample_Origin=='tLung')
rownames(pubcellab)= pubcellab$Index

pubclin=subset(pubclin, Tissue.origins== 'tLung')
rownames(pubclin)=pubclin$Samples
head(pubclin)

length(intersect(rownames(pubcellab), colnames(pubdata)))
pubdata=pubdata[, colnames(pubdata) %in% rownames(pubcellab)]
pubcellab=pubcellab[colnames(pubdata),] 
#######
#create seruat object
library(Seurat)

pubsce= CreateSeuratObject(counts= pubdata, 
	project = "pubsce", assay = "RNA",
    min.cells = 0, min.features = 0,
	meta.data= pubcellab)
#add clincal information into it
Idents(object = pubsce) <- 'Sample'
levels(pubsce)
length(intersect(levels(pubsce),  rownames(pubclin)))
pubclin=pubclin[levels(pubsce), ] #make the order is the same
#age
Idents(object = pubsce) <- 'Sample'
tarid <- pubclin$Age
names(x = tarid) <- levels(x = pubsce)
pubsce <- RenameIdents(object = pubsce, tarid)
pubsce[["age"]] <- Idents(object = pubsce)
#sex
Idents(object = pubsce) <- 'Sample'
tarid <- pubclin$Sex
names(x = tarid) <- levels(x = pubsce)
pubsce <- RenameIdents(object = pubsce, tarid)
pubsce[["Sex"]] <- Idents(object = pubsce)
#smoking
Idents(object = pubsce) <- 'Sample'
tarid <- pubclin$Smoking
names(x = tarid) <- levels(x = pubsce)
pubsce <- RenameIdents(object = pubsce, tarid)
pubsce[["Smoking"]] <- Idents(object = pubsce)
#Stage
Idents(object = pubsce) <- 'Sample'
tarid <- pubclin$Stages
names(x = tarid) <- levels(x = pubsce)
pubsce <- RenameIdents(object = pubsce, tarid)
pubsce[["Stages"]] <- Idents(object = pubsce)
#
Idents(object = pubsce) <- 'Sample'
levels(pubsce)

##
aa=pubsce@meta.data
aa$Cell_type.refined[is.na(aa$Cell_type.refined)]='nan'
aa$Cell_subtype[is.na(aa$Cell_subtype)]='nan'
aa$Cell_subtype=gsub('Undetermined','nan', aa$Cell_subtype)
##
pubsce<- AddMetaData(
  object = pubsce,
  metadata = aa$Cell_type.refined,
  col.name = 'Cell_type.refined'
)

pubsce <- AddMetaData(
  object = pubsce,
  metadata = aa$Cell_subtype,
  col.name = 'Cell_subtype'
)
#
pubsce[["percent.mt.org"]] <- PercentageFeatureSet(pubsce, pattern = "^MT-")
pubsce[["percent.rb.org"]] <- PercentageFeatureSet(pubsce, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
#delete nan cells
#######
Idents(object = pubsce) <- 'Cell_subtype'
levels(pubsce)

idx= unique(pubsce$Cell_subtype)[!(unique(pubsce$Cell_subtype) %in% c('nan'))]

pubsce= subset(pubsce, idents =idx)

#####
#####
#####
#####my sample
####
mysce =readRDS('/data/Zhuxq/young_LC_analysis/my_analysis/adata_orignal.combine8s.rds')
#load('/data/Zhuxq/young_LC_analysis/my_analysis3/output2_main_type/TotalTissue.combined.harmony_mainlienage.defined.RData')
#mysce= TotalTissue.combined.harmony

myclin=read.csv('/data/Zhuxq/young_LC_analysis/sample_info.csv', header=T)
rownames(myclin)=myclin$label
myclin
#add age infomation
Idents(object = mysce) <- 'patient'
levels(mysce)
myclin=myclin[levels(mysce),]

tarid <- myclin$age
names(x = tarid) <- levels(x = mysce)
mysce <- RenameIdents(object = mysce, tarid)
mysce[["age"]] <- Idents(object = mysce)

#
Idents(object = mysce) <- 'patient'
levels(mysce)
#
mysce[["percent.mt.org"]] <- PercentageFeatureSet(mysce, pattern = "^MT-")
mysce[["percent.rb.org"]] <- PercentageFeatureSet(mysce, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
#
#
###
###make the cell label names are the same
###
cel_lab_pub=data.frame(Index= pubsce$Index,
	Sample=pubsce$Sample,
	percent.mt=pubsce$percent.mt.org,
	percent.rb=pubsce$percent.rb.org,
	Cell_type= pubsce$Cell_type,
	Cell_type.refined= pubsce$Cell_type.refined,
	Cell_subtype=pubsce$Cell_subtype,
	age=pubsce$age,
	Sex= pubsce$Sex,
	Smoking= pubsce$Smoking,
	Stages= pubsce$Stages
	)
###
cel_lab_my=data.frame(Index= rownames(mysce@meta.data),
	Sample=mysce$patient,
	percent.mt=mysce$percent.mt.org,
	percent.rb=mysce$percent.rb.org,
	Cell_type= 'nan',
	Cell_type.refined= 'nan',
	Cell_subtype='nan',
	age=mysce$age
	)
#
#
#######
#######
#create sceruat object again

sce_pub= CreateSeuratObject(counts= pubsce[['RNA']]@counts, 
	project = "pubsce", assay = "RNA",
    min.cells = 0, min.features = 0,
	meta.data= cel_lab_pub)
#
sce_my= CreateSeuratObject(counts= mysce[['RNA']]@counts, 
	project = "mysce", assay = "RNA",
    min.cells = 0, min.features = 0,
	meta.data= cel_lab_my)

#
#
#now merge data
#
common.features <- intersect(rownames(sce_pub), rownames(sce_my))
length(common.features)
#
sce.combined <- merge(sce_pub[common.features, ], y = sce_my[common.features, ], project = "Lung.Pri.Tumor.tissue")
sce.combined

saveRDS(sce.combined, file='/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/total.tissue.combined.addMTRPS.percentage.rds')

#
#rownames(sce_pub)[rownames(sce_pub) %in% c('APITD1', 'MHF1')]
#
#library(clusterProfiler)
#eg1 = bitr(rownames(sce_pub), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#head(eg1)
#eg2 = bitr(rownames(sce_my), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#head(eg2)
#length(intersect(unique(eg1$ENTREZID), unique(eg2$ENTREZID)))







































