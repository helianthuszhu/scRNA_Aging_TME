library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(dplyr)
library(Seurat)

setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')
outDir <- file.path("output5_step1_cellchat")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)
	
#load seruat data
sce= readRDS('output2_step2.5.3.filter_junk_scanpy_pipeline/stat_tissue_distribution/v2/sceall_sel.for.cellchat.rds')
Idents(sce)='mainlineage2'
levels(sce)
cname= levels(sce)[!(levels(sce) %in% c('mixed'))]
sce_sel=subset(sce, idents = cname)
#
ageindex= rev(unique(sce_sel$age_group))
ageindex
for (i in 1: length(unique(sce_sel$age_group))){
    #prepare sce data
    aa= sce_sel
    #
    Idents(aa)='mainlineage2'
    levels(aa)
    #aa= subset(aa, idents= c('CD4', 'CD8'))
    #
    #aasel= subset(aa, idents= ageindex[1])
    aasel= subset(aa, subset = age_group %in% ageindex[i])
	#save(aasel, file= paste0(outDir, '/',ageindex[i], '_aasel.rds') )
	
    #
    meta_target= aasel@meta.data
    #
    #cellchat
    data.input= GetAssayData(aasel, assay = "RNA", slot = "data") # normalized data matrix
	#write.table(as.matrix(data.input), file= paste0(outDir, '/',ageindex[i], '_aasel.seurat.data.txt'), quote=F, sep='\t')
	
    cellchat <- createCellChat(object = data.input, meta = meta_target, group.by = "mainlineage2")
    #
    cellchat <- setIdent(cellchat, ident.use = "mainlineage2") # set "labels" as default cell identity
    levels(cellchat@idents) # show factor levels of the cell labels
    groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
    #
    #Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    # use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB
    # set the used database in the object
cellchat@DB <- CellChatDB.use

    
#Preprocessing
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 10) # do parallel
    
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

#Inference of cell-cell communication network
#cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,type = "truncatedMean", trim = 0.1,population.size = TRUE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

#Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

#visualize the aggregated cell-cell communication network
groupSize <- as.numeric(table(cellchat@idents))

    #cell-cell communication network
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
###
    saveRDS(cellchat, file = paste0(outDir, '/',ageindex[i], '_cellchat_output.rds'))
    
}