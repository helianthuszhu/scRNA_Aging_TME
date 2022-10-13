library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(dplyr)
library(Seurat)

setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org/output5_step1_cellchat')

outDir <- file.path("output_cellchat")
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)
#
for (i in 1: 3){
    #prepare sce data
	#ageindex= c('Aged', 'Young', 'intermediated')
	#ageindex= c('Young', 'intermediated')
	ageindex= c( 'intermediated')
    aa= readRDS(paste0(ageindex[i], '_aasel.rds'))
    #
    Idents(aa)='mainlineage2'
    levels(aa)
    #
    meta_target= aa@meta.data
    #
    #cellchat
    data.input= GetAssayData(aa, assay = "RNA", slot = "data") # normalized data matrix
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
future::plan("multiprocess", workers = 20) # do parallel
    
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
    saveRDS(cellchat, file = paste0(outDir, '/',ageindex[i], '_cellchat_output2.rds'))
    
}