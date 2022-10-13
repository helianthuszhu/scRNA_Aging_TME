library(CellChat)
library(patchwork)
#
setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')

ddi= './output8_cellchat_s1/'

data.dir <- 'output8_cellchat_s1/combined_analysis/'
if(!dir.exists(data.dir)) dir.create(data.dir,recursive=TRUE)
#setwd(data.dir)
#
cellchat.aged= readRDS(paste0(ddi,'Aged_cellchat_output_ctypedetailed.rds'))
cellchat.Young= readRDS(paste0(ddi,'intermediated_cellchat_output_ctypedetailed.rds'))
cellchat.intermediated= readRDS(paste0(ddi,'Young_cellchat_output_ctypedetailed.rds'))
#
object.list <- list(cellchat.aged = cellchat.aged, cellchat.intermediated=cellchat.intermediated,cellchat.Young = cellchat.Young)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat
#
#
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
pdf(paste0(data.dir, '/1.No.interactions.pdf'))
gg1 + gg2
dev.off()
#
cid=c(
         
"C1-CD4-CCR7"    ,                          
 "C1-Monocyte-CD14"  ,       "C10-CD8-ZNF683" ,         
"C11-CD8-LAYN"  ,          
  "C12-CD8-TRDV1"  ,         
  "C2-CD4-HSPA1A"   ,        
      "C2-Monocyte-CD16"      ,  
       
 "C3-CD4-CXCR6"   ,                 
"C3-Mac-FABP4"    ,            
 "C4-CD4-CXCL13"      ,          
"C4-Mac-RETN"       ,    
       "C5-CD4-FOXP3"        ,    
        "C5-Mac-MARCO"     ,       
 "C6-CD4-RORC",                
 "C6-Mac-CCL13"  ,       
 "C7-CD8-SLC4A10"    ,            
"C7-Mac-CX3CR1",            "C8-CD8-GZMK",             
 "C8-Mac-SPP1"  ,            "C9-CD8-CX3CR1"        )   
    
#
levels(object.list[[1]]@idents)

tmacid= c(
3, 8, 10,12,15,17,27,31, 22,36,41,45,48,50,

7,24,29,33,38,43,47,49
)

tmacid= c(10,31,36,
          47,49,43)


weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf(paste0(data.dir, '/1.No.interactions.circle.pdf'), width=10)
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count[tmacid, tmacid], weight.scale = T, label.edge= T,
	  edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight[tmacid, tmacid], weight.scale = T, label.edge= T,
	  edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Interaction weights/strength - ", names(object.list)[i]))
}

dev.off()


######heatmap show interaction numbers
######
###############combine pheamap lists
library(pheatmap)
library(RColorBrewer)
hmcols=colorRampPalette(rev(brewer.pal(10, "RdBu")))(256)

## make breaks from combined range
myBreaks <- seq(min(c(object.list[[1]]@net$count[tmacid, tmacid], 
                      object.list[[2]]@net$count[tmacid, tmacid], 
                      object.list[[3]]@net$count[tmacid, tmacid])),
                max(c(object.list[[1]]@net$count[tmacid, tmacid], 
                      object.list[[2]]@net$count[tmacid, tmacid], 
	                  object.list[[3]]@net$count[tmacid, tmacid])), 
                length = 300)


ph1=  pheatmap(object.list[[1]]@net$count[tmacid, tmacid],breaks=myBreaks,cluster_cols=F,cluster_rows=F,color=hmcols, cellwidth=10,cellheight=10,border_color = NA)
ph2=  pheatmap(object.list[[2]]@net$count[tmacid, tmacid],breaks=myBreaks,cluster_cols=F,cluster_rows=F,color=hmcols, cellwidth=10,cellheight=10,border_color = NA)
ph3=  pheatmap(object.list[[3]]@net$count[tmacid, tmacid],breaks=myBreaks,cluster_cols=F,cluster_rows=F,color=hmcols, cellwidth=10,cellheight=10,border_color = NA)
plot_list=list()
plot_list[[1]] = ph1[[4]]
plot_list[[2]] = ph2[[4]]
plot_list[[3]] = ph3[[4]]
#g<-do.call(gridExtra::grid.arrange,plot_list)

g <- gridExtra::grid.arrange(gridExtra::arrangeGrob(grobs= plot_list,limitsize = FALSE,ncol=3))

ggsave(paste0(data.dir, '/1.heatmao.No.interactions.counts.pdf'),g,width =10,height = 10)

#-------------------------
myBreaks <- seq(min(c(object.list[[1]]@net$weight[tmacid, tmacid], 
                      object.list[[2]]@net$weight[tmacid, tmacid], 
                      object.list[[3]]@net$weight[tmacid, tmacid])),
                max(c(object.list[[1]]@net$weight[tmacid, tmacid], 
                      object.list[[2]]@net$weight[tmacid, tmacid], 
	                  object.list[[3]]@net$weight[tmacid, tmacid])), 
                length = 300)


ph1=  pheatmap(object.list[[1]]@net$weight[tmacid, tmacid],breaks=myBreaks,cluster_cols=F,cluster_rows=F,color=hmcols, cellwidth=10,cellheight=10,border_color = NA)
ph2=  pheatmap(object.list[[2]]@net$weight[tmacid, tmacid],breaks=myBreaks,cluster_cols=F,cluster_rows=F,color=hmcols, cellwidth=10,cellheight=10,border_color = NA)
ph3=  pheatmap(object.list[[3]]@net$weight[tmacid, tmacid],breaks=myBreaks,cluster_cols=F,cluster_rows=F,color=hmcols, cellwidth=10,cellheight=10,border_color = NA)
plot_list=list()
plot_list[[1]] = ph1[[4]]
plot_list[[2]] = ph2[[4]]
plot_list[[3]] = ph3[[4]]

ggsave(paste0(data.dir, '/1.heatmao.No.interactions.weight.pdf'),g,width =10,height = 4)
########
#c1= object.list[[1]]@net$count
#c1= c1[grepl('CD|Mac',rownames(c1)), grepl('CD|Mac',colnames(c1))]
#head(c1)
#
###
for (j in 1:length(object.list)) {
groupSize <- as.numeric(table(object.list[[j]]@idents))
mat <- object.list[[j]]@net$weight
#par(mfrow = c(3,4), xpd=TRUE)
    pdf(paste0(data.dir,'/',names(object.list)[j], '.pdf'))
for (i in 1:nrow(mat)) {

  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
    
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
   
}
 dev.off()
}

#############------------------------------------------------
#############------------------------------------------------
#############include Tregs
levels(object.list[[1]]@idents)
pdf(paste0(data.dir, '/1.Tcell2Mac_all_circle.pdf'))

netVisual_chord_gene(object.list[[1]], sources.use = c(10,31,36), title.name='Aged group',
                     targets.use = c(47,49,43), lab.cex = 0.5,legend.pos.y = 30)

#netVisual_chord_gene(object.list[[2]], sources.use = 10, targets.use = c(29,33,38,43,47,49), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(object.list[[3]], sources.use = c(10,31,36), title.name='Young group',
                     targets.use = c(47,49,43), lab.cex = 0.5,legend.pos.y = 30)

dev.off()


pdf(paste0(data.dir, '/1.Tcell2Mac_all_bubble.pdf'))

netVisual_bubble(object.list[[1]], sources.use = c(10,31,36), targets.use = c(47,49,43),remove.isolate = FALSE,
                title.name='Aged group')

netVisual_bubble(object.list[[3]], sources.use = c(10,31,36), targets.use = c(47,49,43),remove.isolate = FALSE,
                title.name='Young group')

dev.off()


pdf(paste0(data.dir, '/1.Mac2Tcell_all_circle.pdf'))

netVisual_chord_gene(object.list[[1]], sources.use = c(47,49,43), title.name='Aged group',
                     targets.use = c(10,31,36), lab.cex = 0.5,legend.pos.y = 30)

#netVisual_chord_gene(object.list[[2]], sources.use = 10, targets.use = c(29,33,38,43,47,49), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(object.list[[3]], sources.use = c(47,49,43), title.name='Young group',
                     targets.use = c(10,31,36), lab.cex = 0.5,legend.pos.y = 30)

dev.off()


pdf(paste0(data.dir, '/1.Mac2Tcell_all_bubble.pdf'))

netVisual_bubble(object.list[[1]], sources.use = c(47,49,43) , targets.use =c(10,31,36),remove.isolate = FALSE,
                title.name='Aged group')

netVisual_bubble(object.list[[3]], sources.use = c(47,49,43) , targets.use =c(10,31,36),remove.isolate = FALSE,
                title.name='Young group')

dev.off()

#####------------------------------------------------------------
#####------------------------------------------------------------
gg1 <- netVisual_bubble(cellchat, sources.use = c(47,49,43), targets.use = c(10,31,36),  comparison = c(1, 3), max.dataset = 1, title.name = "Increased signaling in Aged group", angle.x = 45, remove.isolate = T)+coord_flip() 
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(47,49,43) , targets.use =c(10,31,36),  comparison = c(1, 3), max.dataset = 3, title.name = "Decreased signaling in Aged group", angle.x = 45, remove.isolate = T)+coord_flip() 
#> Comparing communications on a merged object
pdf(paste0(data.dir, '/1.bubble_diff_signal_mac2Tcell.pdf'), width = 12, height=8)
cowplot::plot_grid(gg1,gg2,nrow= 2)
dev.off()

signaling.agedincreased_mac2Tcell=gg1$data
signaling.ageddecreased_mac2Tcell=gg2$data
save(signaling.agedincreased_mac2Tcell, signaling.ageddecreased_mac2Tcell, file = paste0(data.dir, '/1.bubble_diff_signal_mac2Tcell.RData'))
#####-----------------------------

gg1 <- netVisual_bubble(cellchat, sources.use = c(10,31,36), targets.use = c(47,49,43),  comparison = c(1, 3), max.dataset = 1, title.name = "Increased signaling in Aged group", angle.x = 45, remove.isolate = T)+coord_flip() 
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(10,31,36) , targets.use =c(47,49,43),  comparison = c(1, 3), max.dataset = 3, title.name = "Decreased signaling in Aged group", angle.x = 45, remove.isolate = T)+coord_flip() 
#> Comparing communications on a merged object
pdf(paste0(data.dir, '/1.bubble_diff_signal_Tcell2mac.pdf'), width = 12, height=8)
cowplot::plot_grid(gg1,gg2,nrow= 2)
dev.off()

signaling.agedincreased_Tcell2mac=gg1$data
signaling.ageddecreased_Tcell2mac=gg2$data
save(signaling.agedincreased_Tcell2mac, signaling.ageddecreased_Tcell2mac, file = paste0(data.dir, '/1.bubble_diff_signal_Tcell2mac.RData'))
#####
#####
#############------------------------------------------------
#############------------------------------------------------
#############exclude Tregs
levels(object.list[[1]]@idents)
pdf(paste0(data.dir, '/1.Tcell2Mac_all_circle_exclude.Tregs.pdf'))

netVisual_chord_gene(object.list[[1]], sources.use = c(10,31), title.name='Aged group',
                     targets.use = c(47,49,43), lab.cex = 0.5,legend.pos.y = 30)

#netVisual_chord_gene(object.list[[2]], sources.use = 10, targets.use = c(29,33,38,43,47,49), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(object.list[[3]], sources.use = c(10,31), title.name='Young group',
                     targets.use = c(47,49,43), lab.cex = 0.5,legend.pos.y = 30)

dev.off()


pdf(paste0(data.dir, '/1.Tcell2Mac_all_bubble_exclude.Tregs.pdf'))

netVisual_bubble(object.list[[1]], sources.use = c(10,31), targets.use = c(47,49,43),remove.isolate = FALSE,
                title.name='Aged group')

netVisual_bubble(object.list[[3]], sources.use = c(10,31), targets.use = c(47,49,43),remove.isolate = FALSE,
                title.name='Young group')

dev.off()


pdf(paste0(data.dir, '/1.Mac2Tcell_all_circle_exclude.Tregs.pdf'))

netVisual_chord_gene(object.list[[1]], sources.use = c(47,49,43), title.name='Aged group',
                     targets.use = c(10,31), lab.cex = 0.5,legend.pos.y = 30)

#netVisual_chord_gene(object.list[[2]], sources.use = 10, targets.use = c(29,33,38,43,47,49), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(object.list[[3]], sources.use = c(47,49,43), title.name='Young group',
                     targets.use = c(10,31), lab.cex = 0.5,legend.pos.y = 30)

dev.off()


pdf(paste0(data.dir, '/1.Mac2Tcell_all_bubble_exclude.Tregs.pdf'))

netVisual_bubble(object.list[[1]], sources.use = c(47,49,43) , targets.use =c(10,31),remove.isolate = FALSE,
                title.name='Aged group')

netVisual_bubble(object.list[[3]], sources.use = c(47,49,43) , targets.use =c(10,31),remove.isolate = FALSE,
                title.name='Young group')

dev.off()

#####------------------------------------------------------------
#####------------------------------------------------------------
gg1 <- netVisual_bubble(cellchat, sources.use = c(47,49,43), targets.use = c(10,31),  comparison = c(1, 3), max.dataset = 1, title.name = "Increased signaling in Aged group", angle.x = 45, remove.isolate = T)+coord_flip() 
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(47,49,43) , targets.use =c(10,31),  comparison = c(1, 3), max.dataset = 3, title.name = "Decreased signaling in Aged group", angle.x = 45, remove.isolate = T)+coord_flip() 
#> Comparing communications on a merged object
pdf(paste0(data.dir, '/1.bubble_diff_signal_mac2Tcell_exclude.Tregs.pdf'), width = 12, height=8)
cowplot::plot_grid(gg1,gg2,nrow= 2)
dev.off()

signaling.agedincreased_mac2Tcell=gg1$data
signaling.ageddecreased_mac2Tcell=gg2$data
save(signaling.agedincreased_mac2Tcell, signaling.ageddecreased_mac2Tcell, file = paste0(data.dir, '/1.bubble_diff_signal_mac2Tcell_exclude.Tregs.RData'))
#####-----------------------------

gg1 <- netVisual_bubble(cellchat, sources.use = c(10,31), targets.use = c(47,49,43),  comparison = c(1, 3), max.dataset = 1, title.name = "Increased signaling in Aged group", angle.x = 45, remove.isolate = T)+coord_flip() 
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(10,31) , targets.use =c(47,49,43),  comparison = c(1, 3), max.dataset = 3, title.name = "Decreased signaling in Aged group", angle.x = 45, remove.isolate = T)+coord_flip() 
#> Comparing communications on a merged object
pdf(paste0(data.dir, '/1.bubble_diff_signal_Tcell2mac_exclude.Tregs.pdf'), width = 12, height=8)
cowplot::plot_grid(gg1,gg2,nrow= 2)
dev.off()

signaling.agedincreased_Tcell2mac=gg1$data
signaling.ageddecreased_Tcell2mac=gg2$data
save(signaling.agedincreased_Tcell2mac, signaling.ageddecreased_Tcell2mac, file = paste0(data.dir, '/1.bubble_diff_signal_Tcell2mac_exclude.Tregs.RData'))





pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
pdf(paste0(data.dir, '/1.SPP1_circle.pdf'))
#par(mfrow = c(1,3), xpd=TRUE)

for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      sources.use = c(49), 
                      #targets.use = c(10,31,36),
                      signaling = pathways.show, layout = "chord",
                      signaling.name = paste(pathways.show, names(object.list)[i]))
    }
dev.off()


pdf(paste0(data.dir, '/1.SPP1_heatmap.pdf'),width=12)

pathways.show <- c("SPP1") 
par(mfrow = c(1,3), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, 
                               sources.use = c(49), 
                               targets.use = c(10,31,36),
                               color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]], ht_gap = unit(0.5, "cm"))

dev.off()






#####------------------------------

cid=c(
              
"C11-CD8-LAYN"  ,          
 "C4-CD4-CXCL13"      ,          
       "C5-CD4-FOXP3"        ,        
 "C6-Mac-CCL13"  ,       
"C7-Mac-CX3CR1",         
 "C8-Mac-SPP1"  ) 
#
cid= c( "C6-Mac-CCL13"  ,       
"C7-Mac-CX3CR1",         
 "C8-Mac-SPP1"  )
cellchat_tmac= subsetCellChat(cellchat, idents.use=cid)

cellchat_tmac@meta$datasets = factor(cellchat_tmac@meta$datasets, levels = c("cellchat.aged",'cellchat.intermediated', "cellchat.Young")) # set factor level

pdf(paste0(data.dir, '/1.violin.spp1.pdf'),width=4,height=4)
plotGeneExpression(cellchat_tmac, signaling = "SPP1",split.by = "datasets", colors.ggplot = T,
                   color=  c('#F0E716','#47EFFA','#EE2EE8')
                  )
dev.off()
###
###
cellchat.aged= readRDS(paste0(ddi,'Aged_cellchat_output_ctypedetailed.rds'))
#cellchat.Young= readRDS(paste0(ddi,'intermediated_cellchat_output_ctypedetailed.rds'))
cellchat.intermediated= readRDS(paste0(ddi,'Young_cellchat_output_ctypedetailed.rds'))
#
object.list <- list(cellchat.aged = cellchat.aged, 
                    #cellchat.intermediated=cellchat.intermediated,
                    cellchat.Young = cellchat.Young)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#> Merge the following slots: 'data.signaling','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.
cellchat

cid= c( "C6-Mac-CCL13"  ,       
"C7-Mac-CX3CR1",         
 "C8-Mac-SPP1"  )
cellchat_tmac= subsetCellChat(cellchat, idents.use=cid)

cellchat_tmac@meta$datasets = factor(cellchat_tmac@meta$datasets, levels = c("cellchat.aged",
                                                                             #'cellchat.intermediated', 
                                                                             "cellchat.Young")) # set factor level

pdf(paste0(data.dir, '/1.violin.spp1_2.pdf'),width=4,height=4)
plotGeneExpression(cellchat_tmac, signaling = "SPP1",split.by = "datasets", colors.ggplot = T,
                   color=  rev(c('#F0E716',
                             #'#47EFFA',
                             '#EE2EE8'))
                  )
dev.off()


#> The default behaviour of split.by has changed.
#> Separate violin plots are now plotted side-by-side.
#> To restore the old behaviour of a single split violin,
#> set split.plot = TRUE.
#>       
#> This message will be shown once per session.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.
#> Scale for 'y' is already present. Adding another scale for 'y', which will
#> replace the existing scale.



#
#object.list2 <- list(cellchat.Young = cellchat.Young,cellchat.aged = cellchat.aged)
#cellchat2 <- mergeCellChat(object.list2, add.names = names(object.list2))
#cellchat2
#weight.max2 <- getMaxWeight(object.list2, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
#pdf(paste0(data.dir, '/1.No.interactions.circle.AgedVSYoung.pdf'))
#par(mfrow = c(1,2), xpd=TRUE)
#netVisual_diffInteraction(cellchat2, weight.scale = T, measure = "count.merged", label.edge = T)
#netVisual_diffInteraction(cellchat2, weight.scale = T, measure = "weight.merged", label.edge = T)
#dev.off()
#
#
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf(paste0(data.dir, '/2.source_target_2D.pdf'),width=12,	height=4)

patchwork::wrap_plots(plots = gg)
dev.off()
#
#------------------------------------















