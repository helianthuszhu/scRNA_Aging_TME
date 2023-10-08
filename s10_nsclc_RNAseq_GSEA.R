rm(list=ls())
dansclc= read.csv('Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s03_nsclc_cohort/DEG/degs_edgeR.csv',header = T, row.names = 1) %>%
  rownames_to_column('geneName')
head(dansclc)
#
library(clusterProfiler)
library(org.Hs.eg.db)
library(GseaVis)
#
head(dansclc)
########
########
dd=  dansclc
gene <- dd$geneName
## 转换
gene = bitr(gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
## 去重
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
gene_df <- data.frame(logFC=dd$logFC,
                      SYMBOL = dd$geneName)
gene_df <- merge(gene_df,gene,by="SYMBOL")
## geneList 三部曲
## 1.获取基因logFC
geneList <- gene_df$logFC
## 2.命名
names(geneList) = gene_df$ENTREZID
## 3.排序很重要
geneList = sort(geneList, decreasing = TRUE)
head(geneList)
##
########
########
## 读入hallmarks gene set，从哪来？
hallmarks <- read.gmt("~/Documents/project/mmm_aging/analysis_confrim_from_20230821/data/h.all.v2023.1.Hs.entrez.gmt")
# 需要网络
y <- GSEA(geneList,TERM2GENE =hallmarks)
write.csv(y %>% as.data.frame , file= '~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s03_nsclc_cohort/DEG/geneList.gsea.csv')
# all plot
geneSetID = c('HALLMARK_INFLAMMATORY_RESPONSE',
              'HALLMARK_INTERFERON_GAMMA_RESPONSE',
              #'HALLMARK_INTERFERON_ALPHA_RESPONSE',
              'HALLMARK_TNFA_SIGNALING_VIA_NFKB'
)
gseaNb(object = y,subPlot = 3,
       curveCol = c('#8CC2CA','#FFAE34','#77AADD'),
       geneSetID = geneSetID)

pdf('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s03_nsclc_cohort/DEG/gsea_plot.pdf',
    width = 5,height = 5)
gseaNb(object = y,
       geneSetID = geneSetID,
       curveCol = c("#E69F00","#56B4E9","#009E73"),
       subPlot = 3,
       termWidth = 100,
       legend.position = c(0.8,0.8),
       #addGene = gene,
       addPval = T,
       pvalX = 0.05,pvalY = 0.05)
dev.off()
