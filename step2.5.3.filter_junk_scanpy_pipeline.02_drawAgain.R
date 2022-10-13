#########mainlineage ROE
#########
#########
#########
#########

setwd('/data/Zhuxq/young_LC_analysis/my_analysis4_merge.Pubdata.org')

outDir <- "output2_step2.5.3.filter_junk_scanpy_pipeline/stat_tissue_distribution"

meta_non_mal_cell= read.csv(paste0(outDir, '/metadata_combined_cellSubtype_defined.csv'), header=T, row.names=1)

meta_mal_cell= read.csv('output2_step2.5.3.filter_junk_scanpy_pipeline/data_Epithelial/malignant/data/adata_malignantcell.csv',header=T, row.names=1)

#

meta1= meta_non_mal_cell[,c('Sample', 'age', 'Defined_Cell_Type')]

meta2= meta_mal_cell[,c('Sample', 'age')]
meta2= data.frame(meta2, Defined_Cell_Type='Malignant')
#

meta=rbind(meta1, meta2)
head(meta)

meta$age_num= as.numeric(meta$age)

meta$age_group=cut(
  meta$age_num,
  breaks = c(20, 49, 60,100),
  labels = c("Young(20-40)", "intermediated(41-60)", "Aged(61-80)")
)
#
library(dplyr)
options(stringsAsFactors=FALSE)
library(reticulate)
library(tidyr)
library(ggpubr)

ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

#
stat= meta
summary <- table(stat[,c('age_group','Defined_Cell_Type')])
roe <- as.data.frame(ROIE(summary))
roe$Group <- rownames(roe)
roe_long=roe %>% pivot_longer(cols= colnames(roe[,!(colnames(roe) %in% c('Group'))]),
                                                                     names_to= "Cell_type",
                                                                values_to = "ROE")
head(roe_long)

library(Polychrome)
set.seed(723451)
fifth <- createPalette(15, c("#00ffff", "#ff00ff", "#ffff00"), M=1000)

p1=ggdotchart(roe_long, x = "Cell_type", y = "ROE",
           color = "Group",                                # Color by groups
           palette = as.vector(fifth), # Custom color palette
           #palette =FLATUI_CLASS,
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 8,                                 # Large dot size
           group = "Group",
           label = round(roe_long$ROE,2),                        # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+ geom_hline(yintercept = 1, linetype = 2, color = "black")+#theme(legend.position='none')
ylab("Ro/e")+ggtitle("")+theme(axis.text.x = element_text(angle = 45,vjust = 1))

pdf(paste0(outDir, '/ROE.mainlineage.Age_group_malignant_cell.pdf'),width=6,height = 7)
p1
dev.off()



##############
##############
##############
