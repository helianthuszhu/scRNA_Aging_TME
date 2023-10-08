mycolor2_v1= c("#E69F00","#56B4E9")
mycolor2_v2= c("#009E73",'#EF6F6A')
mycolor2_v3= c('#11A579FF',"#E69F00")

mycolor3= c("#E69F00","#56B4E9","#009E73")


mycolor4_v1= c("#E69F00","#56B4E9","#009E73",'#EEDD88')
mycolor4_v2 <- c('#FFAE34','#EF6F6A','#8CC2CA','#77AADD')
mycolor4_v3<- c('#77AADD','#EE8866','#EEDD88','#8CC2CA')
#khroma::colour("okabeito")(8)
mycolor8= c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
#
mycolor24= c("#7F3C8DFF",'#11A579FF','#3969ACFF','#F2B701FF','#E73F74FF','#80BA5AFF','#E68310FF','#008695FF', 
             '#CF1C90FF','#F97B72FF','#4B4B8FFF','#A5AA99FF','#B2DF8AFF','#33A02CFF','#FB9A99FF','#E31A1CFF',
             '#FDBF6FFF','#FF7F00FF','#CAB2D6FF','#6A3D9AFF','#FFFF99FF','#B15928FF','#A6CEE3FF','#1F78B4FF')

mycolor4_v4 = c('#FDBF6FFF','#A6CEE3FF', '#B2DF8AFF',"#F0E442" )

"#ED2891","#00A99D",

"#ED2891","#B2509E","#D49DC7","#C1A72F","#E8C51D","#F9ED32",
"#104A7F","#9EDDF9","#007EB5","#CACCDB","#6E7BA2","#DAF1FC","#00AEEF",
"#F6B667","#D97D25","#FBE3C7","#F89420","#97D1A9","#009444","#754C29","#CEAC8F",
"#3953A4","#BBD642","#00A99D","#D3C3E0","#A084BD","#542C88","#FAD2D9","#ED1C24",
"#F8AFB3","#EA7075","#7E1918","#BE1E2D","#b15928","#a6cee3","#b2df8a","#fb9a99","#fdbf6f",
"#cab2d6","#ffff99","#ff7f00","#33a02c","#6a3d9a","#1f78b4","#e31a1c"
#
#######
#######
#######
rm(list=ls())
#
load('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/data/s01_enterotype_public_cohort_Based.image_20230822.RData')
#
k
max_dim

clust_test <- hkmeans(data[ , 1:max_dim], k = k,
                      hc.metric = c("euclidean", "maximum", "manhattan", "canberra", "binary" ,"minkowski")[3]
                      )

pp1=fviz_dend(clust_test,type = c("rectangle", "circular", "phylogenic")[1],k_colors =mycolor3 ) 

samples$clust_test = as.factor(clust_test$cluster)


table(samples$clust_test,samples$clust)

table(samples$clust_test)
#fviz_dend(clust)


pp2= samples %>%
  #separate(azimut_id, into = c("allozithro_id", NA), sep = "_", remove = FALSE) %>% 
  ggplot(aes(x = PCoA1, y = PCoA2 )) +
  geom_encircle(aes(fill = clust_test, color = clust_test, group = clust_test), expand = 0, alpha = 0.2, size = 2) +
  geom_line(aes(group = Sample_ID), linetype = "dashed", size = 0.2, color = "darkgrey") +
  geom_point(aes(color = clust_test), size = 2) +
  scale_shape_manual(values = c(15, 18) ) +
  ggtitle('PCoA of enterotype clusters') +
  #scale_color_muted(name = "Cluster") +
  #scale_fill_muted(name = "Cluster") +
  #scale_fill_highcontrast(name = "Cluster")+
  #scale_color_highcontrast(name = "Cluster")+
  scale_color_manual(values=  mycolor3)+
  scale_fill_manual(values = mycolor3)+
  theme_minimal() +
  theme(legend.position = "right",
        #legend.justification = c(0,1),
        #legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8),
        aspect.ratio = 1,
        plot.margin = margin(t=5,r=20,b=5,l=20,unit="pt"),
        axis.title.y = element_text(size=8,colour="black"),
        axis.title.x = element_text(size=8,colour="black"),
        #axis.title.x = element_blank(),
        axis.text = element_text(size=8,colour="black",angle=0),
        axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1),
        #axis.text.x = element_blank(),
        plot.title = element_text(size=8,face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(.1, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.5)
  )
out= '~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/'
pdf(paste0(out,"plot_cluster_choose_with_dendrogram.pdf"), height = 5, width = 5)
print(optimal_cluster)
print(pp1)
print(pp2)
dev.off()
