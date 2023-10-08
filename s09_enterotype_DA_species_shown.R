rm(list=ls())
###
### DA of clusters
###
dasp= read.csv('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/DA_analysis_cohendis_wilcox_cluster1_vs_OtherCluster_SpeciesLevel.csv',
               header = T,row.names = 1) 

ids= dasp$metareaction_id[str_length(dasp$metareaction_id) != 23]
ids= ids[!grepl('incertae sedis',ids)]
##
dadraw= dasp %>% as.data.frame %>% 
  filter(metareaction_id %in% ids) %>%
  filter(!grepl('incertae',metareaction_id )) %>%
  filter(p_value < 0.05  ) %>% 
  select(c('metareaction_id', 'cohens_d','adjusted_p_value','Species')) %>%
  mutate(logp= -log10(adjusted_p_value)) %>%
  filter(!duplicated(Species)) %>%
  arrange(cohens_d) %>%  
  mutate(Species= gsub(' ','_',Species )) %>%
  mutate(Species=  factor(Species, levels=Species)) %>%
  mutate(EAE_DA= ifelse(cohens_d > 0, 'EAE_up','EAE_down')) 
head(dadraw)
##
write.csv(dadraw,'~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/test.tmp.csv')

dadraw= read.csv('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/test.tmp.csv', header = T) 
idsh=c('s__Akkermansia_muciniphila_','s__Alistipes_finegoldii_','s__Alistipes_putredinis_','s__Adlercreutzia_equolifaciens','s__Bacteroides_dorei/vulgatus_',
       's__Bacteroides_fragilis','s__Bacteroides_thetaiotaomicron','s__Bacteroides_uniformis','s__Bacteroides_sp._','s__Bifidobacterium_bifidum_',
       's__Bilophila_wadsworthia_','s__Blautia_obeum/wexlerae_','s__Coprobacillus_cateniformis_',#'s__Eggerthella_lenta_',
       's__Faecalibacterium_prausnitzii_','s__Flavonifractor_plautii','s__Parabacteroides_distasonis_','s__Parabacteroides_merdae_',
       's__Roseburia_hominis', 's__Roseburia_inulinivorans_','s__Ruminococcus_bromii_','s__Ruminococcus_bicirculans_','s__Streptococcus_salivarius',
       
       's__Clostridium_sp._CAG:448','s__Faecalibacterium_sp._','s__Prevotella_copri_','s__Clostridium_sp._CAG:433', 's__Clostridium_sp._CAG:798',
       's__Clostridium_sp._CAG:226','s__Clostridium_sp._CAG:127','s__Firmicutes_bacterium_ADurb.Bin467','s__Firmicutes_bacterium_CAG:272',
       's__Phascolarctobacterium_succinatutens',#'s__Hungatella_hathewayi_','s__Prevotella_stercorea_',
       's__Prevotella_lascolaii_','s__Prevotella_sp._CAG:279','s__Butyrivibrio_crossotus_'
      )
dadraw$nameshow= ifelse(dadraw$Species %in% idsh, dadraw$Species,NA)
dadraw= dadraw[order(dadraw$cohens_d,decreasing = T),]
dadraw$Species= factor(dadraw$Species, levels= dadraw$Species)
##
table(dadraw$EAE_DA)
##
p = ggplot(dadraw, aes(y = Species, x = cohens_d, label = nameshow)) +
  geom_col(aes(fill = EAE_DA), width=0.5) +
  coord_flip()+
  #geom_point() +
  geom_point(aes(color=EAE_DA))+
  scale_fill_manual(values = c('gray',"#E69F00"))+ 
  scale_color_manual(values=c('gray',"#E69F00"))+
  #coord_flip()+ 
  theme_classic()+
  #ylim(1, 5.5) +
  theme(
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.title.x = element_blank()
  ) +
  # xlim(1, 1.375) +
  geom_text_repel(
    seed = 233,
    size = 3.5,
    color = 'black',
    min.segment.length = 0,
    force = 2,
    force_pull = 2,
    box.padding = 0.1,
    max.overlaps = Inf,
    segment.linetype = 2, #线段类型,1为实线,2-6为不同类型虚线
    segment.color = 'black', #线段颜色
    segment.alpha = 0.7, #线段不透明度
    #nudge_x = 1 - abs(dadraw$cohens_d), #标签x轴起始位置调整
    nudge_x = 0.3, #标签x轴起始位置调整
    direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
    hjust = 0.5
  ) +
  ggtitle("EAE VS EAD")+
  theme(legend.position = "none",
        aspect.ratio = 0.3,
        #axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        #axis.text.x = element_blank(),
        plot.title = element_text(size=8,face="bold",hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.line = element_line(colour = "black"),
        axis.ticks.length=unit(.1, "cm"),
        axis.ticks = element_line(colour = "black", size = 0.5)
  )

pdf('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/bar_chart_daclusterenrichsigslected.pdf',width = 20,height = 5)
cowplot::plot_grid(p)
dev.off()
#######
####### box plot of selected feature
#######
load('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/data/s01_enterotype_public_cohort_Based.image_20230822.RData')
#######
dim(bact)
head(samples,n=2)
#######
boxdata= log2(bact+0.0001) %>% rownames_to_column('Sample_ID') %>%
  left_join(samples[, !grepl('PCoA', colnames(samples))] %>% select(c('Sample_ID','cluster')) ,by= 'Sample_ID') %>%
  column_to_rownames('Sample_ID')  %>% mutate(group= ifelse(cluster=='cluster1','EAE','EAD')) 
boxdata[1:4,1:4]
table(boxdata$cluster)
table(boxdata$group)

siggene= dadraw %>% filter(!is.na(nameshow))
siggene= siggene$metareaction_id
length(siggene)
boxtheme <-
  list(
    geom_boxplot(color = 'black',  width = 0.4, size = 0.5, fill = NA,outlier.shape = NA),
    geom_violin(color = 'black', scale = 'width',  linewidth = 0.5,trim = TRUE,alpha = 0.8),
    #geom_jitter(shape = 21, color = "gray", fill = "white", width = 0.1, size = 2, lwd = 0.05),
    #stat_compare_means(aes(group = cluster),label.x.npc = c("middle"),label.y.npc=c("top"),label = "p.signif",vjust = 1),
    #geom_violin(color = "black", draw_quantiles = 0.5, lwd = 0.5),
    #geom_jitter(shape = 21, color = "gray", fill = "white", width = 0.1, size = 1, lwd = 0.05),
    xlab(""),
    #scale_fill_muted(),
    #scale_color_muted(), 
    
    scale_fill_manual(values = c("#E69F00",'gray')),
    theme_minimal(),
    theme(legend.position = "none",
          #legend.justification = c(0,1),
          #legend.title = element_blank(),
          legend.text = element_text(size=2),
          legend.title = element_text(size=2),
          aspect.ratio = 1.2,
          #plot.margin = margin(t=5,r=20,b=5,l=20,unit="pt"),
          axis.title.y = element_text(size=2,colour="black"),
          axis.title.x = element_text(size=2,colour="black"),
          #axis.title.x = element_blank(),
          axis.text = element_text(size=2,colour="black",angle=0),
          axis.text.x = element_text(size=2,colour="black",angle=45,hjust=1),
          #axis.text.x = element_blank(),
          plot.title = element_text(size=2,face="bold",hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5),
          axis.line = element_line(colour = "black"),
          axis.ticks.length=unit(.1, "cm"),
          axis.ticks = element_line(colour = "black", size = 0.5)
    )
  )

plot_nsclc <- boxdata %>%
  dplyr::select(group, all_of(siggene)) %>%
  pivot_longer(-group) %>% 
  mutate(name=  fct_relevel(name, siggene)) %>%
  mutate(group= factor(group, levels= c('EAE','EAD'))) %>%
  #mutate(name= sapply(stringr::str_split(name, "[.]"), `[`, 1)) %>%
  ggplot(aes(x = group, y = value, fill = group) ) +
  facet_wrap(. ~ name, ncol = 11,scales = "free")+
  boxtheme#+
  #stat_compare_means(aes(group = group),label.x.npc = c("middle"),label.y.npc=c("top"),label = "p.signif",vjust = 1)

pdf('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/box_da_species_selected.shown.pdf',
      width= 8,height = 5)
cowplot::plot_grid(plot_nsclc)
dev.off()
############-------------------------------------------------------------------------------
############ up features
############-------------------------------------------------------------------------------
dadraw= read.csv('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/test.tmp.csv', header = T) 
idsh=c('s__Akkermansia_muciniphila_','s__Alistipes_finegoldii_','s__Alistipes_putredinis_','s__Bacteroides_dorei/vulgatus_','s__Adlercreutzia_equolifaciens',
       's__Bacteroides_fragilis','s__Bacteroides_thetaiotaomicron','s__Bacteroides_uniformis','s__Bacteroides_sp._','s__Bifidobacterium_bifidum_',
       's__Bilophila_wadsworthia_','s__Blautia_obeum/wexlerae_','s__Coprobacillus_cateniformis_',#'s__Eggerthella_lenta_',
       's__Faecalibacterium_prausnitzii_','s__Flavonifractor_plautii','s__Parabacteroides_distasonis_','s__Parabacteroides_merdae_',
       's__Roseburia_hominis', 's__Roseburia_inulinivorans_','s__Ruminococcus_bromii_','s__Ruminococcus_bicirculans_','s__Streptococcus_salivarius'#,
       
       
)
dadraw$nameshow= ifelse(dadraw$Species %in% idsh, dadraw$Species,NA)
dadraw= dadraw[order(dadraw$cohens_d,decreasing = T),]
dadraw$Species= factor(dadraw$Species, levels= dadraw$Species)

siggene= dadraw %>% filter(!is.na(nameshow))
siggene= siggene$metareaction_id
length(siggene)
##
data3= boxdata %>%
  dplyr::select(group, all_of(siggene)) %>%
  pivot_longer(-group) %>% 
  mutate(name=  fct_relevel(name, siggene)) %>%
  mutate(group= factor(group, levels= c('EAE','EAD'))) %>%
  mutate(sample= name)
  #mutate(sample= sapply(stringr::str_split(name, "[.]"), `[`, 1)) %>%
  #mutate(sample= gsub(' ','_',sample )) %>% as.data.frame
head(data3)
unique(data3$name)
####
####draw 1
####
mean2 <- summarySE(data3, measurevar = "value", groupvars = c("sample", "group"))
mean2
#mean2$sample <- factor(mean2$sample,levels = rev(data2$sample))
mean2$facet <- rep("log2(Rel Abu)",times=dim(mean2)[1])
head(mean2)
#绘图
p2 <- ggplot(mean2, aes(sample,value, color = group)) + 
  geom_errorbar(aes(ymin = value- se, ymax = value + se), 
                width = 0.1,position = position_dodge(0.8),linewidth=0.5) + 
  geom_point(position = position_dodge(0.8),shape=18,size=3)+
  #转变x轴与y轴位置
  #coord_flip()+
  #y轴范围设置
  #scale_y_continuous(limits = c(30, 110))+
  #主题设置
  theme_bw()+
  theme(legend.position = c(0.8,0.95),
        #aspect.ratio = 0.5,
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        #axis.text.x = element_text(color = "black",size=10),
        axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=10))+
  #标题位置
  #labs(y="Functional dispersion values",x=NULL,color=NULL)+
  scale_color_manual(values = c("#E69F00",'gray'))+
  geom_segment(x = 10.5, xend = 11.5, y = 78, yend = 78,color = "black", size = 0.8)+
  #添加分组矩形
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e")+
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6")+
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 7.5, xmax = 8.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 8.5, xmax = 9.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e")+
  annotate("rect", xmin = 9.5, xmax = 10.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 10.5, xmax = 11.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6")+
  annotate("rect", xmin = 12.5, xmax = 13.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 13.5, xmax = 14.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 14.5, xmax = 15.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e")+
  annotate("rect", xmin = 15.5, xmax = 16.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 16.5, xmax = 17.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 17.5, xmax = 18.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6")+
  annotate("rect", xmin = 18.5, xmax = 19.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 19.5, xmax = 20.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 20.5, xmax = 21.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e")+
annotate("rect", xmin = 21.5, xmax = 22.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") 
#facet_grid(~ facet)
#p2
pdf('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/box_da_species_selected.shown2_up.pdf',
    width= 15,height = 5)
cowplot::plot_grid(p2)
dev.off()
#####--------
##### draw 2
#####--------
daorder= dadraw %>% filter(!is.na(nameshow)) %>% arrange(desc(cohens_d))
dim(daorder)
unique(data3$sample)
head(data3)

data3$group= factor(data3$group, levels= c('EAD','EAE'))
#data3$sample=gsub('s__Bacteroides_sp','s__Bacteroides_sp._' ,data3$sample)

#data3$sample = factor(data3$sample, levels=aggregate(value~sample,data3,median) %>% arrange(desc(value)) %>% .[,'sample'])

data3$sample = factor(data3$sample, levels=daorder$metareaction_id)
head(data3)

p3=ggplot(data3, aes(x=group, y=value, group=group)) + 
  geom_boxplot(aes(fill=group),outlier.colour = "black",outlier.size = 0.01,width=0.7,size=0.3)+
  scale_y_continuous(breaks = seq(0, -15, by = -3))+
  stat_compare_means(label = "p.signif")+
  facet_grid(sample~. )+
  theme(aspect.ratio = 0.2,
        axis.text.y = element_blank(),
    strip.text.y = element_text(size=8,colour = "black", angle = 0),
        strip.background = element_rect(colour="black", fill="#9ecae1"))+
  scale_fill_manual(values= c('gray',"#E69F00"))+coord_flip()#+rotate()

pdf('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/box_da_species_selected.shown3_up.pdf',
    width= 8,height = 5)
cowplot::plot_grid(p3)
dev.off()
####
####
####
############-------------------------------------------------------------------------------
############ down features
############-------------------------------------------------------------------------------
dadraw= read.csv('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/test.tmp.csv', header = T) 
idsh=c(#'s__Akkermansia_muciniphila_','s__Alistipes_finegoldii_','s__Alistipes_putredinis_','s__Bacteroides_dorei/vulgatus_',
       #'s__Bacteroides_fragilis','s__Bacteroides_thetaiotaomicron','s__Bacteroides_uniformis','s__Bacteroides_sp._','s__Bifidobacterium_bifidum_',
       #'s__Bilophila_wadsworthia_','s__Blautia_obeum/wexlerae_','s__Coprobacillus_cateniformis_',#'s__Eggerthella_lenta_',
       #'s__Faecalibacterium_prausnitzii_','s__Flavonifractor_plautii','s__Parabacteroides_distasonis_','s__Parabacteroides_merdae_',
       #'s__Roseburia_hominis', 's__Roseburia_inulinivorans_','s__Ruminococcus_bromii_','s__Ruminococcus_bicirculans_','s__Streptococcus_salivarius'#,
       
       's__Clostridium_sp._CAG:448','s__Faecalibacterium_sp._','s__Prevotella_copri_','s__Clostridium_sp._CAG:433', 's__Clostridium_sp._CAG:798',
       's__Clostridium_sp._CAG:226','s__Clostridium_sp._CAG:127','s__Firmicutes_bacterium_ADurb.Bin467','s__Firmicutes_bacterium_CAG:272',
       's__Phascolarctobacterium_succinatutens',#'s__Hungatella_hathewayi_','s__Prevotella_stercorea_',
       's__Prevotella_lascolaii_','s__Prevotella_sp._CAG:279','s__Butyrivibrio_crossotus_'
)
dadraw$nameshow= ifelse(dadraw$Species %in% idsh, dadraw$Species,NA)
dadraw= dadraw[order(dadraw$cohens_d,decreasing = T),]
dadraw$Species= factor(dadraw$Species, levels= dadraw$Species)

siggene= dadraw %>% filter(!is.na(nameshow))
siggene= siggene$metareaction_id
length(siggene)
##
data3= boxdata %>%
  dplyr::select(group, all_of(siggene)) %>%
  pivot_longer(-group) %>% 
  mutate(name=  fct_relevel(name, siggene)) %>%
  mutate(group= factor(group, levels= c('EAE','EAD'))) %>%
  mutate(sample= name)
#mutate(sample= sapply(stringr::str_split(name, "[.]"), `[`, 1)) %>%
#mutate(sample= gsub(' ','_',sample )) %>% as.data.frame
head(data3)
unique(data3$name)
####
####draw 1
####
mean2 <- summarySE(data3, measurevar = "value", groupvars = c("sample", "group"))
mean2
#mean2$sample <- factor(mean2$sample,levels = rev(data2$sample))
mean2$facet <- rep("log2(Rel Abu)",times=dim(mean2)[1])
head(mean2)
#绘图
p2 <- ggplot(mean2, aes(sample,value, color = group)) + 
  geom_errorbar(aes(ymin = value- se, ymax = value + se), 
                width = 0,position = position_dodge(0.8),linewidth=0.5) + 
  geom_point(position = position_dodge(0.8),shape=18,size=3)+
  #转变x轴与y轴位置
  #coord_flip()+
  #y轴范围设置
  #scale_y_continuous(limits = c(30, 110))+
  #主题设置
  theme_bw()+
  theme(legend.position = c(0.8,0.95),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        #axis.text.x = element_text(color = "black",size=10),
        axis.text.x = element_text(size=8,colour="black",angle=45,hjust=1),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=10))+
  #标题位置
  #labs(y="Functional dispersion values",x=NULL,color=NULL)+
  scale_color_manual(values = c("#E69F00",'gray'))+
  geom_segment(x = 10.5, xend = 11.5, y = 78, yend = 78,color = "black", size = 0.8)+
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e")+
  annotate("rect", xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 5.5, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6")+
  annotate("rect", xmin = 6.5, xmax = 7.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 7.5, xmax = 8.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 8.5, xmax = 9.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e")+
  annotate("rect", xmin = 9.5, xmax = 10.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 10.5, xmax = 11.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 11.5, xmax = 12.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6")+
  annotate("rect", xmin = 12.5, xmax = 13.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") 
#facet_grid(~ facet)
#p2
pdf('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/box_da_species_selected.shown2_down.pdf',
    width= 8,height = 4)
cowplot::plot_grid(p2)
dev.off()
#####--------
##### draw 2
#####--------
daorder= dadraw %>% filter(!is.na(nameshow)) %>% arrange(desc(cohens_d))
dim(daorder)
unique(data3$sample)
head(data3)

data3$group= factor(data3$group, levels= c('EAD','EAE'))
#data3$sample=gsub('s__Bacteroides_sp','s__Bacteroides_sp._' ,data3$sample)

#data3$sample = factor(data3$sample, levels=aggregate(value~sample,data3,median) %>% arrange(desc(value)) %>% .[,'sample'])

data3$sample = factor(data3$sample, levels=daorder$metareaction_id)
head(data3)
#data3 %>% group_by(group,name) %>% summarise_all(mean)
aggregate(value~group +name,data3,mean)
#
p3=ggplot(data3, aes(x=group, y=value, group=group)) + 
  geom_boxplot(aes(fill=group),outlier.colour = "black",outlier.size = 0.01,width=0.7,size=0.3)+
  #scale_y_continuous(breaks = seq(0, -15, by = -3))+
  stat_compare_means(label = "p.signif")+
  facet_grid(sample~. )+
  theme(aspect.ratio = 0.2,
        axis.text.y = element_blank(),
        strip.text.y = element_text(size=8,colour = "black", angle = 0),
        strip.background = element_rect(colour="black", fill="#9ecae1"))+
  scale_fill_manual(values= c('gray',"#E69F00"))+coord_flip()#+rotate()

pdf('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/box_da_species_selected.shown3_down.pdf',
    width= 8,height = 5)
cowplot::plot_grid(p3)
dev.off()


#
save.image(file= '~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s01_entertype_public_cohort_based/bar_chart_daclusterenrichsigslected.RData')












###

p2 <- ggplot(mean2, aes(sample,value, color = group)) + 
  geom_errorbar(aes(ymin = value- se, ymax = value + se), 
                width = 0,position = position_dodge(0.8),linewidth=0.5) + 
  geom_point(position = position_dodge(0.8),shape=18,size=3)+
  #转变x轴与y轴位置
  coord_flip()+
  #y轴范围设置
  scale_y_continuous(limits = c(30, 110))+
  #主题设置
  theme_bw()+
  theme(legend.position = c(0.8,0.95),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(color = "black",size=10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "grey", color = "transparent"),
        strip.text = element_text(color="black",size=10))+
  #标题位置
  labs(y="Functional dispersion values",x=NULL,color=NULL)+
  scale_color_manual(values = c("#7fc190","#efb684"))+
  #添加分组矩形
  annotate("rect", xmin = 0, xmax = 6.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e") +
  annotate("rect", xmin = 6.5, xmax = 20.5, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#bcced6") +
  annotate("rect", xmin = 20.5, xmax = 23, ymin = -Inf, ymax = Inf, alpha = 0.2,fill="#01665e")+
  #手动添加显著性标记，图中列出部分，具体根据个人数据进行调整
  annotate('text', label = '**', x =11, y =80, angle=-90, size =5,color="black")+
  geom_segment(x = 10.5, xend = 11.5, y = 78, yend = 78,color = "black", size = 0.8)+
  facet_grid(~ facet)#基于分面函数添加图顶部标题
p2


########
boxdata %>%  ggplot(
    aes(x = group, y = `s__Akkermansia muciniphila .[ref_mOTU_v3_03591]`)) +
  labs(y= 's__Akkermansia muciniphila') +
  #geom_boxplot(color = 'black',  width = 0.4, size = 0.8, fill = NA,outlier.shape = NA)+
  geom_violin(aes(fill = group),color = 'black', draw_quantiles = 0.5,scale = 'width',  linewidth = 0.5,trim = TRUE,alpha = 0.8)+
  #geom_violin(color = "black", draw_quantiles = 0.5, lwd = 0.5),
  #geom_jitter(shape = 21, color = "gray", fill = "white", width = 0.1, size = 3, lwd = 0.05)+
  #geom_violin(aes(fill = group),color = "black", draw_quantiles = 0.5, lwd = 0.5)+
  #geom_point( position=position_jitter(width=0.1),color='grey4',position=position_jitterdodge(),,alpha=1,size=3)+
  #geom_point(position = 'jitter',color = 'grey',size = 2,alpha = 0.8) +
  #geom_jitter(shape = 21, color = "black", fill = "white", width = 0.15, size = 3, lwd = 0.3 )+
  #geom_jitter(shape = 21, color = "black", fill = "white", width = 0.1, size = 0.5, lwd = 0.05)+
  #theme_base()+ 
  theme_minimal()+theme_bw()+
  stat_compare_means()+
  scale_fill_manual(values = )+
  
  theme(
    aspect.ratio = 1.2,
    #plot.margin = margin(t=5,r=20,b=5,l=20,unit="pt"),
    axis.ticks = element_line(),
    axis.line = element_line(),
    legend.position = "none",
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text = element_text(color = "black")
  )
#
#



themeDivClust <-
  list(
    geom_violin(color = "black", draw_quantiles = 0.5, lwd = 0.5),
    #geom_jitter(shape = 21, color = "black", fill = "white", width = 0.15, size = 0.75, lwd = 0.3 ),
    xlab("Clusters"),
    scale_fill_muted(),
    scale_color_muted(), 
    theme_minimal(),
    theme(axis.ticks = element_line(),
          axis.line = element_line(),
          legend.position = "none",
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text = element_text(color = "black")
    )
  )

boxdata %>%
  ggplot(aes(x = group, y = `s__Akkermansia muciniphila .[ref_mOTU_v3_03591]`, fill = group) ) +
  ylab("AKK") +
  themeDivClust




p = ggplot(dadraw, aes(y = Species, x = cohens_d, label = nameshow)) +
  geom_col(aes(fill = EAE_DA), width=0.5) +
  #coord_flip()+
  #geom_point() +
  geom_point(aes(color=EAE_DA))+
  scale_fill_manual(values = c('gray',"#E69F00"))+ 
  scale_color_manual(values=c('gray',"#E69F00"))+
  #coord_flip()+ 
  theme_classic()+
  #ylim(1, 5.5) +
  theme(
    axis.line.x  = element_blank(),
    axis.ticks.x = element_blank(),
    #axis.text.x  = element_blank(),
    axis.text.y  = element_blank(),
    axis.title.x = element_blank(),
    aspect.ratio = 1.2
  ) +
  # xlim(1, 1.375) +
  geom_text_repel(
    seed = 233,
    size = 3.5,
    color = 'black',
    min.segment.length = 0,
    force = 2,
    force_pull = 2,
    box.padding = 0.1,
    max.overlaps = Inf,
    segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
    segment.color = 'black', #线段颜色
    segment.alpha = 0.5, #线段不透明度
    #nudge_x = 1 - abs(dadraw$cohens_d), #标签x轴起始位置调整
    nudge_x = 0.5, #标签x轴起始位置调整
    direction = "x", #按y轴调整标签位置方向，若想水平对齐则为x
    hjust = 0
  ) +
  ggtitle("EAE VS EAD")

pdf('~/Desktop/test_volcano.pdf',width = 10,height = 5)
cowplot::plot_grid(p)
dev.off()


ggplot(dadraw, aes(x = Species, y = cohens_d)) +
  geom_col() +
  geom_text(aes(label = Species), vjust = -0.2)+
  coord_flip()

# Basic barplot
p<-ggplot(data=dadraw, aes(x=cohens_d, y=Species)) +
  #geom_bar(stat="identity")+
  geom_bar(mapping = aes(x = Species, y = cohens_d, fill = EAE_DA), 
           show.legend = FALSE,
           stat = "identity", 
           position = "dodge")+
  scale_fill_manual(values = c('gray',"#E69F00"))+ 
  coord_flip()+ theme_classic()+
  geom_text(data = dadraw,
            aes(x = Species, y = cohens_d, label = nameshow, color = EAE_DA)
  )
p



p <- ggplot(dadraw, aes(factor(Species), cohens_d)) +
  geom_point(aes(color=EAE_DA))+
  scale_fill_manual(values = c('gray',"#E69F00"))+ 
  scale_color_manual(values=c('gray',"#E69F00"))+
  #coord_flip()+ 
  theme_classic()+
  stat_summary(
    fill = "gray90",
    colour = "black", 
    fun = "mean",
    geom = "col"
  )+
  stat_summary(
    aes(label = nameshow), 
    fun = "mean",
    geom = "text_repel",
    min.segment.length = 0, # always draw segments
    position = position_nudge_repel(y = -2)
  ) +
  labs(title = "position_nudge_repel()")
pdf('~/Desktop/test_volcano.pdf',width = 15,height = 5)
cowplot::plot_grid(p)
dev.off()



####
####--- volcano
####
library(dplyr)
library(ggplot2)
library(ggrepel)
#
ids= dasp$metareaction_id[str_length(dasp$metareaction_id) != 23]
ids= ids[!grepl('incertae sedis',ids)]

df= dasp  %>% 
  filter(metareaction_id %in% ids) %>%
  filter(!grepl('incertae',metareaction_id )) %>%
  filter(!duplicated(Species)) 

df$log2FoldChange=  df$cohens_d
df$pvalue =  df$p_value
#df$Symbol =  df$metareaction_id
df$Symbol =  df$Species

head(df)
#
#阈值确定：
pvalue = 0.05
log2FC = 0.1
#根据阈值添加上下调分组标签：
df$group <- case_when(
  df$log2FoldChange > log2FC & df$pvalue < pvalue ~ "up",
  df$log2FoldChange < -log2FC & df$pvalue < pvalue ~ "down",
  TRUE ~ 'none'
)
table(df$group)
#
#转换为因子，指定绘图顺序；
df$'-log10(pvalue)' <- -log10(df$pvalue) #新增-log10p列
df$group <- factor(df$group, levels = c("up","down","none"))
head(df)
#
#自定义颜色：
mycol <- c("#EB4232","#2DB2EB","#d8d8d8")
#自定义主题：
mytheme <- theme_classic() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.margin = margin(15,5.5,5.5,5.5))
#
#ggplot2绘制火山图：
p <-  ggplot(data = df,
             aes(x = log2FoldChange,
                 y = -log10(pvalue),
                 color = group)) + #建立映射
  geom_point(size = 2.2) + #绘制散点
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) + #自定义散点颜色
  scale_x_continuous(limits = c(-2, 2),
                     breaks = seq(-2, 2, by = 0.2)) + #x轴限制
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 40),
                     breaks = seq(0, 40, by = 10)) + #y轴限制
  geom_hline(yintercept = c(-log10(pvalue)),size = 0.7,color = "black",lty = "dashed") + #水平阈值线
  geom_vline(xintercept = c(-log2FC, log2FC),size = 0.7,color = "black",lty = "dashed") + #垂直阈值线
  mytheme
p
#分别筛选上下调中显著性top20，作为本次测试需要添加的目标标签（共40个标签）：
up <- filter(df, group == 'up') %>% distinct(Symbol, .keep_all = T) %>%
  top_n(20, -log10(pvalue))
up = df %>% filter(group == 'up') %>% filter(Symbol %in% c('s__Alistipes obesi',
                                                           's__Bacteroides uniformis',
                                                           's__Akkermansia muciniphila ',
                                                           's__Bacteroides sp. ',
                                                           's__Adlercreutzia equolifaciens'
                                                           )
                                             )
up
  
down <- filter(df, group == 'down') %>% distinct(Symbol, .keep_all = T) %>%
  top_n(20, -log10(pvalue))
head(up);head(down)
###
p5 <- p +
  geom_point(data = up,
             aes(x = log2FoldChange, y = -log10(pvalue)),
             color = '#EB4232', size = 4.5, alpha = 0.2) +
  geom_text_repel(data = up,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol),
                  seed = 233,
                  size = 3.5,
                  color = 'black',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.1,
                  max.overlaps = Inf,
                  segment.linetype = 3, #线段类型,1为实线,2-6为不同类型虚线
                  segment.color = 'black', #线段颜色
                  segment.alpha = 0.5, #线段不透明度
                  nudge_x = 1 - up$log2FoldChange, #标签x轴起始位置调整
                  direction = "y", #按y轴调整标签位置方向，若想水平对齐则为x
                  hjust = 0 #对齐标签：0右对齐，1左对齐，0.5居中
  )
p5
pdf('~/Desktop/test_volcano.pdf',width = 15,height = 5)
cowplot::plot_grid(p5)
dev.off()


###
###
library(ggplot2)
library(ggrepel)
library(stringr) 
library(dplyr)
Cost=structure(list(Row.Labels = structure(c(1L, 2L, 6L, 9L, 4L, 3L, 
                                        5L, 7L, 8L), .Label = c("Change the way P is applied", "Improve management of manure", 
                                                                "In channel measures to slow flow", "Keep stock away from watercourses", 
                                                                "No till trial ", "Reduce runoff from tracks and gateways", "Reversion to different vegetation", 
                                                                "Using buffer strips to intercept pollutants", "Water features to intercept pollutants"
                                        ), class = "factor"), Average.of.FS_Annual_P_Reduction_Kg = c(0.11, 
                                                                                                      1.5425, 1.943, 3.560408144, 1.239230769, 18.49, 0.091238043, 
                                                                                                      1.117113762, 0.11033263), Average.of.FS_._Change = c(0.07, 0.975555556, 
                                                                                                                                                           1.442, 1.071692763, 1.212307692, 8.82, 0.069972352, 0.545940711, 
                                                                                                                                                           0.098636339), Average.of.Cost_Per_Kg_P_Removal.undiscounted..LOW_Oncost = c(2792.929621, 
                                                                                                                                                                                                                                      2550.611429, 964.061346, 9966.056875, 2087.021801, 57.77580744, 
                                                                                                                                                                                                                                       165099.0425, 20682.62962, 97764.80805), Sum.of.Total_._Cost = c(358.33, 
                                                                                                                                                                                                                                                                                                       114310.49, 19508.2, 84655, 47154.23, 7072, 21210, 106780.34, 
                                                                                                                                                                                                                                                                                                       17757.89), Average.of.STW_Treatment_Cost_BASIC = c(155.1394461, 
                                                                                                                                                                                                                                                                                                                                                          155.1394461, 155.1394461, 155.1394461, 155.1394461, 155.1394461, 
                                                                                                                                                                                                                                                                                                                                                          155.1394461, 155.1394461, 155.1394461), Average.of.STW_Treatment_Cost_HIGH = c(236.4912345, 
                                                                                                                                                                                                                                                                                                                                                                                                                                         236.4912345, 236.4912345, 236.4912345, 236.4912345, 236.4912345, 
                                                                                                                                                                                                                                                                                                                                                                                                                                         236.4912345, 236.4912345, 236.4912345), Average.of.STW_Treatment_Cost_INTENSIVE = c(1023.192673, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             1023.192673, 1023.192673, 1023.192673, 1023.192673, 1023.192673, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             1023.192673, 1023.192673, 1023.192673)), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          -9L))



#

Cost %>%
  arrange(desc(Average.of.Cost_Per_Kg_P_Removal.undiscounted..LOW_Oncost)) %>%
  mutate(Row.Labels = forcats::fct_inorder(Row.Labels),
         cuml_reduc = cumsum(Average.of.FS_Annual_P_Reduction_Kg),
         bar_start  = cuml_reduc - Average.of.FS_Annual_P_Reduction_Kg,
         bar_center = cuml_reduc - 0.5*Average.of.FS_Annual_P_Reduction_Kg) 
  ggplot(aes(xmin = bar_start, xmax = cuml_reduc,
             ymin = 0, ymax = Average.of.Cost_Per_Kg_P_Removal.undiscounted..LOW_Oncost)) +
  geom_rect(fill = "grey", colour = "black") +
  geom_text_repel(aes(x = bar_center, 
                      y = Average.of.Cost_Per_Kg_P_Removal.undiscounted..LOW_Oncost,
                      label = str_wrap(Row.Labels, 15)), 
                  ylim = c(0, NA), xlim = c(0, NA),  ## EDIT
                  size = 3, nudge_y = 1E4, nudge_x = 2, lineheight = 0.7, 
                  segment.alpha = 0.3) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Measure code and average P reduction (kg/P/yr)",
       y = "Mean annual TOTEX (£/kg) of P removal (thousands)")+
  coord_flip()+ theme_classic()
#
#
cost <- `rownames<-`(Cost[-1], Cost[,1])
# logarithmize values
w <- log(w1 <- cost$Average.of.Cost_Per_Kg_P_Removal.undiscounted..LOW_Oncost)
# define vector labels inside / outside, at best by hand
inside <- as.logical(c(0, 1, 0, 1, 1, 0, 1, 1, 1))
# calculate `x0` values of labels
x0 <- w / 2 + c(0, cumsum(w)[- length(w)])
# define y values o. labels
y0 <- ifelse(inside, colSums(t(cost)) / 2, 1.5e5)
# make labels using 'strwrap' 
labs <- mapply(paste, strwrap(rownames(cost), 15, simplify=F), collapse="\n")
# define nine colors
colores <- hcl.colors(9, "Spectral", alpha=.7)

# the actual plot
b <- barplot(cs <- colSums(t(cost)), width=w, space=0, ylim=c(1, 2e5), 
             xlim=c(-1, 80), xaxt="n", xaxs="i", col=colores, border=NA,
             xlab="Measure code and average P reduction (kg/P/yr)",
             ylab="Mean annual TOTEX (£/kg) of P removal (thousands)")

# place lables, leave out # 6
text(x0[-6], y0[-6], labels=labs[-6], cex=.7)
# arrows
arrows(x0[c(1, 3)], 1.35e5, x0[c(1, 3)], cs[c(1, 3)], length=0)
# label # 6
text(40, 1e5, labs[6], cex=.7)
# arrow # 6
arrows(40, 8.4e4, x0[6], cs[6], length=0)
# make x axis
axis(1, c(0, cumsum(log(seq(0, 1e5, 1e4)[-1]))), 
     labels=format(c(0, cumsum(seq(0, 1e5, 1e4)[-1])), format="d"), tck=-.02)
# put it in a box
box()

