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
