rm(list=ls())
#
statrf= read.csv('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/data/data_allFeature.csv',header = T, row.names = 1)
statrf=  statrf[,1:32]

head(statrf)

statrf= statrf[statrf$atb %in% c('no','yes'),]
statrf$age_group= ifelse(statrf$age> 60, 'old','young')
statrf$atb_age=paste(statrf$age_group, statrf$atb,sep='_')
table(statrf$atb_age)
###
table(statrf$atb_age,statrf$response_bin) %>% .[c(1,3),] %>% prop.table(margin = 1)


pval <- chisq.test(statrf$atb_age, statrf$response_bin)$p.value
g1=ggplot(statrf, aes(x=atb_age, fill=response_bin)) + 
  geom_bar(position = "fill")+scale_fill_manual(values = c("#E69F00","#56B4E9"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  #coord_flip()+
  ylab('Percentage (%)')+
  annotate("text", x=3, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)

g1
###
###
### 1
statrf= read.csv('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/data/data_allFeature.csv',header = T, row.names = 1)

statrf= statrf %>% select(c('motuID','age','atb','response_bin','Study')) %>% 
  filter(atb %in% c('no','yes'))
head(statrf)
### 2
load('~/Documents/data_from_macintel/pub.data.cellection/1.metagenome-DB/pub_data_list/190_PD1_all_wilcox.RData')
head(meta)
table(meta$Study)
###
ah= meta %>% filter(!(Study == 'McCulloch')) %>%
  select(c('Sample_ID', 'Age', 'Antibiotics', 'ICB_response','Study')) %>% 
  setNames(c('motuID','age','atb','response_bin','Study')) %>%
  mutate(response_bin=  ifelse(response_bin== 'R','responder','non_responder'))
head(ah)
dim(ah)
### 3
attty= read.csv('~/Documents/data_from_macintel/pub.data.cellection/1.metagenome-DB/pub_data_list/ICB-antibiotics-meta_ttty.csv',header = T)
head(attty)
table(attty$ProjectID)
attty= attty %>% filter(ProjectID %in% c('c770295','i541981')) %>% 
  select(c('SampleID', 'Age','Antibiotics','ResponseGroup','ProjectID')) %>%
  setNames(c('motuID','age','atb','response_bin','Study')) %>% 
  mutate(atb= ifelse(atb=='1','yes','no'),
         response_bin=  ifelse(response_bin== 'R','responder','non_responder')) %>%
  mutate(age= gsub('[>]','',age))
###
table(attty$response_bin)
###
statdata= rbind(statrf, ah,attty) %>% mutate(atbuse= sapply(stringr::str_split(atb, " "), `[`, 1)) %>%
  mutate(atbuse= ifelse(grepl('es',atbuse), 'yes','no')) %>%
  mutate(age= as.numeric(age)) %>%
  mutate(age_group= ifelse(age> 60, 'old','young') ) %>%
  mutate(atb_age= paste(age_group,atbuse,sep = '_')) %>%
  filter((Study %in% c('Derosa_etal_NM_2022','McCulloch_etal_NM_2022',
                       'Routy_etal_Science_2017.NSCLC','Routy_etal_Science_2017.RCC'#,
                       #'i541981','Zhang'
                       #'i541981','Hakozaki','c770295'
                       ))) #c770295, Hakozaki, i541981, Zhang

head(statdata)
table(statdata$Study, statdata$atbuse)
table(statdata$atb_age)
dim(statdata)
###
library(ggplot2)
atbyes= subset(statdata, atb=='yes')
atbno= subset(statdata, atb=='no')

pval <- chisq.test(atbyes$age_group, atbyes$response_bin)$p.value
g1=ggplot(atbyes, aes(x=age_group, fill=response_bin)) + 
  geom_bar(position = "fill")+scale_fill_manual(values = c("#E69F00","#56B4E9"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  #coord_flip()+
  ylab('Percentage (%)')+
  annotate("text", x=1, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)+ggtitle('ABT yes')#+theme(legend.position = "")

pval <- chisq.test(atbno$age_group, atbno$response_bin)$p.value
g2=ggplot(atbno, aes(x=age_group, fill=response_bin)) + 
  geom_bar(position = "fill")+scale_fill_manual(values = c("#E69F00","#56B4E9"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  #coord_flip()+
  ylab('Percentage (%)')+
  annotate("text", x=1, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)+ggtitle('ABT no')

#
out <- file.path(paste0('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s08_aging_atb_response/'))
if(!dir.exists(out)) dir.create(out,recursive=TRUE)

pdf(paste0(out, 'stacked_atb_age_response.pdf'),width = 7,height = 5)
cowplot::plot_grid(g2,g1,nrow=1)
dev.off()

write.csv(statdata, file= paste0(out, 'stacked_atb_age_response.csv'))
#



