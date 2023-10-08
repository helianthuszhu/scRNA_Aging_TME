#########
#########
load('~/Documents/shanghai/chen_lungcancer/clin_age_related/tmpdata20210628.RData')
head(clindata_2065s)
length(unique(clindata_2065s$Cohort))
#set age group
qqq= clindata_2065s
qqq$agebin=ifelse(qqq$Age>60, 'old','young')
summary(qqq$Age)
head(qqq)
dim(qqq)
###
library(ggplot2)
#pdf('~/Documents/shanghai/chen_lungcancer/clin_age_related/meta_age/stacked_plot_response_bin.pdf',width=4,height = 5)
table(clindata$age_group)
table(clindata_2065s$age_group)

table(clindata$Response_bin_t)
table(clindata_2065s$Response_bin_t)

clindata_2065s$age_group=cut(
  clindata_2065s$Age,
  breaks = c(5,20, 40, 60,80,100),
  labels = c("group1(<=20)",'group2(21-40)', "group3(41-60)", "group4(61-80)",'group5(>=81)')
)
table(clindata_2065s$age_group)
table(clindata_2065s$age_group, clindata_2065s$clinical_benefit)
pdf('~/Documents/shanghai/chen_lungcancer/reanalysis_afterNM_20230506/clin_meta_analysis//stacked_plot_response_bin_2065s.pdf',width=4,height = 5)

pval <- fisher.test(clindata_2065s$age_group, clindata_2065s$clinical_benefit)$p.value
ggplot(clindata_2065s, aes(x=age_group, fill=clinical_benefit)) + 
  ggtitle(paste0('No.',dim(clindata_2065s)[1]))+
  geom_bar(position = "fill")+scale_fill_manual(values = c("#41ab5d","#fe9929"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  #coord_flip()+
  ylab('percentage')+
  annotate("text", x=3, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)
dev.off()
####################
tt= clindata_2065s
tt$age_group=cut(
    tt$Age,
    breaks = c(1, 40, 60,80,100),
    labels = c('group2(<40)', "group3(41-60)", "group4(61-80)",'group5(>=81)')
  )
pval <- chisq.test(tt$age_group, tt$Response_bin_t)$p.value
ggplot(tt, aes(x=age_group, fill=Response_bin_t)) + 
  ggtitle(paste0('No.',dim(tt)[1]))+
  geom_bar(position = "fill")+scale_fill_manual(values = c("#41ab5d","#fe9929"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  #coord_flip()+
  ylab('percentage')+
  annotate("text", x=3, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)

####################
pdf('~/Documents/shanghai/chen_lungcancer/reanalysis_afterNM_20230506/clin_meta_analysis/Boxplot_response_age.pdf',width=4,height = 5)
ggplot(clindata_2065s)+
  aes(x = clinical_benefit, y = Age, fill = clinical_benefit) + # add color to boxes with fill
  geom_boxplot(varwidth = FALSE) + # vary boxes width according to n obs.
  #geom_jitter(alpha = 0.8, width = 0,size=0.5) + # adds random noise and limit its width
  #facet_wrap(~year) + # divide into 2 panels
  theme(legend.position = "none") +# remove legend
  scale_fill_manual(values =  c("#41ab5d","#fe9929"))+
  stat_compare_means(aes(group = clinical_benefit))+
  ggtitle(paste0('Age','\n','No.',dim(clindata_2065s)[1]))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())
dev.off()
###
###
cutdata= clindata_2065s
#cutdata= clindata
cutdata$age_group=cut(
  cutdata$Age,
  breaks = c(5,20, 40, 60,80,100),
  labels = c("group1(<=20)",'group2(21-40)', "group3(41-60)", "group4(61-80)",'group5(>=81)')
)
cutdata$age_group=cut(
  cutdata$Age,
  breaks = c(5, 40, 60,100),
  labels = c("group1(<=40)", "group3(41-60)", 'group5(>60)')
)

table(cutdata$age_group)
cutdata$Smoking= gsub('current','Current', cutdata$Smoking)
cutdata$Smoking= gsub('former','Former', cutdata$Smoking)
cutdata$Smoking= gsub('never','Never', cutdata$Smoking)
cutdata$Smoking=factor(cutdata$Smoking, levels = c('Never','Former', 'Current' ))


cutdata$Drug2= cutdata$Drug
cutdata$Drug2= gsub('PD1_CTLA-4', 'Combo', cutdata$Drug2)
cutdata$Drug2= gsub('CTLA-4', 'CTLA4', cutdata$Drug2)
cutdata$Drug2= gsub('PD-L1', 'PD-1/PDL-1', cutdata$Drug2)
cutdata$Drug2= gsub('PD1', 'PD-1/PDL-1', cutdata$Drug2)
cutdata$Drug2=factor(cutdata$Drug2, levels = c('CTLA4','PD-1/PDL-1','Combo'))
table(cutdata$Drug2)


cutdata$CancerType2=cutdata$CancerType
cutdata$CancerType2= gsub('GC', 'Esophagogastric', cutdata$CancerType2)

cutdata$CancerType3= gsub('Bladder', 'Others', cutdata$CancerType2)
cutdata$CancerType3= gsub('Brain', 'Others', cutdata$CancerType3)
cutdata$CancerType3= gsub('HNSCC', 'Others', cutdata$CancerType3)
cutdata$CancerType3= gsub('Others', 'Others_Bladder/Brain/HNSCC', cutdata$CancerType3)

cutdata$CancerType3=factor(cutdata$CancerType3, levels = c('Melanoma','NSCLC','RCC','Esophagogastric','Others_Bladder/Brain/HNSCC'))
table(cutdata$CancerType3)

table(cutdata$Smoking)
cutdata$TMB=log2(cutdata$Mutations)
#####
head(cutdata)
library(survival)
my.surv=Surv(cutdata$PFS_month, cutdata$PFS_status)
#multicox=coxph(my.surv ~ age_group+ Sex+  Drug2+ CancerType3+Smoking+TMB, data =  cutdata)
multicox=coxph(my.surv ~ age_group+ Sex+  Drug2+ CancerType3,#+Smoking,#+TMB, 
               data =  cutdata)
library(survminer)
pdf('~/Documents/shanghai/chen_lungcancer/reanalysis_afterNM_20230506/clin_meta_analysis/RF_clin_2065s_3.pdf',width = 6,height = 5)
ggforest(multicox, data = cutdata,main='meta-cohort(21)_2065s')
dev.off()
######
######Chowell data multicox
#######
val1= readxl::read_xlsx('~/Documents/shanghai/chen_lungcancer/clinical_validation_Chowell_2021/41587_2021_1070_MOESM3_ESM.xlsx', sheet = 1)
val2= readxl::read_xlsx('~/Documents/shanghai/chen_lungcancer/clinical_validation_Chowell_2021/41587_2021_1070_MOESM3_ESM.xlsx', sheet = 2)
#
val1= as.data.frame(val1)
val2=as.data.frame(val2)
val1$cohort= 'discovery'
val2$cohort= 'validation'
head(val1)
#
Chowell= rbind(val1,val2)
head(Chowell)
#
colnames(Chowell)[5]='chemo'
colnames(Chowell)[7]='gender'
colnames(Chowell)[9]='stage'
colnames(Chowell)[21]='MSI'

Chowell$chemo= ifelse(Chowell$chemo=='0','no','yes')
Chowell$gender= ifelse(Chowell$gender=='0', 'Female','Male')
Chowell$stage= ifelse(Chowell$stage=='0', 'I-III','IV')
Chowell$MSI= ifelse(Chowell$MSI=='0', 'MSS','MSI')
table(Chowell$MSI)

Chowell$Drug_class= factor(Chowell$Drug_class, levels= c("CTLA4", "PD1/PDL1", "Combo"))

Chowell$age_group=cut(
  Chowell$Age,
  breaks = c(5, 40, 60,80,100),
  labels = c("group1(<=40)", "group3(41-60)", 'group5(60-80)','group5(>80)')
)

Chowell$age_group=cut(
  Chowell$Age,
  breaks = c(5, 40, 60,100),
  labels = c("group1(<=40)", "group3(41-60)", 'group5(>60)')
)

my.surv=Surv(Chowell$PFS_Months, Chowell$PFS_Event)
multicox=coxph(my.surv ~ age_group+TMB+ chemo+ Albumin+NLR+ 
                 + Platelets+ FCNA+ BMI+ HED+  HGB + Cancer_Type+ HLA_LOH+gender+ Drug_class+stage+MSI, data =  Chowell)

multicox=coxph(my.surv ~ age_group+ gender+TMB+chemo+Cancer_Type+Drug_class+stage+MSI+ BMI+ Albumin#+NLR
               +Platelets+FCNA#+ HED+HGB+HLA_LOH
               , 
               data =  Chowell)

outDir <- file.path('~/Documents/shanghai/chen_lungcancer/reanalysis_afterNM_20230506/Chowell/')
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

pdf('~/Documents/shanghai/chen_lungcancer/reanalysis_afterNM_20230506/Chowell/RF_Chowell3.pdf',width = 6,height = 8)
ggforest(multicox, data = Chowell,main = 'Hazard ratio---Chowell_etal_1479s')
dev.off()
#cc=summary(multicox)
#cc

###############################
###############################
#加载所需要的包
library(ggplot2)
#install.packages('rms')
library(rms)  
# 对数据进行打包，整理
data= cutdata
dd <- datadist(data) #为后续程序设定数据环境
options(datadist='dd') #为后续程序设定数据环境
# 拟合模型
fit<- cph(Surv(PFS_month,PFS_status) ~ rcs(Age,4) +Sex+ Drug2+  CancerType3 ,data=data)  # 节点数设为4
# 非线性检验
# P<0.05为存在非线性关系
anova(fit)

# 查看HR预测表
# 看一下预测的HR所对因的age
HR<-Predict(fit, Age,fun=exp,ref.zero = TRUE)
head(HR)
#
# 绘图
g1=ggplot()+
  geom_line(data=HR, aes(Age,yhat),
            linetype="solid",size=1,alpha = 0.7,colour="#0070b9")+
  geom_ribbon(data=HR, 
              aes(Age,ymin = lower, ymax = upper),
              alpha = 0.1,fill="#0070b9")+
  theme_classic()+
  geom_hline(yintercept=1, linetype=2,size=1)+
  geom_vline(xintercept=60,size=1,color = '#d40e8c',linetype=2)+ #查表HR=1对应的age
  labs(title = "ICI data-2065", x="Age", y="HR (95%CI)")+
  scale_x_continuous(limits = c(20, 80), breaks = seq(20, 80, by = 20))
###
########nonICI data
########
# 对数据进行打包，整理
data= nonicidata
dim(nonicidata)
dd <- datadist(data) #为后续程序设定数据环境
options(datadist='dd') #为后续程序设定数据环境
# 拟合模型
#fit<- cph(Surv(OS_month,OS_status) ~ rcs(Age,4) +Sex.x+Smoking.Status+Race+Sample.Type+Cancer.Type.defined,data=data)  # 节点数设为4
fit<- cph(Surv(OS_month,OS_status) ~ rcs(Age,4) +Sex.x+Smoking.Status+Race+Sample.Type,data=data)  # 节点数设为4
# 非线性检验
# P<0.05为存在非线性关系
anova(fit)

# 查看HR预测表
# 看一下预测的HR所对因的age
HR<-Predict(fit, Age,fun=exp,ref.zero = TRUE)
head(HR)
#
# 绘图
g2=ggplot()+
  geom_line(data=HR, aes(Age,yhat),
            linetype="solid",size=1,alpha = 0.7,colour="#0070b9")+
  geom_ribbon(data=HR, 
              aes(Age,ymin = lower, ymax = upper),
              alpha = 0.1,fill="#0070b9")+
  theme_classic()+
  geom_hline(yintercept=1, linetype=2,size=1)+
  geom_vline(xintercept=60,size=1,color = '#d40e8c',linetype=2)+ #查表HR=1对应的age
  labs(title = paste0("nonICI data, n= ", dim(data)[1]), x="Age", y="HR (95%CI)")+
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, by = 20))
#
pdf('~/Documents/shanghai/chen_lungcancer/reanalysis_afterNM_20230506/clin_meta_analysis/RCS_age_vs_ICB_vs_nonICB.pdf', width = 10,height=5)
cowplot::plot_grid(g1,g2,ncol=2)
dev.off()
#################
#################
#################
####
head(clindata_2065s)
dim(clindata_2065s)
clindata_2065s$clinical_benefit2=ifelse(clindata_2065s$clinical_benefit=='clinical_nenefit','Yes','No')
#####response bin
a1= clindata_2065s[,c('Ptid','Age','clinical_benefit2')]
a2= Chowell[,c('SAMPLE_ID','Age','Response_bin_t')]
colnames(a2)=c('Ptid','Age','clinical_benefit2')
aaaaaa=rbind(a1,a2)
head(aaaaaa)
#ageindex= unique(clindata$Age)[order(unique(clindata$Age), decreasing = F)][-c(1,2,71:72)]
ageindex= seq(6,99,1)
orlist = list()
for (i in 1: length(ageindex)) {
  tmp1= clindata_2065s
  tmp1$group=cut(
    tmp1$Age,
    breaks = c(5,ageindex[i],100),
    labels = c("1low",'2high')
  )
  #
  library(questionr)
  #
  mm= table(tmp1$group, tmp1$Response_bin_t)
  #mm= table(tmp1$group, tmp1$clinical_benefit2)
  #mm=mm[,seq(dim(mm)[2],1)]
  #class(mm)
  mm
  o1= odds.ratio(mm)
  o1=as.data.frame(o1)
  #
  ss= data.frame(paste0(mm[1,2],'/',rowSums(mm)[1] ), paste0(mm[2,2],'/',rowSums(mm)[2] ),  mm[1,1], mm[1,2], mm[2,1],mm[2,2])
  colnames(ss)=c('RR_lowcutoff','RR_highcutoff', 'NR.low','R.low','NR.high','R.high')
  ss
  ##
  dd= cbind(ss, o1)
  rownames(dd)= ageindex[i]
  dd
  orlist[[i]]= dd
  
}
or_output= do.call(rbind, orlist)
or_output$agecutoff=as.numeric(paste( rownames(or_output)))
or_output$sig= ifelse(or_output$p < 0.05, 'significant','non significant')
or_output$OR=as.numeric(paste(gsub('Inf', 'NA', or_output$OR)))
or_output=subset(or_output, NR.low >4 & R.low>4 & NR.high >4 & R.high >4)
or_output_response_bin=or_output
#######
#pdf('~/Documents/shanghai/chen_lungcancer/clin_age_related/cutoff_explore/cor_cutoff_age_OR_Response_bin.pdf',width = 5,height = 5)
ggscatter(or_output, x = "agecutoff", y = "OR", color = "sig",add = "reg.line", 
          conf.int = T,cor.coef = T, cor.method = "pearson",
          palette = c("#41ab5d","#fe9929","#00AFBB", "#E7B800", "#FC4E07"),
          add.params = list(color = "black",fill = "lightgray"))+
  scale_x_continuous(breaks=seq(20, 90, 10), limits = c(20, 90))+
  scale_y_continuous(breaks=seq(0, round(max(or_output_response_bin$OR), digits = 0), 1),
                     limits = c(0, round(max(or_output_response_bin$OR),digits = 0)))+
  xlab('Age cutoff')+ ylab('OR')+ggtitle('cor_cutoff_age_OR_Response_bin')+
  geom_vline(xintercept = 46, linetype="dotted", color = "blue", size=1.5)+
  annotate(geom="text", x=46, y=1, label="Age=46",color="red")
#dev.off()
#######
tt=aaaaaa
tt$age_group=cut(
  tt$Age,
  breaks = c(5, 40, 60,80,100),
  labels = c("group1(<=40)", "group3(41-60)", 'group5(60-80)','group5(>80)')
)
pval <- chisq.test(tt$age_group, tt$clinical_benefit2)$p.value
ggplot(tt, aes(x=age_group, fill=clinical_benefit2)) + 
  ggtitle(paste0('No.',dim(tt)[1]))+
  geom_bar(position = "fill")+scale_fill_manual(values = c("#41ab5d","#fe9929"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Change axis line
        axis.line = element_line(colour = "black"),
        #panel.border = element_blank(),
        panel.background = element_blank())+
  #coord_flip()+
  ylab('percentage')+
  annotate("text", x=3, y=0.8, label=paste0("p-value: ","\n" ,signif(pval,4)),size=6)
#######
#######
#######pie chart for 2065s
#######
dim(cutdata)
#
table(cutdata$CancerType3)
#
df=as.data.frame(table(cutdata$CancerType))
colnames(df)=c('Product','Value')
library(dplyr)
df <- df %>% 
  mutate(end = 2 * pi * cumsum(Value)/sum(Value),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))
df$Label <- paste(df$Product, paste0(round(((df$Value/sum(df$Value))*100),2),"%"), sep="-")
df
library(ggforce) # for 'geom_arc_bar'
g1=ggplot(df) + 
  geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                   start = start, end = end, fill = Product)) +
  geom_text(aes(x = 1 * sin(middle), y = 1 * cos(middle), label = Label,
                hjust = hjust, vjust = vjust)) +
  coord_fixed() +
  scale_x_continuous(limits = c(-1.5, 1.5),  # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL) +
  scale_y_continuous(limits = c(-1, 1.1),    # Adjust so labels are not cut off
                     name = "", breaks = NULL, labels = NULL)+
  scale_fill_brewer(palette="Set3")
pdf('~/Documents/shanghai/chen_lungcancer/reanalysis_afterNM_20230506/clin_meta_analysis/pie_chart_2065s.pdf',width = 7,height = 7)
g1
dev.off()
