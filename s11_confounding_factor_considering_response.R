##
## confounding factors
##
library("tidyverse")
library("coin")
##
load('~/Documents/project/mmm_aging/analysis_confrim_from_20230821/data/merged_clean_otumetatax_fiveICB_cohort.RData')
##
# general
set.seed(2022)
alpha.meta <- 1e-05
#
feat.all <- feat.rel.red
meta <- meta.crc

stopifnot(all(meta$Sample_ID %in% colnames(feat.all)))

studies <- meta %>% pull(Study) %>% unique

feat.all <- feat.all[,meta$Sample_ID]
#fn.pval <- paste0('../files/', tag, '/p_adj.tsv')
#if (!file.exists(fn.pval)){
#  stop("Please run the marker analysis script first. Exiting...\n")
#}
#p.vals.adj <- read.table(fn.pval, sep='\t', check.names = FALSE)
meta$age_bin= ifelse(meta$age> 60, 'old','young')
table(meta$Group)
meta$Group= meta$response_bin
##
table(meta$response_bin)
##
# ##############################################################################
# Group is the variable interested
#  variance explained by disease status
ss.disease <- apply(feat.all, 1, FUN=function(x, label){
  rank.x <- rank(x)/length(x)
  ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
  ss.o.i <- sum(vapply(unique(label), function(l){
    sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
  }, FUN.VALUE = double(1)))/length(rank.x)
  return(1-ss.o.i/ss.tot)
}, label=meta %>% pull(Group))

# calculate trimmed mean abundance
t.mean <- apply(feat.all, 1, mean, trim=0.1)

df.plot.all <- tibble(
  species=rownames(feat.all),
  disease=ss.disease,
  t.mean=t.mean#,
  #adj.p.val=p.vals.adj[rownames(feat.all), 'all'],
  #meta.significance=p.vals.adj[rownames(feat.all), 'all'] < 1e-05
)
#
meta$drug2= ifelse(meta$Study=='Lee_etal_NM_2022', 'comb','pd1')
#

# ##############################################################################
# Test all possible confounder variables
df.list <- list()
for (meta.var in c('age_bin','gender', 'bmi_bin', 'Study','atb','response_bin','PPI_use','country','drug2')){
  
  cat('###############################\n', meta.var, '\n')
  meta.c <- meta %>%
    filter(!is.na(eval(parse(text=meta.var))))
  
  cat('After filtering, the distribution of variables is:\n')
  print(table(meta.c$Group, meta.c %>% pull(meta.var)))
  print(table(meta.c$Study))
  feat.red <- feat.all[,meta.c$Sample_ID]
  
  cat('Calculating variance explained by meta-variable...\n')
  ss.var <- apply(feat.red, 1, FUN=function(x, label){
    rank.x <- rank(x)/length(x)
    ss.tot <- sum((rank.x - mean(rank.x))^2)/length(rank.x)
    ss.o.i <- sum(vapply(unique(label), function(l){
      sum((rank.x[label==l] - mean(rank.x[label==l]))^2)
    }, FUN.VALUE = double(1)))/length(rank.x)
    return(1 - ss.o.i/ss.tot)
  }, label=meta.c %>% pull(meta.var))
  df.plot.all[[meta.var]] <- ss.var
  
  cat('Calculating association with the meta-variable...\n')
  if (meta.c %>% pull(meta.var) %>% unique %>% length > 2){
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      kruskal.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  } else {
    meta.significance <- apply(feat.red, 1, FUN=function(x, var){
      wilcox.test(x~as.factor(var))$p.value
    }, var=meta.c %>% pull(meta.var))
  }
  meta.significance <- p.adjust(meta.significance, method='fdr')
  df.plot.all[[paste0(meta.var, '.significance')]] <- meta.significance
  cat('\n')
}
###
###
aa=as.data.frame( df.plot.all)
rownames(aa)= df.plot.all$species
###
summary(aa$disease)
g0=ggplot(aa,
          aes(x=disease, y=age_bin, size=t.mean+1e-08)) +
  #geom_point(shape=19) +
  geom_point(aes(size=t.mean+1e-08),alpha = 0.2,color = "steelblue")+
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  #facet_wrap(~facet, ncol=3) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))  +ggtitle('Age')

g1=ggplot(aa,
          aes(x=disease, y=gender, size=t.mean+1e-08)) +
  #geom_point(shape=19) +
  geom_point(aes(size=t.mean+1e-08),alpha = 0.2,color = "steelblue")+
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  #facet_wrap(~facet, ncol=3) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))  +ggtitle('Sex')

g2=ggplot(aa,
          aes(x=disease, y=atb, size=t.mean+1e-08)) +
  #geom_point(shape=19) +
  geom_point(aes(size=t.mean+1e-08),alpha = 0.2,color = "steelblue")+
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  #facet_wrap(~facet, ncol=3) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +ggtitle('atb use')

g3=ggplot(aa,
          aes(x=disease, y=bmi_bin, size=t.mean+1e-08)) +
  #geom_point(shape=19) +
  geom_point(aes(size=t.mean+1e-08),alpha = 0.2,color = "steelblue")+
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  #facet_wrap(~facet, ncol=3) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01))  +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))  +ggtitle('BMI')

g4=ggplot(aa,
          aes(x=disease, y=Study, size=t.mean+1e-08)) +
  #geom_point(shape=19) +
  geom_point(aes(size=t.mean+1e-08),alpha = 0.2,color = "steelblue")+
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  #facet_wrap(~facet, ncol=3) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01))  +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))  +ggtitle('Study')

g5=ggplot(aa,
          aes(x=disease, y=response_bin, size=t.mean+1e-08)) +
  #geom_point(shape=19) +
  geom_point(aes(size=t.mean+1e-08),alpha = 0.2,color = "steelblue")+
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  #facet_wrap(~facet, ncol=3) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01))  +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +ggtitle('ICB response')

g6=ggplot(aa,
          aes(x=disease, y=PPI_use, size=t.mean+1e-08)) +
  #geom_point(shape=19) +
  geom_point(aes(size=t.mean+1e-08),alpha = 0.2,color = "steelblue")+
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  #facet_wrap(~facet, ncol=3) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))  +ggtitle('PPI use')
g7=ggplot(aa,
          aes(x=disease, y=country, size=t.mean+1e-08)) +
  #geom_point(shape=19) +
  geom_point(aes(size=t.mean+1e-08),alpha = 0.2,color = "steelblue")+
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  #facet_wrap(~facet, ncol=3) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))  +ggtitle('country')

g8=ggplot(aa,
          aes(x=disease, y=drug2, size=t.mean+1e-08)) +
  #geom_point(shape=19) +
  geom_point(aes(size=t.mean+1e-08),alpha = 0.2,color = "steelblue")+
  xlab('Variance explained by disease status') +
  ylab('Variance explained by metadata variable') +
  theme_bw() +
  #facet_wrap(~facet, ncol=3) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0, 0.04), breaks = seq(0, 0.04, 0.01)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))  +ggtitle('drug')
##
out <- file.path("~/Documents/project/mmm_aging/analysis_confrim_from_20230821/output/s09_pub_metagdata_confoundingFactor_considering_response/")
if(!dir.exists(out)) dir.create(out,recursive=TRUE)
##
pdf(paste0(out, '/variance.metastudy_considering_response.pdf'),width = 12,height = 6)
cowplot::plot_grid(g0,g1,g2,g3,g4,g5,g6,g7,g8,ncol=3)
dev.off()
