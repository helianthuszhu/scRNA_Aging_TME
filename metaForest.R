BiocManager::install('metasens')

library(meta)
library(metasens)
#########连续变量
data(Fleiss93cont)
Fleiss93cont
meta1 <-metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=Fleiss93cont, studlab=paste(study,year),sm="SMD")
meta1
forest(meta1)
metainf(meta1)

########二分类变量
data(Olkin95)
head(Olkin95)
meta1<- metabin(event.e, n.e, event.c, n.c, data=Olkin95[1:10,], sm="RR",studlab=paste(author, year))
summary(meta1)
forest(meta1)
########
########
library(metafor)

### copy BCG vaccine meta-analysis data into 'dat'
dat <- dat.bcg
dat
### calculate log risk ratios and corresponding sampling variances (and use
### the 'slab' argument to store study labels as part of the data frame)
dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat,
              slab=paste(author, year, sep=", "))

### fit random-effects model
res <- rma(yi, vi, data=dat)

### forest plot with extra annotations
forest(res, atransf=exp, at=log(c(.05, .25, 1, 4)), xlim=c(-16,6),
       ilab=cbind(tpos, tneg, cpos, cneg), ilab.xpos=c(-9.5,-8,-6,-4.5),
       cex=.75, header="Author(s) and Year", mlab="")
op <- par(cex=.75, font=2)
text(c(-9.5,-8,-6,-4.5), res$k+2, c("TB+", "TB-", "TB+", "TB-"))
text(c(-8.75,-5.25),     res$k+3, c("Vaccinated", "Control"))
par(op)

### add text with Q-value, dfs, p-value, and I^2 statistic
text(-16, -1, pos=4, cex=0.75, bquote(paste("RE Model (Q = ",
                                            .(formatC(res$QE, digits=2, format="f")), ", df = ", .(res$k - res$p),
                                            ", p = ", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, " = ",
                                            .(formatC(res$I2, digits=1, format="f")), "%)")))
##########
##########
head(Olkin95)
###########################
##########

load('~/Documents/shanghai/chen_lungcancer/clin_age_related/tmpdata20210628.RData')
head(clindata_2065s)
length(unique(clindata_2065s$Cohort))
#set age group
qqq= clindata_2065s
qqq= clindata
qqq$agebin=ifelse(qqq$Age>60, 'old','young')
summary(qqq$Age)
#
studyid= unique(qqq$Cohort)
studyid
countlist= list()
for (i in 1:length(studyid)) {
  w= subset(qqq, Cohort==studyid[i])
  tt= table(w$agebin,w$Response_bin_t) %>% as.vector %>% t %>% as.data.frame
  ee= cbind(studyid[i],tt)
  colnames(ee)= c('author','age_old_response_No','age_young_response_No',
                  'age_old_response_Yes','age_young_response_Yes')
  ee$totalOld=ee$age_old_response_No+ee$age_old_response_Yes
  ee$totalYoung=ee$age_young_response_No+ee$age_young_response_Yes

  countlist[[i]]= ee
}

rr= do.call(rbind, countlist)
rr$Author=sapply(str_split(rr$author, "_"), `[`, 1)
rr$CancerType=sapply(str_split(rr$author, "_"), `[`, 2)
rr$year=sapply(str_split(rr$author, "_"), `[`, 3)
rr= rr[order(rr$CancerType,decreasing = T),]
rr
metaICB<- metabin(event.e= age_old_response_Yes, 
                  n.e= totalOld, 
                  event.c= age_young_response_Yes, 
                  n.c= totalYoung,# byvar = CancerType,
                  #test.subgroup = T,
                  data=rr, sm="RR",studlab=paste(Author,year,CancerType)#,
                  #test.effect.subgroup.random = TRUE
                  )
summary(metaICB)
#pdf('~/Desktop/tmpforest.pdf',height = 20,width = 20)
forest(metaICB,layout = 'RevMan5',
       common = T,random = F)
#dev.off()
metabias(metaICB,method.bias="linreg")
funnel(metaICB)
####
####add Chowell cohort
#######
#######
val1= readxl::read_xlsx('~/Documents/shanghai/chen_lungcancer/clinical_validation_Chowell_2021/41587_2021_1070_MOESM3_ESM.xlsx', sheet = 1)
val2= readxl::read_xlsx('~/Documents/shanghai/chen_lungcancer/clinical_validation_Chowell_2021/41587_2021_1070_MOESM3_ESM.xlsx', sheet = 2)
#
val1= as.data.frame(val1)
val2=as.data.frame(val2)
val1$cohort= 'discovery'
val2$cohort= 'validation'
Chowell= rbind(val1,val2)
head(Chowell)
table(Chowell$`Response (1:Responder; 0:Non-responder)`)
Chowell$Response_bin_t= ifelse(Chowell$`Response (1:Responder; 0:Non-responder)`=='1','Yes','No')
Chowell$agebin=ifelse(Chowell$Age>60, 'old','young')
#summary(Chowell$Age)
tt= table(Chowell$agebin,Chowell$Response_bin_t) %>% as.vector %>% t %>% as.data.frame
ee= cbind('Chowell_etal_PanCancer_2018',tt)
ee
colnames(ee)= c('author','age_old_response_No','age_young_response_No',
                'age_old_response_Yes','age_young_response_Yes')
ee$totalOld=ee$age_old_response_No+ee$age_old_response_Yes
ee$totalYoung=ee$age_young_response_No+ee$age_young_response_Yes
ee
rr22=rbind(rr,ee)
rr22

metaICB22<- metabin(event.e= age_old_response_Yes, 
                  n.e= totalOld, 
                  event.c= age_young_response_Yes, 
                  n.c= totalYoung, 
                  data=rr22, sm="RR",studlab=paste(author))
summary(metaICB22)
forest(metaICB22,layout = 'RevMan5')
metabias(metaICB22,method.bias="linreg")
funnel(metaICB22)
#############
#############
#############
#https://www.metafor-project.org/doku.php/plots:forest_plot_bmj
library(metafor)
dat=rr
dat <- escalc(measure="RR", 
              ai=age_old_response_Yes, 
              n1i=totalOld, 
              ci=age_young_response_Yes, 
              n2i=totalYoung, data=dat,
              slab=paste(" ", author), addyi=FALSE)
### fit random-effects model (using DL estimator as RevMan does)
res <- rma(yi, vi, data=dat, method="FE")
res
forest(res)
############################################################################

### colors to be used in the plot
colp <- "#6b58a6"
coll <- "#a7a9ac"

### total number of studies
k <- nrow(dat)

### generate point sizes
psize <- weights(res)
psize <- 1.2 + (psize - min(psize)) / (max(psize) - min(psize))

### get the weights and format them as will be used in the forest plot
weights <- formatC(weights(res), format="f", digits=1)

### adjust the margins
par(mar=c(2.7,3.2,2.3,1.3), mgp=c(3,0,0), tcl=0.15)

### forest plot with extra annotations
sav <- forest(dat$yi, dat$vi, xlim=c(-3.4,2.1), ylim=c(-0.5,k+3),# alim=c(-1,1), 
              cex=0.88,
              pch=18, psize=psize, efac=0, refline=NA, lty=c(1,0), xlab="",
              ilab=cbind(paste(dat$ai, "/", dat$n1i), paste(dat$ci, "/", dat$n2i), weights),
              ilab.xpos=c(-1.9,-1.3,1.2), annosym=c(" (", " to ", ")"),
              rowadj=-.07)

### add vertical reference line at 0
segments(0, -1, 0, k+1.6, col=coll)

### add vertical reference line at the pooled estimate
segments(coef(res), 0, coef(res), k, col=colp, lty="33", lwd=0.8)

### redraw the CI lines and points in the chosen color

#segments(summary(dat)$ci.lb, k:1, summary(dat)$ci.ub, k:1, col=colp, lwd=1.5)

points(dat$yi, k:1, pch=18, cex=psize*1.15, col="white")
points(dat$yi, k:1, pch=18, cex=psize, col=colp)

### add the summary polygon
addpoly(res, row=0, mlab="Total (95% CI)", efac=2, col=colp, border=colp)

### add horizontal line at the top
abline(h=k+1.6, col=coll)

### redraw the x-axis in the chosen color
axis(side=1, at=seq(-1,1,by=0.5), col=coll, labels=FALSE)

### now we add a bunch of text; since some of the text falls outside of the
### plot region, we set xpd=NA so nothing gets clipped
par(xpd=NA)

### adjust cex as used in the forest plot and use a bold font
par(cex=sav$cex, font=2)

### add headings
text(sav$xlim[1], k+2.5, pos=4, "Study or\nsubgroup")
text(sav$ilab.xpos[1:2], k+2.3, c("Experimental","Control"))
text(mean(sav$ilab.xpos[1:2]), k+3.4, "No of events / total")
text(0, k+2.7, "Risk difference, IV,\nrandom (95% CI)")
segments(sav$ilab.xpos[1]-0.22, k+2.8, sav$ilab.xpos[2]+0.13, k+2.8)
text(c(sav$ilab.xpos[3],sav$xlim[2]-0.35), k+2.7, c("Weight\n(%)","Risk difference, IV,\nrandom (95% CI)"))

### add 'Favours caffeine'/'Favours decaf' text below the x-axis
text(c(-1,1), -2.5, c("Favors control","Favors experimental"), pos=c(4,2), offset=-0.3)

### use a non-bold font for the rest of the text
par(cex=sav$cex, font=1)

### add the 100.0 for the sum of the weights
text(sav$ilab.xpos[3], 0, "100.0")

### add the column totals for the counts and sample sizes
text(sav$ilab.xpos[1:2], 0, c(paste(sum(dat$ai), "/", sum(dat$n1i)), paste(sum(dat$ci), "/", sum(dat$n2i))))

### add text with heterogeneity statistics
text(sav$xlim[1], -1, pos=4, bquote(paste("Test for heterogeneity: ", tau^2, "=",
                                          .(formatC(res$tau2, digits=2, format="f")), "; ", chi^2, "=",
                                          .(formatC(res$QE, digits=2, format="f")), ", df=", .(res$k - res$p),
                                          ", P=", .(formatC(res$QEp, digits=2, format="f")), "; ", I^2, "=",
                                          .(formatC(res$I2, digits=0, format="f")), "%")))

### add text for test of overall effect
text(sav$xlim[1], -2, pos=4, bquote(paste("Test for overall effect: Z=",
                                          .(formatC(res$zval, digits=2, format="f")),
                                          ", P", .(ifelse(res$pval<.001, "<0.001",
                                                          paste0("=",formatC(res$pval, digits=2, format="f")))))))
pdf('~/Documents/shanghai/chen_lungcancer/reanalysis_afterNM_20230506/clin_meta_analysis/metaForest.pdf',width = 10,height = 10)
forest(metaICB,layout = 'RevMan5',common = T,random = F)
OEtext= bquote(paste("Test for overall effect: Z=",
                     .(formatC(res$zval, digits=2, format="f")),
                     ", P", .(ifelse(res$pval<.001, "<0.001",
                                     paste0("=",formatC(res$pval, digits=2, format="f"))))))
grid.text(OEtext, .275, .17, gp=gpar(cex=0.88))
dev.off()