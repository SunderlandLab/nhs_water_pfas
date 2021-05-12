library(visreg)
library(mgcv)
library(ggplot2)
library(reshape2)
library(NADA)
library(dplyr)

#source(main.R)
#gam model
df_gam$tapbin<-as.factor(df_gam$tapbin)
df_gam<-df_imput[!is.na(df_imput$snpfos),]
fit_pfoa <-gam(log(spfoab)~ s(log(wpfoa_imput),fx=FALSE,bs='cr', k=4, by=tapbin),data=df_gam)
visreg(fit_pfoa, 'wpfoa_imput',by='tapbin',trans = exp, partial=T,overlay=TRUE)
summary(fit_pfoa) 

##remove influential point #64
fit_pfoa2 <-gam(log(spfoab)~ s(log(wpfoa_imput),fx=FALSE,bs='cr', k=4,by=tapbin),data=df_gam[-64,])
summary(fit_pfoa2)
visreg(fit_pfoa2, 'wpfoa_imput',by='tapbin',trans = exp, partial=T,overlay=TRUE)

fit_pfna <-gam(log(spfnab)~ s(log(wpfna_imput),fx=FALSE,bs='cr',k=4,by=tapbin),data=df_gam)
plot(fit_pfna)
summary(fit_pfna)
visreg(fit_pfna, 'wpfna_imput',by='tapbin', trans=exp, partial=T, overlay=TRUE)

fit_npfos <-gam(log(snpfosb)~ s(log(wnpfos_imput),fx=FALSE,bs='cr',k=4, by=tapbin),data=df_gam)
plot(fit_npfos)
summary(fit_npfos)
visreg(fit_npfos, 'wnpfos_imput',by='tapbin', trans=exp, partial=T,overlay=TRUE)

fit_brpfos <-gam(log(sbrpfosb)~ s(log(wbrpfos_imput),fx=FALSE,bs='cr',k=4,by=tapbin),data=df_gam)
plot(fit_brpfos)
summary(fit_brpfos)
visreg(fit_brpfos, 'wbrpfos_imput',by='tapbin', trans=exp, partial=T,overlay=TRUE)

fit_pfhxs <-gam(log(spfhxsb)~s(log(wpfhxs_imput),fx=FALSE,bs='cr',k=4,by=tapbin),data=df_gam)
plot(fit_pfhxs)
summary(fit_pfhxs)
visreg(fit_pfhxs, 'wpfhxs_imput',by='tapbin', trans=exp, partial=T,overlay=TRUE)

#Table 4 make predictions using unadjusted gam model
ChemList<-c('pfoa','pfoa2','pfna','npfos','brpfos','pfhxs')
ucmrList<-c(20,20,20,28,12,30)
for (i in 1:6){
  if (!i==2)
  {
    newd<-data.frame(median(df_gam[[paste0('w',ChemList[i],"_imput")]], na.rm=T),tapbin=1)
    colnames(newd)[1]<-paste0('w',ChemList[i],'_imput')
    newd2<-data.frame(ucmrList[i],tapbin=1)
    colnames(newd2)[1]<-paste0('w',ChemList[i],'_imput')}
  else
  {
    newd<-data.frame(wpfoa_imput=median(df_gam[['wpfoa_imput']], na.rm=T),tapbin=1)
    newd2<-data.frame(wpfoa_imput=ucmrList[i],tapbin=1)
  }
  p1<-predict.gam(eval(parse(text=paste0('fit_',ChemList[i]))),newd,type='link', se.fit=TRUE)
  p2<-predict.gam(eval(parse(text=paste0('fit_',ChemList[i]))),newd2,type='link', se.fit=TRUE)
  predicted<-exp(p2$fit-p1$fit)
  upr <- exp(p2$fit-p1$fit + (1.96 * sqrt(p1$se.fit^2+p2$se.fit^2)))
  lwr <- exp(p2$fit-p1$fit - (1.96 * sqrt(p1$se.fit^2+p2$se.fit^2)))
  cat(paste('The estimated change for', ChemList[i], 'is', round((predicted-1)*100,1), "%.\n",
            'The lower limit is', round((lwr-1)*100,1), "%.\n",
            'The upper limit is', round((upr-1)*100,1), "%.\n"))
}

#generate number for line 295
ChemList<-c('pfoa','pfoa2')
for (i in 1:2){
  if (!i==2)
  { 
    print(i)
    newd<-data.frame(median(df_gam[[paste0('w',ChemList[i],"_imput")]], na.rm=T),tapbin=1)
    colnames(newd)[1]<-paste0('w',ChemList[i],'_imput')
    newd2<-data.frame(10*newd$wpfoa_imput,tapbin=1)
    colnames(newd2)[1]<-paste0('w',ChemList[i],'_imput')}
  else
  {
    print(i)
    newd<-data.frame(wpfoa_imput=median(df_gam[['wpfoa_imput']], na.rm=T),tapbin=1)
    newd2<-data.frame(wpfoa_imput=10*newd$wpfoa_imput,tapbin=1)
  }
  p1<-predict.gam(eval(parse(text=paste0('fit_',ChemList[i]))),newd,type='link', se.fit=TRUE)
  p2<-predict.gam(eval(parse(text=paste0('fit_',ChemList[i]))),newd2,type='link', se.fit=TRUE)
  predicted<-exp(p2$fit-p1$fit)
  upr <- exp(p2$fit-p1$fit + (1.96 * sqrt(p1$se.fit^2+p2$se.fit^2)))
  lwr <- exp(p2$fit-p1$fit - (1.96 * sqrt(p1$se.fit^2+p2$se.fit^2)))
  cat(paste('The estimated change for', ChemList[i], 'is', round((predicted-1)*100,1), "%.\n",
            'The lower limit is', round((lwr-1)*100,1), "%.\n",
            'The upper limit is', round((upr-1)*100,1), "%.\n"))
}

pdf("plasma_water_gam_041218.pdf",height = 5, width = 6.5, onefile=TRUE, pointsize=12)
par(mfrow=c(2,3))
visreg(fit_pfoa, 'wpfoa_imput',by='tapbin',trans = exp, partial=T, overlay=T, legend=T,
       points=list(cex=1, pch=16), 
       main='PFOA', xlab='Tap water concentrations (ng/L)', ylab='Plasma concentrations (ng/mL)')
visreg(fit_pfoa2, 'wpfoa_imput',by='tapbin',trans = exp, partial=T, overlay=T,legend=F,
       points=list(cex=1, pch=16),
       main='PFOA (excluding one)', xlab='Tap water concentrations (ng/L)', ylab='Plasma concentrations (ng/mL)')
visreg(fit_pfna, 'wpfna_imput', by='tapbin',trans=exp, partial=T, overlay=T,legend=F,
       points=list(cex=1, pch=16),
       main='PFNA', xlab='Tap water concentrations (ng/L)', ylab='Plasma concentrations (ng/mL)')
visreg(fit_npfos, 'wnpfos_imput',by='tapbin', trans=exp, partial=T, overlay=T,legend=F,
       points=list(cex=1, pch=16),
       main='nPFOS', xlab='Tap water concentrations (ng/L)', ylab='Plasma concentrations (ng/mL)')
visreg(fit_brpfos, 'wbrpfos_imput', by='tapbin',trans=exp, partial=T, overlay=T,legend=F,
       points=list(cex=1, pch=16),
       main='brPFOS', xlab='Tap water concentrations (ng/L)', ylab='Plasma concentrations (ng/mL)')
visreg(fit_pfhxs, 'wpfhxs_imput',by='tapbin', trans=exp, partial=T, overlay=T,legend=F,
       points=list(cex=1, pch=16),
       main='PFHxS', xlab='Tap water concentrations (ng/L)', ylab='Plasma concentrations (ng/mL)')
dev.off()

#adjusted models, consider covariates
#PFOA
fit_pfoa3 <-gam(log(spfoab)~ s(log(wpfoa_imput),fx=FALSE,bs='cr',k=4, by=tapbin)+res_yr,data=df_gam[-64,])
summary(fit_pfoa3)
#visreg(fit_pfoa3, 'wpfoa_imput', by='res_yr', trans=exp, partial=T, overlay=T)
df_gam$parity[df_gam$npar90>0]<-TRUE
df_gam$parity[df_gam$npar90==0]<-FALSE
df_gam$brfd2[df_gam$brfd86=='none']<-FALSE
df_gam$brfd2[!df_gam$brfd86=='none']<-TRUE
df_gam$tapnum<-df_gam$taph88num+df_gam$tapn88num
fit_pfoa4 <-gam(log(spfoab)~ s(log(wpfoa_imput),fx=FALSE,bs='cr',k=4, by=tapbin)+res_yr+men90+parity+brfd2+age+white+wt90f+seafood+popc90dn,
                data=df_gam)
summary(fit_pfoa4)
visreg(fit_pfoa4, 'wpfoa_imput',by='tapbin',trans = exp, partial=T,overlay=TRUE)


fit_pfoa5 <-gam(log(spfoab)~ s(log(wpfoa_imput),fx=FALSE,bs='cr',k=4, by=tapbin)+res_yr+men90+parity+brfd2+age+white+wt90f+seafood+popc90dn,
                data=df_gam[-64,])
summary(fit_pfoa5)
# succint model, but for consistency with other compounds, not use
# fit_pfoa5<-gam(log(spfoab)~ s(log(wpfoa_imput),fx=FALSE,bs='cr',k=4)+res_yr+tapnum,
#               data=df_gam)
# plot.gam(fit_pfoa5,shade=T, seWithMean = T, residuals = T, pch=16, cex=0.8)
# summary(fit_pfoa5) 


#PFNA
fit_pfna2 <-gam(log(spfnab)~ s(log(wpfna_imput),fx=FALSE,bs='cr',k=4, by=tapbin)+res_yr+men90+parity+brfd2+age+white+wt90f+seafood+popc90dn,
                data=df_gam)
plot(fit_pfna2)
summary(fit_pfna2)
# fit_pfna3 <-gam(log(spfnab)~ s(log(wpfna_imput)),
#                data=df_gam)
# summary(fit_pfna3)
# visreg(fit_pfna3, 'wpfna_imput',trans = exp, partial=T, main='PFNA (adjusted)', xlab='Tap water concentrations (ng/L)', ylab='Serum concentrations (ng/mL)')

#nPFOS
fit_npfos2 <-gam(log(snpfosb)~ s(log(wnpfos_imput),fx=FALSE,bs='cr',k=4, by=tapbin)+res_yr+men90+parity+brfd2+age+white+wt90f+seafood+popc90dn,
                 data=df_gam)
summary(fit_npfos2)

#brpfos
fit_brpfos2 <-gam(log(sbrpfosb)~ s(log(wbrpfos_imput),fx=FALSE,bs='cr',k=4, by=tapbin)+res_yr+men90+parity+brfd2+age+white+wt90f+seafood+popc90dn,
                  data=df_gam)
summary(fit_brpfos2)

#PFHxS
fit_pfhxs2 <-gam(log(spfhxsb)~ s(log(wpfhxs_imput),fx=FALSE,bs='cr',k=4, by=tapbin)+res_yr+men90+parity+brfd2+age+white+wt90f+seafood+popc90dn,
                 data=df_gam)
summary(fit_pfhxs2)

#table 4 adjusted model
ChemList<-c('pfoa4','pfoa5','pfna2','npfos2','brpfos2','pfhxs2')
ucmrList<-c(20,20,20,28,12,30)
for (i in 1:6){
  chemname<-substr(ChemList[i],1,nchar(ChemList[i])-1)
  newd<-data.frame(median(df_gam[[paste0('w',chemname,"_imput")]], na.rm=T),
                   res_yr='4-14',
                   men90=1,
                   parity=FALSE,
                   brfd2=FALSE,
                   age=median(df_gam$age),
                   white=1,
                   wt90f=median(df_gam$wt90f),
                   tapbin=1,
                   seafood=median(df_gam$seafood),
                  popc90dn=median(df_gam$popc90dn, na.rm=T))
  colnames(newd)[1]<-paste0('w',chemname,'_imput')
  newd2<-data.frame(ucmrList[i],
                    res_yr='4-14',
                    men90=1,
                    parity=FALSE,
                    brfd2=FALSE,
                    age=median(df_gam$age),
                    white=1,
                    wt90f=median(df_gam$wt90f),
                    tapbin=1,
                    seafood=median(df_gam$seafood),
                    popc90dn=median(df_gam$popc90dn, na.rm = T))
  colnames(newd2)[1]<-paste0('w',chemname,'_imput')
  p1<-predict.gam(eval(parse(text=paste0('fit_',ChemList[i]))),newd,type='link', se.fit=TRUE)
  p2<-predict.gam(eval(parse(text=paste0('fit_',ChemList[i]))),newd2,type='link', se.fit=TRUE)
  predicted<-exp(p2$fit-p1$fit)
  upr <- exp(p2$fit-p1$fit + (1.96 * sqrt(p1$se.fit^2+p2$se.fit^2)))
  lwr <- exp(p2$fit-p1$fit - (1.96 * sqrt(p1$se.fit^2+p2$se.fit^2)))
  cat(paste('The estimated change for', ChemList[i], 'is', round((predicted-1)*100,1), "%.\n",
            'The lower limit is', round((lwr-1)*100,1), "%.\n",
            'The upper limit is', round((upr-1)*100,1), "%.\n"))
}

