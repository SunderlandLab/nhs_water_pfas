library(visreg)
library(mgcv)
library(ggplot2)
library(reshape2)
library(NADA)
library(dplyr)

one<-read.csv('/Users/cindyhu/Documents/Research/PFAS/NHSdrinking water/NHS data/Serum-water pair/one_230_022818.csv',header=T, stringsAsFactors = F)
water_master_clean<-read.csv("/Users/cindyhu/Documents/Research/PFAS/NHSdrinking water/NHS data/Redone_addBranched/water_master_clean.csv",
                             header=T,
                             stringsAsFactors = F)
colnames(water_master_clean)
water_master_clean_r<-water_master_clean[,c(1,2,45,46,53)] # subset relevant columns
df<-left_join(water_master_clean_r,one,by='id')
df<-df[!is.na(df$id),]
df<-df[!df$id=='111914',] #does not have any covariates
names(df)<-tolower(names(df))
colnames(df)[colnames(df)=='pfos_flags']<-'npfos_flags'
colnames(df)[colnames(df)=='pfos_br_flags']<-'brpfos_flags'
df$notmove8890<-(df$gdtlat1988==df$gdtlat1990) & (df$gdtlong1988==df$gdtlong1990)
table(df$notmove8890) #14 ppl moved
df$res_yr[!df$notmove8890]<-'<2'
df$notmove8688<-(df$gdtlat1986==df$gdtlat1988) & (df$gdtlong1986==df$gdtlong1988)
table(df$notmove8688)
df$res_yr[!df$notmove8688&df$notmove8890]<-'2-4'
table(df$notmove7686)
df$res_yr[!df$notmove7686&df$notmove8688&df$notmove8890]<-'4-14'
df$res_yr[df$notmove7686&df$notmove8688&df$notmove8890]<-'>14'

#ROS imputation
ChemList<-c('pfoa','pfna',
            'npfos','brpfos','pfhxs')
df_imput<-df
for (i in 1:5)
{
  print(paste0('impute water concentration for ',ChemList[i]))
  df[[paste0('w',ChemList[i])]][is.na(df[[paste0('w',ChemList[i])]])]<-0
  chem_censor_model<-ros(obs=df[[paste0('w',ChemList[i])]],
                         censored=(!(df[[paste0(ChemList[i],'_flags')]]=="")|is.na(df[[paste0('w',ChemList[i])]])))
  imput<-as.data.frame(chem_censor_model$modeled)
  id<-df$id[order(df[[paste0('w',ChemList[i])]])]
  new.df<-as.data.frame(cbind(id,imput))
  df_imput<-merge(df_imput,new.df,by='id')
  colnames(df_imput)[dim(df_imput)[2]]<-paste0('w',ChemList[i],'_imput')
}

df_imput$seafood<-rowSums(cbind(df_imput$ctuna90dn,df_imput$dkfsh90dn,df_imput$ofish90dn,df_imput$shrim90dn),na.rm=T)
df_imput$res_yr<-factor(df_imput$res_yr, levels=c('<2','2-4','4-14','>14'))

#dichotomize tap water consumption
df_imput$tapbin<-0
df_imput$tapbin[df_imput$tapnum>=8]<-1 #median is 8 cups a day
#dichotomize residential duration
df_imput$resbin<-0
df_imput$resbin[df_imput$res_yr=='4-14' | df_imput$res_yr=='>14']<-1
#fill in missing body weight
df_imput$wt90f[df_imput$wt90f==0]<-df_imput$wt88f[df_imput$wt90f==0]
df_imput$wt90f[df_imput$wt90f==0]<-df_imput$wt86f[df_imput$wt90f==0]
df_imput$tapnum<-df_imput$taph88num+df_imput$tapn88num

