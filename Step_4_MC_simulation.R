#source('main.R')
#Monte Carlo simulations
nsim<-300
ChemList<-c('pfoa','pfna',#'pfda',
            'npfos','brpfos','pfhxs')
#GI tract absorption factor, unitless
AFList<-c(0.91,0.91,#0.9,
          0.91,0.91,0.91)
#Volume of distribution, unit: mL/kg
VdList_miu<-c(170, 243.1, #source (Thompson, 2010, Env Intl)
              #441.1, #https://echa.europa.eu/documents/10162/3e08ce73-4b1a-444f-ac15-2c130616341c
              230, 230, 213 #source (Thompson, 2010, Env Intl)
)
VdList_sigma<-c(1.7, 48.9, #source (Thompson, 2010, Env Intl)
                #441.1, #https://echa.europa.eu/documents/10162/3e08ce73-4b1a-444f-ac15-2c130616341c
                2.25, 2.25, 28 #source (Thompson, 2010, Env Intl)
)

#Avanasi, 2016
#miu=5.2453, sigma=0.3540

#half-life, unit:day
ThalfList_miu<-c(4.7,2.7,4.8,4.8,7.3)
ThalfList_sigma<-c(1.2,2.0,1.1,1.1,1.1)
df_mc<-list()
for (sim in 1:nsim){
  df_mc[[sim]]<-df_imput
  print(paste('MC simulation',sim))
  for (i in 1:5)
  {
    #print(ChemList[i])
    df_mc[[sim]][[paste0('Vd_',ChemList[i])]]<-rep(0,110)
    for (j  in 1:110){
      df_mc[[sim]][[paste0('Vd_',ChemList[i])]][j]<-rlnormTrunc(1,meanlog = log(VdList_miu[i]),sdlog = log(VdList_sigma[i]),
                                                                min = (VdList_miu[i]-3*VdList_sigma[i]),
                                                                max = (VdList_miu[i]+3*VdList_sigma[i]))
      #truncated at miu+/-3*sigma
      df_mc[[sim]][[paste0('Thalf_',ChemList[i])]][j]<-rlnormTrunc(1,meanlog = log(ThalfList_miu[i]),sdlog = log(ThalfList_sigma[i]),
                                                                   min = max(0.34,(ThalfList_miu[i]-3*ThalfList_sigma[i])),#spefcify for PFNA otherwise min is negative
                                                                   max = (ThalfList_miu[i]+3*ThalfList_sigma[i]))*365#unit from years to days
      #tap water ingestion rate
      if (df_mc[[sim]]$taph88[j]=='none'){
        df_mc[[sim]]$taph88num[j]<-0
      }  else if (df_mc[[sim]]$taph88[j]=='1 - 2'){
        df_mc[[sim]]$taph88num[j]<-runif(1,1,2)
      }  else if (df_mc[[sim]]$taph88[j]=='3 - 5'){
        df_mc[[sim]]$taph88num[j]<-runif(1,3,5)
      }  else if (df_mc[[sim]]$taph88[j]=='6 - 9'){
        df_mc[[sim]]$taph88num[j]<-runif(1,6,9)
      } else if (df_mc[[sim]]$taph88[j]=='10 or more'){
        df_mc[[sim]]$taph88num[j]<-11
      } else {
        df_mc[[sim]]$taph88num[j]<-0
      }
      
      if (df_mc[[sim]]$tapn88[j]=='none'){
        df_mc[[sim]]$tapn88num[j]<-0
      } else if (df_mc[[sim]]$tapn88[j]=='1 - 2'){
        df_mc[[sim]]$tapn88num[j]<-runif(1,1,2)
      } else if (df_mc[[sim]]$tapn88[j]=='3 - 5'){
        df_mc[[sim]]$tapn88num[j]<-runif(1,3,5)
      } else if (df_mc[[sim]]$tapn88[j]=='6 - 9'){
        df_mc[[sim]]$tapn88num[j]<-runif(1,6,9)
      } else if (df_mc[[sim]]$tapn88[j]=='10 or more'){
        df_mc[[sim]]$tapn88num[j]<-11
      } else {
        df_mc[[sim]]$tapn88num[j]<-0
      }
      
      df_mc[[sim]]$tapnum<-df_mc[[sim]]$taph88num+df_mc[[sim]]$tapn88num
      
      #print(summary(df_mc[[sim]]$tapnum))
      
      if (df_mc[[sim]][j,paste0(ChemList[i],"_flags")]==""){
        C_water<-df_mc[[sim]][j,paste0('w',ChemList[i],'_imput')]
        df_mc[[sim]][j,paste0('sm',ChemList[i])]<-AFList[i]*df_mc[[sim]]$tapnum[j]*lpcup*C_water/(df_mc[[sim]]$wt90f[j]*kgplb*df_mc[[sim]][[paste0('Vd_',ChemList[i])]][j]*log(2)/df_mc[[sim]][[paste0('Thalf_',ChemList[i])]][j])
        if (!is.na(df_mc[[sim]][j,paste0('s',ChemList[i],'b')])) {
          df_mc[[sim]][j,paste0('r_',ChemList[i])]<-df_mc[[sim]][j,paste0('sm',ChemList[i])]/df_mc[[sim]][j,paste0('s',ChemList[i],'b')]
        }}
      else {
        df_mc[[sim]][j,paste0('sm',ChemList[i])]<-NA
        df_mc[[sim]][j,paste0('r_',ChemList[i])]<-NA
      }
    }
    #print(summary(df_mc[[sim]][paste0('r_',ChemList[i])]))
  }
}
#reference: http://www.milanor.net/blog/how-to-sort-a-list-of-dataframes-in-r/
sort_list_df <- function(df_l)
{
  out <- do.call(cbind, df_l)
  number_of_vars <- ncol(df_l[[1]])
  name_of_dfs <- names(out)
  out_l <- list()
  
  for(i in 1:number_of_vars)
  {
    index <- seq(i,ncol(out),number_of_vars)
    tempdf = out[, names(out) %in% name_of_dfs[index]]
    names(tempdf) <- names(df_l)
    out_l[[i]] = tempdf
  }
  
  names(out_l) <- names(df_l[[1]])
  return(out_l)
}

df_mc_sorted<-sort_list_df(df_mc)

#contribution to variability from three inout parameters, tapnum, Vd, and Thalf
MC_var_contrib<-{}
for (i in 1:5){
  chemname<-ChemList[i]
  print(chemname)
  
  cor1<-cor(unlist(df_mc_sorted[[paste0('r_',chemname)]]), unlist(df_mc_sorted[['tapnum']]), use='complete.obs')
  cor2<-cor(unlist(df_mc_sorted[[paste0('r_',chemname)]]), unlist(df_mc_sorted[[paste0('Vd_',chemname)]]), use='complete.obs')
  cor3<-cor(unlist(df_mc_sorted[[paste0('r_',chemname)]]), unlist(df_mc_sorted[[paste0('Thalf_',chemname)]]), use='complete.obs')
  MC_var_contrib$chem[i]<-chemname
  MC_var_contrib$tap_contrib[i]<-cor1^2/(cor1^2+cor2^2+cor3^2)
  MC_var_contrib$Vd_contrib[i]<-cor2^2/(cor1^2+cor2^2+cor3^2)
  MC_var_contrib$Thalf_contrib[i]<-cor3^2/(cor1^2+cor2^2+cor3^2)
}
MC_var_contrib

postscript('MC_TK_020218.eps',horizontal = FALSE, onefile = FALSE, paper = "special",
           width = 6.5,height = 4, pointsize = 9)
par(mfrow=c(2,3))
hist(sapply(df_mc_sorted[["r_pfoa"]],median,na.rm=TRUE), 
     main='PFOA',col='dodgerblue',border ="white",
     xlab='Median RSC', xlim=c(0,0.5))
hist(sapply(df_mc_sorted[["r_pfna"]],median,na.rm=TRUE), 
     main='PFNA',col='dodgerblue',border ="white",
     xlab='Median RSC', xlim=c(0,0.5))
hist(sapply(df_mc_sorted[["r_npfos"]],median,na.rm=TRUE), 
     main='nPFOS',col='dodgerblue',border =NA,
     xlab='Median RSC', xlim=c(0,0.5))
hist(sapply(df_mc_sorted[["r_brpfos"]],median,na.rm=TRUE), 
     main='brPFOS',col='dodgerblue',border =NA,
     xlab='Median RSC', xlim=c(0,0.5))
hist(sapply(df_mc_sorted[["r_pfhxs"]],median,na.rm=TRUE), 
     main='PFHxS',col='dodgerblue',border ="white",
     xlab='Median RSC', xlim=c(0,0.5))
dev.off()

#generating Table S6
tableS6<-matrix(0,nrow=10,ncol=5)
for (i in (1:5))
{
  chemname<-ChemList[i]
  print(chemname)
  tableS6[2*i-1,1]<-chemname
  tableS6[2*i-1,2]<-median(df_imput[[paste0('r_',chemname)]], na.rm=T)*100
  tableS6[2*i,1]<-chemname
  t1<-mean(sapply(df_mc_sorted[[paste0("r_",chemname)]],median,na.rm=TRUE))*100
  t2<-quantile(sapply(df_mc_sorted[[paste0("r_",chemname)]],median,na.rm=TRUE)*100,0.025)
  t3<-quantile(sapply(df_mc_sorted[[paste0("r_",chemname)]],median,na.rm=TRUE)*100,0.975)
  tableS6[2*i,2]<-paste0(round(t1,1)," (", round(t2,1), ", ", round(t3,1), ")")
  
  tableS6[2*i-1,3]<-mean(df_imput[[paste0('r_',chemname)]], na.rm=T)*100
  t1<-mean(sapply(df_mc_sorted[[paste0("r_",chemname)]],mean,na.rm=TRUE))*100
  t2<-quantile(sapply(df_mc_sorted[[paste0("r_",chemname)]],mean,na.rm=TRUE)*100,0.025)
  t3<-quantile(sapply(df_mc_sorted[[paste0("r_",chemname)]],mean,na.rm=TRUE)*100,0.975)
  tableS6[2*i,3]<-paste0(round(t1,1)," (", round(t2,1), ", ", round(t3,1), ")")
  
  tableS6[2*i-1,4]<-quantile(df_imput[[paste0('r_',chemname)]], 0.25, na.rm=T)*100
  t1<-mean(sapply(df_mc_sorted[[paste0("r_",chemname)]],quantile,0.25,na.rm=TRUE))*100
  t2<-quantile(sapply(df_mc_sorted[[paste0("r_",chemname)]],quantile,0.25,na.rm=TRUE)*100,0.025)
  t3<-quantile(sapply(df_mc_sorted[[paste0("r_",chemname)]],quantile,0.25,na.rm=TRUE)*100,0.975)
  tableS6[2*i,4]<-paste0(round(t1,1)," (", round(t2,1), ", ", round(t3,1), ")")
  
  tableS6[2*i-1,5]<-quantile(df_imput[[paste0('r_',chemname)]], 0.75, na.rm=T)*100
  t1<-mean(sapply(df_mc_sorted[[paste0("r_",chemname)]],quantile,0.75,na.rm=TRUE))*100
  t2<-quantile(sapply(df_mc_sorted[[paste0("r_",chemname)]],quantile,0.75,na.rm=TRUE)*100,0.025)
  t3<-quantile(sapply(df_mc_sorted[[paste0("r_",chemname)]],quantile,0.75,na.rm=TRUE)*100,0.975)
  tableS6[2*i,5]<-paste0(round(t1,1)," (", round(t2,1), ", ", round(t3,1), ")")
}

write.csv(tableS6,'tableS6.csv')

postscript('MC_contrib_var_020218.eps',horizontal = FALSE, onefile = FALSE, paper = "special",
           width = 6.5,height = 4, pointsize = 9)
par(mfrow=c(2,3),las=2) 
for (i in 1:5){
  y<-as.numeric(as.data.frame(MC_var_contrib)[i,-1])
  barplot(y, main=MC_var_contrib$chem[i], horiz=TRUE,
          names.arg=c("DW", "V_D", "Half-life"),col='grey50',border =NA,  xlab='Contribution to variability', xlim=c(0,1))
}
dev.off()
