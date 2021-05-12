library(reshape2)
library(ggplot2)
library(visreg)
library(mgcv)
library(lambda.tools)
library(EnvStats)

#setwd("Serum-water pair")
#source('main.R')
#TK model
#C_serum_m=AF*V_water*C_water/(bw*V_d*ln(2)/t_half)

ChemList<-c('pfoa','pfna',#'pfda',
            'npfos','brpfos','pfhxs')
#GI tract absorption factor, unitless
AFList<-c(0.91,0.91,#0.9,
          0.91,0.91,0.91)
#Volume of distribution, unit: mL/kg
VdList<-c(170, 243.1, #source (Thompson, 2010, Env Intl)
          #441.1, #https://echa.europa.eu/documents/10162/3e08ce73-4b1a-444f-ac15-2c130616341c
          230, 230, 213 #source (Thompson, 2010, Env Intl)
          )

#half-life, unit:day
ThalfList<-c(4.7*365,2.7*365,4.8*365,4.8*365,7.3*365)
            #1278, 913, #1642, 
             #1752, 1752, 2665)
#ThalfList<-c(2055, 1035, 3722, 1971, 1971, 2665)

lpcup<-0.237 # liter per cup=0.237L/cup
kgplb<-0.454 # kg per lb = 0.454kg/lb

# illustrate number in line 217
#pfoa average LOD is 0.6 ng/L
0.91*1*0.6/(70*VdList[1]*log(2)/ThalfList[1])/median(df_imput$spfoab, na.rm=T)
#illustrate number in line 326
#pfoa
0.91*1*median(df_imput$wpfoa, na.rm=T)/(70*VdList[1]*log(2)/ThalfList[1])
#pfna
0.91*1*median(df_imput$wpfna, na.rm=T)/(70*VdList[2]*log(2)/ThalfList[2])
#npfos
0.91*1*median(df_imput$wnpfos, na.rm=T)/(70*VdList[3]*log(2)/ThalfList[3])
#brpfos
0.91*1*median(df_imput$wbrpfos, na.rm=T)/(70*VdList[4]*log(2)/ThalfList[4])
#pfhxs
0.91*1*median(df_imput$wpfhxs, na.rm=T)/(70*VdList[5]*log(2)/ThalfList[5])


#1990 TK
for (i in 1:5)
{
  print(ChemList[i])
  print(table(df_imput[paste0(ChemList[i],"_flags")]))
  for (j  in 1:dim(df_imput)[1]){
    #print(j)
    if (df_imput[j,paste0(ChemList[i],"_flags")]==""){
      #print('detect')
      C_water<-df_imput[j,paste0('w',ChemList[i],'_imput')]
      df_imput[j,paste0('sm',ChemList[i])]<-AFList[i]*df_imput$tapnum[j]*lpcup*C_water/(df_imput$wt90f[j]*kgplb*VdList[i]*log(2)/ThalfList[i])
      if (!is.na(df_imput[j,paste0('s',ChemList[i],'b')])) {
        df_imput[j,paste0('r_',ChemList[i])]<-df_imput[j,paste0('sm',ChemList[i])]/df_imput[j,paste0('s',ChemList[i],'b')]
      }}
    else {
      #print('non detect')
      df_imput[j,paste0('sm',ChemList[i])]<-NA
      df_imput[j,paste0('r_',ChemList[i])]<-NA
    }
  }
  print(summary(df_imput[paste0('r_',ChemList[i])]))
}

#visualize RSC (relative source contribution)
df2<-df_imput[,c('id',paste0('r_',ChemList))]
sapply(df2,IQR, na.rm=T)
sapply(df2,cv, na.rm=T)
quantile(df2$r_pfoa,0.25, na.rm=T)
quantile(df2$r_pfoa,0.5, na.rm=T)
quantile(df2$r_pfoa,0.75, na.rm=T)
df2l<-melt(df2,"id")
png("water_rsc_012318.png",height = 3.2, width = 3.4, pointsize = 9, unit='in', res=500) 
p<-ggplot(df2l, aes(x=variable, y=value)) +
  geom_boxplot(outlier.shape=NA,fill='skyblue1') +
  geom_jitter(aes(variable,value),
              position=position_jitter(width=0.1,height=0),
              alpha=0.6,
              colour='grey30',
              size=0.7,
              show_guide=FALSE) +
  geom_hline(yintercept = 0.2, colour='navyblue', linetype='dashed')+
    scale_y_continuous(labels = scales::percent)+
    coord_cartesian(ylim=c(0,1.5))+
    labs(x='',y='')+
    scale_x_discrete(labels=c('PFOA','PFNA',#'PFDA',
                              'nPFOS','brPFOS','PFHxS'))+
    theme(panel.background = element_blank(),axis.line = element_line(colour = "black"))
p
dev.off()

df3<-df_imput[,c('id',paste0('r_',ChemList),'res_yr')]
df3l<-melt(df3,c("id",'res_yr'))
df3l$variable_n[df3l$variable == "r_pfoa"] <- 1
df3l$variable_n[df3l$variable == "r_pfna"] <- 2
df3l$variable_n[df3l$variable == "r_npfos"] <- 3
df3l$variable_n[df3l$variable == "r_brpfos"] <- 4
df3l$variable_n[df3l$variable == "r_pfhxs"] <- 5

df3l$scat_adj[df3l$res_yr == "<2"] <- -0.3
df3l$scat_adj[df3l$res_yr == "2-4"] <- -0.1
df3l$scat_adj[df3l$res_yr == "4-14"] <- 0.1
df3l$scat_adj[df3l$res_yr == ">14"] <- 0.3
df3l$scat_adj[df3l$res_yr == "4-14"&& df3l$variable=='r_pfna'] <- 0
png("water_rsc_by_res_yr_012318.png",height = 3.4, width = 3.4, pointsize = 9, unit='in', res=500) 
p<-ggplot(df3l, aes(x=variable, y=value, fill=res_yr)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(variable_n+scat_adj,value),
              position=position_jitter(width=0.05,height=0),
              alpha=0.5,
              colour='grey30',
              size=0.4,
              show_guide=FALSE) +
  geom_hline(yintercept = 0.2, colour='navyblue', linetype='dashed')+
  scale_y_continuous(labels = scales::percent)+
  coord_cartesian(ylim=c(0,1.5))+
  labs(x='',y='')+
  scale_x_discrete(labels=c('PFOA','PFNA',#'PFDA',
                            'nPFOS','brPFOS','PFHxS'))+
  scale_fill_brewer(palette='Blues',name="Years in current residence",
                    guide = guide_legend(
                      direction = "horizontal",
                      title.position = "top"
                    )) +
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
        legend.position="top")
p
dev.off()

