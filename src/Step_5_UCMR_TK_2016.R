#purpose of this script: compare county-specific PFAS concentrations in 1989/1990 and 2013-2015
#focus on PFOS, PFOA, PFNA, PFHxS
library(dplyr)
library(visreg)
library(mgcv)
library(ggplot2)
library(reshape2)
ucmr3_by_county<-read.csv('ucmrByCounty_final_030218_orig.csv')
PFAS_2015<-cbind(ucmr3_by_county[,c('COUNTY','countyname','state')],
                 ucmr3_by_county[,grep('PFBS',colnames(ucmr3_by_county))],
                 ucmr3_by_county[,grep('PFHpA',colnames(ucmr3_by_county))],
                 ucmr3_by_county[,grep('PFHxS',colnames(ucmr3_by_county))],
                 ucmr3_by_county[,grep('PFNA',colnames(ucmr3_by_county))],
                 ucmr3_by_county[,grep('PFOA',colnames(ucmr3_by_county))],
                 ucmr3_by_county[,grep('PFOS',colnames(ucmr3_by_county))])
PFAS_2015<-cbind(PFAS_2015[,c('COUNTY','countyname','state')],
                 PFAS_2015[,grep('AnalyticalResult*',colnames(PFAS_2015))])
colnames(PFAS_2015)<-gsub('AnalyticalResult','',colnames(PFAS_2015))
colnames(PFAS_2015)<-gsub('Value.','',colnames(PFAS_2015))
colnames(PFAS_2015)<-gsub('sSign.','Freq.',colnames(PFAS_2015))
PFAS_2015<-PFAS_2015[!is.na(PFAS_2015$count.PFBS),]
#change unit from ppb in UCMR3 to ppt for easier comparison with NHS
PFAS_2015[,c(5,8,11,14,17,20)]<-1000*PFAS_2015[,c(5,8,11,14,17,20)]
tail(PFAS_2015)
colnames(PFAS_2015)[4:21]<-paste0(colnames(PFAS_2015)[4:21],'.2015')
temp1<-left_join(df_imput[!df_imput$stct%in%PFAS_2015$COUNTY,], PFAS_2015,
            by=c('state90'='state','county90'='countyname'))
df_imput_ucmr<-left_join(df_imput,PFAS_2015,by=c('stct'='COUNTY'))
key<-intersect(colnames(temp1),colnames(df_imput_ucmr))
temp1<-temp1[,key]
df_imput_ucmr<-df_imput_ucmr[,key]
df_imput_ucmr[!df_imput_ucmr$stct%in%PFAS_2015$COUNTY,]<-temp1
df_imput_ucmr<-unique(df_imput_ucmr)
#water_master_clean[water_master_clean$County90!=water_master_clean$countyname,]
colnames(df_imput_ucmr)

#tabulate detection in 2015
key<-Reduce(union, list(which(df_imput_ucmr$Freq.mean.PFBS.2015>0),
      which(df_imput_ucmr$Freq.mean.PFHpA.2015>0),
      which(df_imput_ucmr$Freq.mean.PFHxS.2015>0),
      which(df_imput_ucmr$Freq.mean.PFNA.2015>0),
      which(df_imput_ucmr$Freq.mean.PFOS.2015>0),
      which(df_imput_ucmr$Freq.mean.PFOA.2015>0)))
length(key)/dim(df_imput_ucmr)[1]
#46% homes are from UCMR3 positive counties
df_imput_ucmr_detect<-df_imput_ucmr[key,]
length(unique(df_imput_ucmr_detect$state90))
length(unique(df_imput_ucmr_detect$stct))
length(unique(df_imput_ucmr$state90))
length(unique(df_imput_ucmr$stct))

#compare 1990 with 2015
cor.test(df_imput_ucmr$Freq.mean.PFHxS.2015,df_imput_ucmr$wpfhxs_imput,
         use='complete.obs',method='spearman') #n.s.
cor.test(df_imput_ucmr$max.PFHxS.2015,(df_imput_ucmr$wpfhxs_imput),
         use='complete.obs',method='spearman') #n.s.

cor.test(df_imput_ucmr$Freq.mean.PFNA.2015,(df_imput_ucmr$wpfna_imput),
         use='complete.obs',method='spearman') #n.s.
cor.test(df_imput_ucmr$max.PFNA.2015,(df_imput_ucmr$wpfna_imput),
         use='complete.obs',method='spearman') #n.s.

cor.test(df_imput_ucmr$Freq.mean.PFOA.2015,(df_imput_ucmr$wpfoa_imput),
         use='complete.obs',method='spearman') #n.s.
cor.test(df_imput_ucmr$max.PFOA.2015,(df_imput_ucmr$wpfoa_imput), 
         use='complete.obs',method='spearman') #n.s. 
#sig
cor.test(df_imput_ucmr$Freq.mean.PFOS.2015,(df_imput_ucmr$wtotpfos),
         use='complete.obs',method='spearman') #rho=0.20, p=0.006
cor.test(df_imput_ucmr$max.PFOS.2015,(df_imput_ucmr$wtotpfos),
         use='complete.obs',method='spearman') #rho=0.19, p=0.008

#evaluate within county variability
df_imput_ucmr<-df_imput_ucmr[order(df_imput_ucmr$stct),]
key<-unique(df_imput_ucmr$stct)[table(df_imput_ucmr$stct)>1]
county_sub<-df_imput_ucmr[df_imput_ucmr$stct%in%key,c('stct','wpfhxs_imput','wpfna_imput',
'wtotpfos','wpfoa_imput')]
county_sub[is.na(county_sub)]<-0
ag <- aggregate(. ~ stct, county_sub, function(x) c(max= max(x),cv = sd(x)/mean(x)))
sapply(ag,summary)
df_imput_ucmr[which(df_imput_ucmr$stct=='25009'),] #PFOS range from non detect to 67 ng/L
#CV are moderate, therefore the within county variability is unlikely to explain the temp diff


#TK 2016
ChemList<-c('pfoa','pfna',#'pfda',
            'npfos','brpfos','pfhxs')
ChemList2015<-c('PFOA','PFNA','PFOS','PFOS','PFHxS')
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

for (i in 1:5)
{
  print(ChemList[i])
  print(table(df_imput_ucmr[paste0(ChemList[i],"_flags")]))
  for (j  in 1:dim(df_imput_ucmr)[1]){
    #print(j)
    if (df_imput_ucmr[j,paste0(ChemList[i],"_flags")]==""){
      #print('detect')
      C_water<-df_imput_ucmr[j,paste0('max.',ChemList2015[i],'.2015')]
      df_imput_ucmr[j,paste0('sm2015',ChemList[i])]<-AFList[i]*df_imput_ucmr$tapnum[j]*lpcup*C_water/(df_imput_ucmr$wt90f[j]*kgplb*VdList[i]*log(2)/ThalfList[i])
    }
    else {
      #print('non detect')
      df_imput_ucmr[j,paste0('sm2015',ChemList[i])]<-NA
    }
  }
}

#assume linear PFOS: branched PFOS = 70:30
df_imput_ucmr$sm2015brpfos<-df_imput_ucmr$sm2015brpfos*0.3
df_imput_ucmr$sm2015npfos<-df_imput_ucmr$sm2015npfos*0.7

summary(df_imput_ucmr$sm2015brpfos)
summary(df_imput_ucmr$smbrpfos)
summary(df_imput_ucmr$sm2015npfos)
summary(df_imput_ucmr$smnpfos)
summary(df_imput_ucmr$sm2015pfoa)
summary(df_imput_ucmr$smpfoa)
summary(df_imput_ucmr$sm2015pfhxs)
summary(df_imput_ucmr$smpfhxs)
summary(df_imput_ucmr$sm2015pfna)
summary(df_imput_ucmr$smpfna)

#the 8 locations in 2016
idList<-c('108067','134102','126013','112072','246454','334458','253856','131456')
df_imput_ucmr[df_imput_ucmr$id %in% idList,c('id','wpfoa_imput','mean.PFOA.2015',
                                             'wpfna_imput','mean.PFOA.2015',
                                             'wnpfos_imput','mean.PFOS.2015',
                                             "wbrpfos_imput",
                                             'wpfhxs_imput','mean.PFHxS.2015')]
#visualization
#subset
compare_1990_2015<-df_imput_ucmr[,c('id','smpfoa','smpfna','smnpfos','smbrpfos','smpfhxs',
                                    'sm2015pfoa','sm2015pfna','sm2015npfos','sm2015brpfos','sm2015pfhxs')]
#wide to long
compare_1990_2015_l<-melt(compare_1990_2015, id.vars=c("id"))
compare_1990_2015_l$chem<-gsub('sm','',compare_1990_2015_l$variable)
compare_1990_2015_l$year<-'1990'
compare_1990_2015_l$year[grep('2015',compare_1990_2015_l$chem)]<-'2015'
compare_1990_2015_l$chem<-gsub('2015','',compare_1990_2015_l$chem)
compare_1990_2015_l$chem<-factor(compare_1990_2015_l$chem,levels = c('pfoa','pfna','npfos','brpfos','pfhxs'))
#test plot
ggplot(aes(y = value, x = year, fill = chem), data = compare_1990_2015_l) +
  geom_boxplot()+ylim(c(0,30))+ggtitle("orig")

compare_1990_2015_l$scat_adj[compare_1990_2015_l$chem == "pfoa"] <- -0.3
compare_1990_2015_l$scat_adj[compare_1990_2015_l$chem == "pfna"] <- -0.15
compare_1990_2015_l$scat_adj[compare_1990_2015_l$chem == "npfos"] <- 0
compare_1990_2015_l$scat_adj[compare_1990_2015_l$chem == "brpfos"] <- 0.15
compare_1990_2015_l$scat_adj[compare_1990_2015_l$chem == "pfhxs"] <- 0.3

compare_1990_2015_l$variable_n[compare_1990_2015_l$year == "1990"] <- 1
compare_1990_2015_l$variable_n[compare_1990_2015_l$year == "2015"] <- 2

png("sm_1990_2015_030118_orig.png",height = 3.4, width = 3.4, pointsize = 9, unit='in', res=500) 
p<-ggplot(compare_1990_2015_l, aes(x=year, y=value, fill=chem)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(variable_n+scat_adj,value),
              position=position_jitter(width=0.05,height=0),
              alpha=0.3,
              colour='grey30',
              size=0.4,
              show.legend=FALSE) +
  coord_cartesian(ylim=c(0,30))+
  labs(x='',y='Modeled serum PFAS (ng/mL)')+
  scale_x_discrete(labels=c('1989/1990','2013-2015'))+
  scale_fill_brewer(palette='Pastel2',name="",
                    labels = c("PFOA", "PFNA","nPFOS","brPFOS","PFHxS"),
                    guide = guide_legend(
                      direction = "horizontal",
                      title.position = "top"
                    )) +
  theme(panel.background = element_blank(),axis.line = element_line(colour = "black"),
        legend.position="top")
p
dev.off()
