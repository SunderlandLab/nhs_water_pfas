library(dplyr)
library(qwraps2)
library(reshape2)
library(ggplot2)
options(qwraps2_markup = "markdown")
#reference code from https://cran.r-project.org/web/packages/qwraps2/vignettes/summary-statistics.html
source('main.R')

#generating table1, characteristics of study participants, n=225
df_imput$tapnum[is.na(df_imput$tapnum)]<-df_imput$h2o90dn[is.na(df_imput$tapnum)]
df_imput$tapnumc<-cut(df_imput$tapnum, breaks = c(-0.1,2,5,9,17.5))
df_imput$tapnumc<-factor(df_imput$tapnumc)
df_imput$bmi90num<-20
df_imput$bmi90num[df_imput$bmi90==2]<-20.5
df_imput$bmi90num[df_imput$bmi90==3]<-21.5
df_imput$bmi90num[df_imput$bmi90==4]<-22.5
df_imput$bmi90num[df_imput$bmi90==5]<-23.5
df_imput$bmi90num[df_imput$bmi90==6]<-24.5
df_imput$bmi90num[df_imput$bmi90==7]<-26
df_imput$bmi90num[df_imput$bmi90==8]<-28
df_imput$bmi90num[df_imput$bmi90==9]<-29.5
df_imput$bmi90num[df_imput$bmi90==10]<-31
df_imput$bmi90num[df_imput$bmi90==11]<-33.5
df_imput$bmi90num[df_imput$bmi90==12]<-37.5
df_imput$bmi90num[df_imput$bmi90==13]<-40
df_imput$bmi90num[df_imput$bmi90==14]<-NA
#$label 1.bmi < 20.0;\
# 2.bmi 20.0-<21.0;\
# 3.bmi 21.0-<22.0;\
# 4.bmi 22.0-<23.0;\
# 5.bmi 23.0-<24.0;\
# 6.bmi 24.0-<25.0;\
# 7.bmi 25.0-<27.0;\
# 8.bmi 27.0-<29.0;\
# 9.bmi 29.0-<30.0;\
# 10.bmi 30.0-<32.0;\
# 11.bmi 32.0-<35.0;\
# 12.bmi 35.0-<40.0;\
# 13.bmi 40.0 +;\
# 14.missing
our_summary3<-
  list("age" = 
         list("mean (SD)" = ~ qwraps2::mean_sd(age, na_rm=TRUE,digits = getOption("qwraps2_frmt_digits", 1),
                                                     show_n = "never")),
       "white" = 
         list("White" = ~ qwraps2::n_perc(white==1,digits = getOption("qwraps2_frmt_digits", 0))),
       "bmi (kg/m2)" = 
         list("mean (SD)" = ~ qwraps2::mean_sd(bmi90num, na_rm=TRUE,digits = getOption("qwraps2_frmt_digits", 1),
                                               show_n = "never")),
       "weight (lb)" = 
         list("mean (SD)" = ~ qwraps2::mean_sd(wt90f, na_rm=TRUE, digits = getOption("qwraps2_frmt_digits", 1),
                                               show_n = "never")),
       "parity" = 
         list("No birth" = ~ qwraps2::n_perc(npar90==0,digits = getOption("qwraps2_frmt_digits", 0)),
              "1-3 birth" = ~ qwraps2::n_perc(npar90>0 & npar90<=3,digits = getOption("qwraps2_frmt_digits", 0)),
              "3+ births"= ~ qwraps2::n_perc(npar90>3,digits = getOption("qwraps2_frmt_digits", 0))),
       "breastfeeding (months)" = 
         list("Never breastfed" = ~ qwraps2::n_perc(brfdn==0, digits = getOption("qwraps2_frmt_digits", 0)),
              "< 12months" = ~ qwraps2::n_perc(brfdn<12 & brfdn>0, digits = getOption("qwraps2_frmt_digits", 0)),
              ">= 12months" = ~ qwraps2::n_perc(brfdn>=12, digits = getOption("qwraps2_frmt_digits", 0))),
       "menstruration status" = 
         list("Premenopause" = ~ qwraps2::n_perc(men90==1,digits = getOption("qwraps2_frmt_digits", 0)),
              "Postmenopause" = ~ qwraps2::n_perc(men90==2,digits = getOption("qwraps2_frmt_digits", 0))),
       "seafood consumption (servings/day)" = 
         list("mean (SD)" = ~ qwraps2::mean_sd(seafood, na_rm=TRUE, digits = getOption("qwraps2_frmt_digits", 1),
                                               show_n = "never")),
       "popcorn consumption (servings/day)" = 
         list("mean (SD)" = ~ qwraps2::mean_sd(popc90dn, na_rm=TRUE, digits = getOption("qwraps2_frmt_digits", 1),
                                               show_n = "never")),
       "years residing in current location" = 
         list("<2" = ~ qwraps2::n_perc(res_yr=='<2',digits = getOption("qwraps2_frmt_digits", 0)),
              "2-4" = ~ qwraps2::n_perc(res_yr=='2-4',digits = getOption("qwraps2_frmt_digits", 0)),
              "4-14"= ~ qwraps2::n_perc(res_yr=='4-14',digits = getOption("qwraps2_frmt_digits", 0)),
              ">14" = ~ qwraps2::n_perc(res_yr=='>14',digits = getOption("qwraps2_frmt_digits", 0)))) 

nhstable1<-summary_table(dplyr::group_by(df_imput,tapnumc), our_summary3)
write.csv(nhstable1,'nhstable1_225_041218.csv')

#generating table 2, water PFASs in 225 archived samples
serumid<-df_imput$id[!is.na(df_imput$snpfos)]
water_master_clean$inserum<-F
water_master_clean$inserum[water_master_clean$id%in%serumid]<-T
water_master_clean$inserum<-factor(water_master_clean$inserum)
water_master_clean<-water_master_clean[water_master_clean$id%in%df_imput$id,]

args(summary_table)
our_summary1<-
  list("PFPeA" =
         list("% detect" = ~ qwraps2::n_perc0(PFPeA_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFPeA, na_rm=TRUE,
                                                     show_n = "never"),
              "max"=~round(max(PFPeA, na.rm=T),2)),
       "PFHpA" =
         list("% detect" = ~ qwraps2::n_perc0(PFHpA_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFHpA, na_rm = T,show_n = "never"),
              "max"=~round(max(PFHpA, na.rm=T),2)),
       "nPFOA" =
         list("% detect" = ~ qwraps2::n_perc0(PFOA_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFOA, na_rm = T,show_n = "never"),
              "max"=~round(max(PFOA, na.rm=T),2)),
       "brPFOA" =
         list("% detect" = ~ qwraps2::n_perc0(PFOA_br_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFOA_br, na_rm = T,show_n = "never"),
              "max"=~round(max(PFOA_br, na.rm=T),2)),
       "nPFNA" =
         list("% detect" = ~ qwraps2::n_perc0(PFNA_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFNA, na_rm = T,show_n = "never"),
              "max"=~round(max(PFNA, na.rm=T),2)),
       "brPFNA" =
         list("% detect" = ~ qwraps2::n_perc0(PFNA_br_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFNA_br, na_rm = T,show_n = "never"),
              "max"=~round(max(PFNA_br, na.rm=T),2)),
       "PFDA" =
         list("% detect" = ~ qwraps2::n_perc0(PFDA_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFDA, na_rm = T,show_n = "never"),
              "max"=~round(max(PFDA, na.rm=T),2)),
       "PFUnDA" =
         list("% detect" = ~ qwraps2::n_perc0(PFUA_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFUA, na_rm = T,show_n = "never"),
              "max"=~round(max(PFUA, na.rm=T),2)),
       "PFDoDA" =
         list("% detect" = ~ qwraps2::n_perc0(PFDoA_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFDoA, na_rm = T,show_n = "never"),
              "max"=~round(max(PFDoA, na.rm=T),2)),
       "PFBS" =
         list("% detect" = ~ qwraps2::n_perc0(PFBS_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFBS, na_rm = T,show_n = "never"),
              "max"=~round(max(PFBS, na.rm=T),2)),
       "nPFHxS" =
         list("% detect" = ~ qwraps2::n_perc0(PFHxS_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFHxS, na_rm = T,show_n = "never"),
              "max"=~round(max(PFHxS, na.rm=T),2)),
       "brPFHxS" =
         list("% detect" = ~ qwraps2::n_perc0(PFHxS_br_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFHxS_br, na_rm = T,show_n = "never"),
              "max"=~round(max(PFHxS_br, na.rm=T),2)),
       "nPFOS" =
         list("% detect" = ~ qwraps2::n_perc0(PFOS_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFOS, na_rm = T,show_n = "never"),
              "max"=~round(max(PFOS, na.rm=T),2)),
       "brPFOS" =
         list("% detect" = ~ qwraps2::n_perc0(PFOS_br_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFOS_br, na_rm = T,show_n = "never"),
              "max"=~round(max(PFOS_br, na.rm=T),2)),
       "PFDS" =
         list("% detect" = ~ qwraps2::n_perc0(PFDS_Flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(PFDS, na_rm = T,show_n = "never"),
              "max"=~round(max(PFDS, na.rm=T),2))

  )
#result<-summary_table(dplyr::group_by(water_master_clean, inserum), our_summary1)
result<-summary_table(water_master_clean,our_summary1)
result2<-cbind(result,rep(1:15,each=3))
result2<-as.data.frame(result2)
result2$timevar<-rownames(result2)
result3<-reshape(result2, idvar = "V2", timevar = "timevar", direction = "wide")

rownames(result3)<-c('PFPeA', 'PFHpA','nPFOA','brPFOA','nPFNA','brPFNA','PFDA','PFUnDA','PFDoDA',
                     'PFBS','nPFHxS','brPFHxS','nPFOS','brPFOS','PFDS')
write.csv(result3,'watertable2_225_combined_030818.csv')
table(water_master_clean$tPFOS>11)
table(water_master_clean$PFOA+water_master_clean$PFOA_br>14)
table(water_master_clean$PFNA+water_master_clean$PFNA_br>13)

#Generating table3, describe matched water and serum data (raw)
our_summary2<-
  list("spfoa" =
         list("% detect" = ~ qwraps2::n_perc0(spfoa>0.015),
              "median (IQR)" = ~ qwraps2::median_iqr(spfoa, na_rm=TRUE,
                                                     show_n = "never"),
              "max"=~round(max(spfoa, na.rm=T),2)),
       "spfna" =
         list("% detect" = ~ qwraps2::n_perc0(spfna>0.015),
              "median (IQR)" = ~ qwraps2::median_iqr(spfna, na_rm = T,show_n = "never"),
              "max"=~round(max(spfna, na.rm=T),2)),
       "snpfos" =
         list("% detect" = ~ qwraps2::n_perc0(snpfos>0.015),
              "median (IQR)" = ~ qwraps2::median_iqr(snpfos, na_rm = T,show_n = "never"),
              "max"=~round(max(snpfos, na.rm=T),2)),
       "sbrpfos" =
         list("% detect" = ~ qwraps2::n_perc0(sbrpfos>0.015),
              "median (IQR)" = ~ qwraps2::median_iqr(sbrpfos, na_rm = T,show_n = "never"),
              "max"=~round(max(sbrpfos, na.rm=T),2)),
       "spfhxs" =
         list("% detect" = ~ qwraps2::n_perc0(spfhxs>0.015),
              "median (IQR)" = ~ qwraps2::median_iqr(spfhxs, na_rm = T,show_n = "never"),
              "max"=~round(max(spfhxs, na.rm=T),2)),
       "wpfoa" =
         list("% detect" = ~ qwraps2::n_perc0(pfoa_flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(wpfoa, na_rm = T,show_n = "never"),
              "max"=~round(max(wpfoa, na.rm=T),2)),
       "wpfna" =
         list("% detect" = ~ qwraps2::n_perc0(pfna_flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(wpfna, na_rm = T,show_n = "never"),
              "max"=~round(max(wpfna, na.rm=T),2)),
       "wnpfos" =
         list("% detect" = ~ qwraps2::n_perc0(npfos_flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(wnpfos, na_rm = T,show_n = "never"),
              "max"=~round(max(wnpfos, na.rm=T),2)),
       "wbrpfos" =
         list("% detect" = ~ qwraps2::n_perc0(brpfos_flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(wbrpfos, na_rm = T,show_n = "never"),
              "max"=~round(max(wbrpfos, na.rm=T),2)),
       "wpfhxs" =
         list("% detect" = ~ qwraps2::n_perc0(pfhxs_flags==''),
              "median (IQR)" = ~ qwraps2::median_iqr(wpfhxs, na_rm = T,show_n = "never"),
              "max"=~round(max(wpfhxs, na.rm=T),2))
  )
RESULT<-summary_table(df[!is.na(df$snpfos),], our_summary2)
RESULT2<-cbind(RESULT,rep(1:10,each=3))
RESULT2<-as.data.frame(RESULT2)
RESULT2$timevar<-rownames(RESULT2)
RESULT3<-reshape(RESULT2, idvar = "V2", timevar = "timevar", direction = "wide")
rownames(RESULT3)<-c('serum PFOA', 'serum PFNA','serum nPFOS','serum brPFOS','serum PFHxS','water PFOA','water PFNA','water nPFOS','water brPFOS',
                     'water PFHxS')
write.csv(RESULT3,'table3.csv')

#visualize serum data
serumPFASs<-df[!is.na(df$snpfosb),c('id','spfoab','spfnab','snpfosb','sbrpfosb','spfhxsb')]
serumPFASs$spfos<-serumPFASs$sbrpfosb+serumPFASs$snpfosb
serumPFASs$snpfosb <- serumPFASs$sbrpfosb <- NULL
serumPFASs.l<-melt(serumPFASs,id="id")
#median from NHANES 1999-2000
data.NHANES<-data.frame(variable=c('spfoab','spfnab','spfos','spfhxsb'),value=c(4.6,0.5,27.7,1.7))
postscript('NHS_NHANES_median_031318.eps',horizontal = FALSE, onefile = FALSE, paper = "special",
           width = 5.5,height = 5.5, pointsize = 24)
ggplot(data = serumPFASs.l, aes(x = variable, y = value)) +
  geom_boxplot() +
  geom_point(data=data.NHANES, aes(x = variable, y = value), color='red', size=3)+
  labs(x='',y='')+
  scale_x_discrete(labels=c('PFOA','PFNA','PFOS','PFHxS'))+
  scale_y_continuous(name = 'serum PFASs (ng/mL)', trans = 'log10',breaks=c(1,10,100))+
  theme(text=element_text(size=18),panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()
postscript('NHS_noNHANES_median_031518.eps',horizontal = FALSE, onefile = FALSE, paper = "special",
           width = 5.5,height = 5.5, pointsize = 24)
ggplot(data = serumPFASs.l, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(x='',y='')+
  scale_x_discrete(labels=c('PFOA','PFNA','PFOS','PFHxS'))+
  scale_y_continuous(name = 'serum PFASs (ng/mL)', trans = 'log10',breaks=c(1,10,100))+
  theme(text=element_text(size=18),panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()
  # fun.y=mean calculates the mean of the y-axis variable (mpg) for each group.  
  # geom="point" indicates means should be plotted as points over top of the boxplots.
  #stat_summary(fun.y = mean, geom = "point", color = "red", size = 3)
