library(car)
library(readr)
library(lattice)
library(nlme)
library(ggplot2)
library(GGally)
library(foreign)
library(MASS)
library(lme4)
library(multcomp)
library(dplyr)
library(knitr)
library(xtable)
library(kableExtra)
library(glmnet)
library(corrplot)
library(ggpubr)
library(lmerTest)
library("merTools")
library(reshape2)
library(gplots)
library(tidyr)
library(blme)
library(grpreg)
library(gridExtra)
library(ggcorrplot)
library(BWQS)
library(data.table)
library(mice)
library(spatstat)
library(foreach)
library(doParallel)
library(iterators)
library(parallel)
library(anytime)
library(LUCIDus)


##------------------------------------------- q.value

BioMe_proteome_PFAS_wide <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_wide_imputed.txt")

proteome_vs_pfas_bwqs <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/scale/proteome_vs_pfas_bwqs.txt")
proteome_t2d_model <- fread("~/Projects/BioMe/proteome/input/pwas/proteome_vs_t2d_scale.txt")

## date of blood draw
BioMe_proteome_PFAS_wide$date_enrl <- rep(NA_real_, nrow(BioMe_proteome_PFAS_wide))
BioMe_proteome_PFAS_wide$month_yr_enrl<- as.character(BioMe_proteome_PFAS_wide$month_yr_enrl)


for(i in 1:nrow(BioMe_proteome_PFAS_wide)){
  
  x <- anytime::anydate(paste((strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][2]), " 1,", 2000 + as.numeric(strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][1])))
  mydates <- as.Date(c("2011-01-01"))
  BioMe_proteome_PFAS_wide$date_enrl[i] <- as.numeric((x - mydates[1])/365 )
  
}

BioMe_proteome_PFAS_wide$c_date_enrl <- ifelse(BioMe_proteome_PFAS_wide$date_enrl > 0, 1,0)


## preliminary screening: 64 out of 2612
proteome_pfas<- (proteome_vs_pfas_bwqs %>% filter(q.value < 0.1))$OlinkID
proteome_t2d<- (proteome_t2d_model %>% filter(q.value < 0.1))$OlinkID
common<- intersect(proteome_pfas, proteome_t2d)



## LUCID - find latent clusters
omics_selected<- BioMe_proteome_PFAS_wide %>% 
  dplyr::select(all_of(common))

exposure<- BioMe_proteome_PFAS_wide %>% 
  dplyr::select(ends_with("_q"))

outcome<- BioMe_proteome_PFAS_wide[, "td2_case_all"]

covariate<- as.data.frame(dummify(BioMe_proteome_PFAS_wide)) %>% 
  dplyr::select(self_reported_race.African.American, self_reported_race.European.American,
               age_at_enrollment,
               gender.Female,
               c_date_enrl)


set.seed(123)
fit_try1 = lucid(G = exposure, 
                 Z = omics_selected, 
                 Y = outcome, 
                 CoY = covariate,
                 family = "binary",
                 useY = FALSE,
                 K = 2:5)



summary_lucid(fit_try1)
plot_lucid(fit_try1,
           pos_link_color = "red",
           neg_link_color = "blue")


## creat exposure and omics profiles base on LUCID model
#========================================================================#
## 4. Distribution of zBMI score for each cluster predicted by LUCID    ##
#========================================================================#
# prediction of LUCID model
pred_fit_try2 = predict_lucid(model = fit_try1, 
                              lucid_model = "early",
                              G = exposure,
                              Z = omics_selected,
                              CoY = covariate)
# prediction on latent cluster
table(pred_fit_try2$pred.x)


Y_fit_try2 = as.data.frame(cbind(cluster = as.factor(pred_fit_try2$pred.x), BioMe_proteome_PFAS_wide[, "td2_case_all"]))
Y_fit_try2 = melt(Y_fit_try2, id.vars = "cluster")
ggplot(Y_fit_try2, aes(x = as.factor(cluster), 
                       fill = as.factor(value))) +
  geom_bar(stat = "count",
           position = "stack") +
  xlab("cluster") +
  ylab("z-bmi")












