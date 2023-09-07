library(tidyverse)
library(ggpubr)
library(data.table)
library(readr)
library(readxl)
library(formattable)
library(ComplexHeatmap)
library(OlinkAnalyze)
library(gtsummary)
library(flextable)
library(ggplotify)
library(pheatmap)
library(plotly)
library(gapminder)
library(ggfortify)
library(reactable)
library(sandwich)
library(ggcorrplot)
library(corrplot)


##------------------------------------------- import data
BioMe_proteome_PFAS_wide <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_wide.txt")
protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")
BioMe_proteome_PFAS_long <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_long.txt")



##------------------------------------------- data processing
# ## restrict to protein in inflammation panel
# BioMe_proteome_PFAS_allmation<-  BioMe_proteome_PFAS_long %>% 
#                                     filter(Panel == "Inflammation")

## protein in inflammation panel
protein_in_allpanels<- protein_in_panel$protein




## date of blood draw
BioMe_proteome_PFAS_wide$date_enrl <- rep(NA_real_, nrow(BioMe_proteome_PFAS_wide))
BioMe_proteome_PFAS_wide$month_yr_enrl<- as.character(BioMe_proteome_PFAS_wide$month_yr_enrl)


for(i in 1:nrow(BioMe_proteome_PFAS_wide)){
  
  x <- anytime::anydate(paste((strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][2]), " 1,", 2000 + as.numeric(strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][1])))
  mydates <- as.Date(c("2011-01-01"))
  BioMe_proteome_PFAS_wide$date_enrl[i] <- as.numeric((x - mydates[1])/365 )
  
}

BioMe_proteome_PFAS_wide$c_date_enrl <- ifelse(BioMe_proteome_PFAS_wide$date_enrl > 0, 1,0)

## race
BioMe_proteome_PFAS_wide$self_reported_race<- factor(BioMe_proteome_PFAS_wide$self_reported_race,
                                                     levels = c("Hispanic", "African American", "European American"))


## sex
BioMe_proteome_PFAS_wide$sex<- factor(BioMe_proteome_PFAS_wide$sex,
                                      levels = c("Male", "Female"))

## smoking
BioMe_proteome_PFAS_wide$smoking_at_enrollment<- factor(BioMe_proteome_PFAS_wide$smoking_at_enrollment,
                                                        levels = c("No", "Yes"))

## status
BioMe_proteome_PFAS_wide$status<- factor(ifelse(BioMe_proteome_PFAS_wide$status=="case", 1, 0),
                                         levels = c(0,1))


## analysis dataset
# analysis_data<- BioMe_proteome_PFAS_wide %>% 
#                 select(starts_with("OID"), self_reported_race, sex, age_at_enrollment, smoking_at_enrollment, bmi_at_enrollment, date_enrl, PFDA_Aug21_bi)

###### before dummify, we need to handle missing carefully in categorical variables
# analysis_data<- as.data.frame(DataExplorer::dummify(analysis_data))


# # correlation plot
# #--- corrplot
# corr_fun<- function(data){
#   corrplot(data, 
#            method="color" ,
#            col = colorRampPalette(c("steelblue", "white", "darkred"))(100),cl.lim=c(0,1),
#            type="full",
#            # order="hclust" ,
#            # addrect = 2,
#            tl.pos = 'lt',
#            tl.srt=30,
#            tl.col = "black",
#   )
# }
# 
# 
# 
# cor_data<- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_allpanels)), use = 'complete.obs')
# 
# cor_plot<- corr_fun(cor_data) 

#----------------------- Binary PFDA
##---------------------- Uadjusted
PFDA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFDA_Aug21_bi")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
            Estimate = coef(s_lm)
            , "Robust SE" = std.err
            , z = (coef(s_lm)/std.err)
            , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
          
  
  PFDA_all_lm <- rbind(PFDA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}
 
 

PFDA_all_lm <- PFDA_all_lm[-1,]
PFDA_all_lm$z.value <- as.numeric(PFDA_all_lm$z.value)
PFDA_all_lm$p.value <- as.numeric(PFDA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFDA_all_lm$p.value), lambda=0)
PFDA_all_lm$q.value <-  q$qvalues


PFDA_all_unlm_results<- PFDA_all_lm %>% 
                           left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
                           distinct()


write.csv(PFDA_all_unlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_unlm.csv",
          row.names = F)



##---------------------- Adjusted
PFDA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFDA_Aug21_bi + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment+ c_date_enrl + smoking_at_enrollment")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFDA_all_lm <- rbind(PFDA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFDA_all_lm <- PFDA_all_lm[-1,]
PFDA_all_lm$z.value <- as.numeric(PFDA_all_lm$z.value)
PFDA_all_lm$p.value <- as.numeric(PFDA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFDA_all_lm$p.value), lambda=0)
PFDA_all_lm$q.value <-  q$qvalues


PFDA_all_adlm_results<- PFDA_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFDA_all_adlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_adlm.csv",
          row.names = F)


#----------------------- Continuous PFDA
##---------------------- Uadjusted
PFDA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFDA_Aug21")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFDA_all_lm <- rbind(PFDA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFDA_all_lm <- PFDA_all_lm[-1,]
PFDA_all_lm$z.value <- as.numeric(PFDA_all_lm$z.value)
PFDA_all_lm$p.value <- as.numeric(PFDA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFDA_all_lm$p.value), lambda=0)
PFDA_all_lm$q.value <-  q$qvalues


PFDA_all_unlm_results<- PFDA_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFDA_all_unlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_unlm_cont.csv",
          row.names = F)



##---------------------- Adjusted
PFDA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFDA_Aug21 + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment+ c_date_enrl + smoking_at_enrollment")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFDA_all_lm <- rbind(PFDA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFDA_all_lm <- PFDA_all_lm[-1,]
PFDA_all_lm$z.value <- as.numeric(PFDA_all_lm$z.value)
PFDA_all_lm$p.value <- as.numeric(PFDA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFDA_all_lm$p.value), lambda=0)
PFDA_all_lm$q.value <-  q$qvalues


PFDA_all_adlm_results<- PFDA_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFDA_all_adlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_adlm_cont.csv",
          row.names = F)

#----------------------- Continuous PFDA - tertile
##---------------------- Uadjusted
PFDA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFDA_Aug21_q")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFDA_all_lm <- rbind(PFDA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFDA_all_lm <- PFDA_all_lm[-1,]
PFDA_all_lm$z.value <- as.numeric(PFDA_all_lm$z.value)
PFDA_all_lm$p.value <- as.numeric(PFDA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFDA_all_lm$p.value), lambda=0)
PFDA_all_lm$q.value <-  q$qvalues


PFDA_all_unlm_results<- PFDA_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFDA_all_unlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_unlm_q.csv",
          row.names = F)



##---------------------- Adjusted
PFDA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFDA_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment+ c_date_enrl + smoking_at_enrollment")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFDA_all_lm <- rbind(PFDA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFDA_all_lm <- PFDA_all_lm[-1,]
PFDA_all_lm$z.value <- as.numeric(PFDA_all_lm$z.value)
PFDA_all_lm$p.value <- as.numeric(PFDA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFDA_all_lm$p.value), lambda=0)
PFDA_all_lm$q.value <-  q$qvalues


PFDA_all_adlm_results<- PFDA_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFDA_all_adlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_adlm_q.csv",
          row.names = F)


#----------------------- Continuous PFOA - tertile
##---------------------- Uadjusted
PFOA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFOA_Aug21_q")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFOA_all_lm <- rbind(PFOA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFOA_all_lm <- PFOA_all_lm[-1,]
PFOA_all_lm$z.value <- as.numeric(PFOA_all_lm$z.value)
PFOA_all_lm$p.value <- as.numeric(PFOA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFOA_all_lm$p.value), lambda=0)
PFOA_all_lm$q.value <-  q$qvalues


PFOA_all_unlm_results<- PFOA_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFOA_all_unlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFOA_all_unlm_q.csv",
          row.names = F)



##---------------------- Adjusted
PFOA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFOA_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment+ c_date_enrl + smoking_at_enrollment")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFOA_all_lm <- rbind(PFOA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFOA_all_lm <- PFOA_all_lm[-1,]
PFOA_all_lm$z.value <- as.numeric(PFOA_all_lm$z.value)
PFOA_all_lm$p.value <- as.numeric(PFOA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFOA_all_lm$p.value), lambda=0)
PFOA_all_lm$q.value <-  q$qvalues


PFOA_all_adlm_results<- PFOA_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFOA_all_adlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFOA_all_adlm_q.csv",
          row.names = F)


#----------------------- Continuous PFOS - tertile
##---------------------- Uadjusted
PFOS_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFOS_Aug21_q")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFOS_all_lm <- rbind(PFOS_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFOS_all_lm <- PFOS_all_lm[-1,]
PFOS_all_lm$z.value <- as.numeric(PFOS_all_lm$z.value)
PFOS_all_lm$p.value <- as.numeric(PFOS_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFOS_all_lm$p.value), lambda=0)
PFOS_all_lm$q.value <-  q$qvalues


PFOS_all_unlm_results<- PFOS_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFOS_all_unlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFOS_all_unlm_q.csv",
          row.names = F)



##---------------------- Adjusted
PFOS_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFOS_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment+ c_date_enrl + smoking_at_enrollment")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFOS_all_lm <- rbind(PFOS_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFOS_all_lm <- PFOS_all_lm[-1,]
PFOS_all_lm$z.value <- as.numeric(PFOS_all_lm$z.value)
PFOS_all_lm$p.value <- as.numeric(PFOS_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFOS_all_lm$p.value), lambda=0)
PFOS_all_lm$q.value <-  q$qvalues


PFOS_all_adlm_results<- PFOS_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFOS_all_adlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFOS_all_adlm_q.csv",
          row.names = F)



#----------------------- Continuous PFHpA - tertile
##---------------------- Uadjusted
PFHpA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFHpA_Aug21_q")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFHpA_all_lm <- rbind(PFHpA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFHpA_all_lm <- PFHpA_all_lm[-1,]
PFHpA_all_lm$z.value <- as.numeric(PFHpA_all_lm$z.value)
PFHpA_all_lm$p.value <- as.numeric(PFHpA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFHpA_all_lm$p.value), lambda=0)
PFHpA_all_lm$q.value <-  q$qvalues


PFHpA_all_unlm_results<- PFHpA_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFHpA_all_unlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHpA_all_unlm_q.csv",
          row.names = F)



##---------------------- Adjusted
PFHpA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFHpA_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment+ c_date_enrl + smoking_at_enrollment")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFHpA_all_lm <- rbind(PFHpA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFHpA_all_lm <- PFHpA_all_lm[-1,]
PFHpA_all_lm$z.value <- as.numeric(PFHpA_all_lm$z.value)
PFHpA_all_lm$p.value <- as.numeric(PFHpA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFHpA_all_lm$p.value), lambda=0)
PFHpA_all_lm$q.value <-  q$qvalues


PFHpA_all_adlm_results<- PFHpA_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFHpA_all_adlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHpA_all_adlm_q.csv",
          row.names = F)



#----------------------- Continuous PFHxS - tertile
##---------------------- Uadjusted
PFHxS_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFHxS_Aug21_q")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFHxS_all_lm <- rbind(PFHxS_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFHxS_all_lm <- PFHxS_all_lm[-1,]
PFHxS_all_lm$z.value <- as.numeric(PFHxS_all_lm$z.value)
PFHxS_all_lm$p.value <- as.numeric(PFHxS_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFHxS_all_lm$p.value), lambda=0)
PFHxS_all_lm$q.value <-  q$qvalues


PFHxS_all_unlm_results<- PFHxS_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFHxS_all_unlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHxS_all_unlm_q.csv",
          row.names = F)



##---------------------- Adjusted
PFHxS_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFHxS_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment+ c_date_enrl + smoking_at_enrollment")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFHxS_all_lm <- rbind(PFHxS_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFHxS_all_lm <- PFHxS_all_lm[-1,]
PFHxS_all_lm$z.value <- as.numeric(PFHxS_all_lm$z.value)
PFHxS_all_lm$p.value <- as.numeric(PFHxS_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFHxS_all_lm$p.value), lambda=0)
PFHxS_all_lm$q.value <-  q$qvalues


PFHxS_all_adlm_results<- PFHxS_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFHxS_all_adlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHxS_all_adlm_q.csv",
          row.names = F)




#----------------------- Continuous PFNA - tertile
##---------------------- Uadjusted
PFNA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFNA_Aug21_q")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFNA_all_lm <- rbind(PFNA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFNA_all_lm <- PFNA_all_lm[-1,]
PFNA_all_lm$z.value <- as.numeric(PFNA_all_lm$z.value)
PFNA_all_lm$p.value <- as.numeric(PFNA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFNA_all_lm$p.value), lambda=0)
PFNA_all_lm$q.value <-  q$qvalues


PFNA_all_unlm_results<- PFNA_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFNA_all_unlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFNA_all_unlm_q.csv",
          row.names = F)



##---------------------- Adjusted
PFNA_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFNA_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment+ c_date_enrl + smoking_at_enrollment")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFNA_all_lm <- rbind(PFNA_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFNA_all_lm <- PFNA_all_lm[-1,]
PFNA_all_lm$z.value <- as.numeric(PFNA_all_lm$z.value)
PFNA_all_lm$p.value <- as.numeric(PFNA_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFNA_all_lm$p.value), lambda=0)
PFNA_all_lm$q.value <-  q$qvalues


PFNA_all_adlm_results<- PFNA_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFNA_all_adlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFNA_all_adlm_q.csv",
          row.names = F)






#----------------------- Continuous PFHpS - tertile
##---------------------- Uadjusted
PFHpS_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFHpS_Aug21_q")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFHpS_all_lm <- rbind(PFHpS_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFHpS_all_lm <- PFHpS_all_lm[-1,]
PFHpS_all_lm$z.value <- as.numeric(PFHpS_all_lm$z.value)
PFHpS_all_lm$p.value <- as.numeric(PFHpS_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFHpS_all_lm$p.value), lambda=0)
PFHpS_all_lm$q.value <-  q$qvalues


PFHpS_all_unlm_results<- PFHpS_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFHpS_all_unlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHpS_all_unlm_q.csv",
          row.names = F)



##---------------------- Adjusted
PFHpS_all_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein_in_allpanels)){
  
  s_lm <- (lm(as.formula(paste0(protein_in_allpanels[i], "~ PFHpS_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment+ c_date_enrl + smoking_at_enrollment")), 
              data = BioMe_proteome_PFAS_wide))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  PFHpS_all_lm <- rbind(PFHpS_all_lm, c(protein_in_allpanels[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
}



PFHpS_all_lm <- PFHpS_all_lm[-1,]
PFHpS_all_lm$z.value <- as.numeric(PFHpS_all_lm$z.value)
PFHpS_all_lm$p.value <- as.numeric(PFHpS_all_lm$p.value)


q <-qvalue::qvalue(as.numeric(PFHpS_all_lm$p.value), lambda=0)
PFHpS_all_lm$q.value <-  q$qvalues


PFHpS_all_adlm_results<- PFHpS_all_lm %>% 
  left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()


write.csv(PFHpS_all_adlm_results,
          "~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHpS_all_adlm_q.csv",
          row.names = F)
















