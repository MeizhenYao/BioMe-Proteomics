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
BioMe_proteome_PFAS_wide <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_wide_imputed.txt")
protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")
BioMe_proteome_PFAS_long <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_long.txt")



##------------------------------------------- data processing
# ## restrict to protein in inflammation panel
# BioMe_proteome_PFAS_allmation<-  BioMe_proteome_PFAS_long %>% 
#                                     filter(Panel == "Inflammation")

## protein in inflammation panel
protein_in_allpanels<- protein_in_panel$OlinkID

## PFHpA binary formate
BioMe_proteome_PFAS_wide$PFHpA_Aug21_bi<- factor(BioMe_proteome_PFAS_wide$PFHpA_Aug21_bi,
                                          levels = c("Lower PFHpA", "Higher PFHpA"))


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
BioMe_proteome_PFAS_wide$sex<- factor(BioMe_proteome_PFAS_wide$gender,
                                      levels = c("Male", "Female"))

## smoking
BioMe_proteome_PFAS_wide$smoking_at_enrollment<- factor(BioMe_proteome_PFAS_wide$smoking_at_enrollment,
                                                        levels = c("No", "Yes"))

## status
BioMe_proteome_PFAS_wide$status<- factor(ifelse(BioMe_proteome_PFAS_wide$status=="case", 1, 0),
                                         levels = c(0,1))

## stratified by status
BioMe_proteome_PFAS_wide_case<- BioMe_proteome_PFAS_wide %>% 
                                filter(td2_case_all == 1)

BioMe_proteome_PFAS_wide_control<- BioMe_proteome_PFAS_wide %>% 
                                filter(td2_case_all == 0)




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



##########################
## function for model fit
##########################

lm_fit_info<- function(protein, data, data_in_long, covariates, path){
  
  PFAS_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)
  
  for(i in 1:length(protein)){
    
    s_lm <- (lm(as.formula(paste0(protein[i], covariates)), 
                data = data, weights = ipw))
    
    cov.m1 <- vcovHC(s_lm, type = "HC3")
    
    std.err <- sqrt(diag(cov.m1))
    
    r.est <- cbind(
      Estimate = coef(s_lm)
      , "Robust SE" = std.err
      , z = (coef(s_lm)/std.err)
      , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
    
    
    PFAS_lm <- rbind(PFAS_lm, c(protein[i], as.numeric(r.est[2, c(1, 2, 3, 4)])))
  }
  
  
  
  PFAS_lm <- PFAS_lm[-1,]
  PFAS_lm$z.value <- as.numeric(PFAS_lm$z.value)
  PFAS_lm$p.value <- as.numeric(PFAS_lm$p.value)
  PFAS_lm$fdr <- as.numeric(p.adjust(PFAS_lm$p.value, method="fdr"))
  
  
  q <- qvalue::qvalue(as.numeric(PFAS_lm$p.value), lambda=0)
  PFAS_lm$q.value <-  q$qvalues
  
  
  PFAS_lm_results<- PFAS_lm %>% 
    left_join(data_in_long[,c("OlinkID", "Protein_name", "UniProt")], by="OlinkID") %>% 
    distinct()
  
  write.csv(PFAS_lm_results,
            path,
            row.names = F)
}


##########################

#-------------------------------------------------- whole sample
#----------------------- Continuous PFDA - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFDA_Aug21_q + self_reported_race + age_at_enrollment + sex  + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFDA_allpanel_adlm_q.csv")


#----------------------- Continuous PFOA - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFOA_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOA_allpanel_adlm_q.csv")


#----------------------- Continuous PFOS - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFOS_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOS_allpanel_adlm_q.csv")



#----------------------- Continuous PFHxS - quartile
#---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFHxS_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHxS_allpanel_adlm_q.csv")



#----------------------- Continuous PFNA - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFNA_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFNA_allpanel_adlm_q.csv")



#----------------------- Continuous PFHpS - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFHpS_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpS_allpanel_adlm_q.csv")

#----------------------- Continuous PFHpA - binary
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFHpA_Aug21_bi + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpA_allpanel_adlm_bi.csv")

#-------------------------------------------------- in case
#----------------------- Continuous PFDA - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_case, BioMe_proteome_PFAS_long, 
            "~ PFDA_Aug21_q + self_reported_race + age_at_enrollment + sex  + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFDA_allpanel_adlm_q_case.csv")


#----------------------- Continuous PFOA - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_case, BioMe_proteome_PFAS_long, 
            "~ PFOA_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOA_allpanel_adlm_q_case.csv")


#----------------------- Continuous PFOS - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_case, BioMe_proteome_PFAS_long, 
            "~ PFOS_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOS_allpanel_adlm_q_case.csv")



#----------------------- Continuous PFHxS - quartile
#---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_case, BioMe_proteome_PFAS_long, 
            "~ PFHxS_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHxS_allpanel_adlm_q_case.csv")



#----------------------- Continuous PFNA - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_case, BioMe_proteome_PFAS_long, 
            "~ PFNA_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFNA_allpanel_adlm_q_case.csv")



#----------------------- Continuous PFHpS - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_case, BioMe_proteome_PFAS_long, 
            "~ PFHpS_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpS_allpanel_adlm_q_case.csv")

#----------------------- Continuous PFHpA - binary
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_case, BioMe_proteome_PFAS_long, 
            "~ PFHpA_Aug21_bi + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpA_allpanel_adlm_bi_case.csv")

#-------------------------------------------------- in control
#----------------------- Continuous PFDA - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_control, BioMe_proteome_PFAS_long, 
            "~ PFDA_Aug21_q + self_reported_race + age_at_enrollment + sex  + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFDA_allpanel_adlm_q_control.csv")


#----------------------- Continuous PFOA - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_control, BioMe_proteome_PFAS_long, 
            "~ PFOA_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOA_allpanel_adlm_q_control.csv")


#----------------------- Continuous PFOS - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_control, BioMe_proteome_PFAS_long, 
            "~ PFOS_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOS_allpanel_adlm_q_control.csv")



#----------------------- Continuous PFHxS - quartile
#---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_control, BioMe_proteome_PFAS_long, 
            "~ PFHxS_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHxS_allpanel_adlm_q_control.csv")



#----------------------- Continuous PFNA - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_control, BioMe_proteome_PFAS_long, 
            "~ PFNA_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFNA_allpanel_adlm_q_control.csv")



#----------------------- Continuous PFHpS - quartile
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_control, BioMe_proteome_PFAS_long, 
            "~ PFHpS_Aug21_q + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpS_allpanel_adlm_q_control.csv")

#----------------------- Continuous PFHpA - binary
##---------------------- Adjusted

lm_fit_info(protein_in_allpanels, BioMe_proteome_PFAS_wide_control, BioMe_proteome_PFAS_long, 
            "~ PFHpA_Aug21_bi + self_reported_race + age_at_enrollment + sex + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpA_allpanel_adlm_bi_control.csv")



