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

## restrict to protein in inflammation panel
BioMe_proteome_PFAS_inflammation<-  BioMe_proteome_PFAS_long %>% 
                                    filter(Panel == "Inflammation")

## protein in inflammation panel
protein_in_inflammation<- (protein_in_panel %>% filter(Panel == "Inflammation"))$protein



## Binary PFDA
BioMe_proteome_PFAS_wide$PFDA_Aug21_bi<- factor(BioMe_proteome_PFAS_wide$PFDA_Aug21_bi,
                                                levels = c( "Lower PFDA", "Higher PFDA"))

## date of blood draw
BioMe_proteome_PFAS_wide$date_enrl <- rep(NA_real_, nrow(BioMe_proteome_PFAS_wide))
BioMe_proteome_PFAS_wide$month_yr_enrl<- as.character(BioMe_proteome_PFAS_wide$month_yr_enrl)


for(i in 1:nrow(BioMe_proteome_PFAS_wide)){
  
  x <- anytime::anydate(paste((strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][2]), " 1,", 2000 + as.numeric(strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][1])))
  mydates <- as.Date(c("2011-01-01"))
  BioMe_proteome_PFAS_wide$date_enrl[i] <- as.numeric((x - mydates[1])/365 )
  BioMe_proteome_PFAS_wide$year_enrl[i] <- round(2000 + as.numeric(strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][1]), 0)
  
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



## Table
demo_table<- BioMe_proteome_PFAS_wide %>% 
             tbl_summary(
               include = c(self_reported_race, sex, age_at_enrollment, smoking_at_enrollment, bmi_at_enrollment, year_enrl),
               type = list(smoking_at_enrollment ~ "categorical"),
               digits = year_enrl~ 0
             )

demo_table


## analysis dataset
# analysis_data<- BioMe_proteome_PFAS_wide %>% 
#                 select(starts_with("OID"), self_reported_race, sex, age_at_enrollment, smoking_at_enrollment, bmi_at_enrollment, date_enrl, PFDA_Aug21_bi)

###### before dummify, we need to handle missing carefully in categorical variables
# analysis_data<- as.data.frame(DataExplorer::dummify(analysis_data))




################################
# function for correlation plot
################################

#--- corrplot
corr_fun<- function(data){
  corrplot(data, 
           method="color" ,
           col = colorRampPalette(c("steelblue", "white", "darkred"))(100),cl.lim=c(0,1),
           type="full",
           # order="hclust" ,
           # addrect = 2,
           tl.pos = 'lt',
           tl.srt=30,
           tl.col = "black",
  )
}

################################


# cor_data<- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')

# cor_plot<- corr_fun(cor_data) 

##########################
## function for model fit
##########################

lm_fit_info<- function(protein, data, data_in_long, covariates, path){

PFAS_lm <- data.frame(OlinkID = NA_character_, Value = NA_real_, Std.Error = NA_real_, z.value = NA_real_ , p.value = NA_real_)

for(i in 1:length(protein)){
  
  s_lm <- (lm(as.formula(paste0(protein[i], covariates)), 
              data = data))
  
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


q <-qvalue::qvalue(as.numeric(PFAS_lm$p.value), lambda=0)
PFAS_lm$q.value <-  q$qvalues


PFAS_lm_results<- PFAS_lm %>% 
  left_join(data_in_long[,c("OlinkID", "Protein_name")], by="OlinkID") %>% 
  distinct()

write.csv(PFAS_lm_results,
          path,
          row.names = F)
}

##########################



#----------------------- Binary PFDA
##---------------------- Uadjusted


lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFDA_Aug21_bi", "~/Projects/BioMe/proteome/input/exwas/check.csv")

##---------------------- Adjusted

lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFDA_Aug21_bi + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_adlm.csv")


#----------------------- Continuous PFDA
##---------------------- Uadjusted


lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFDA_Aug21", "~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_unlm_cont.csv")

##---------------------- Adjusted

lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFDA_Aug21 + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_adlm_cont.csv")


#----------------------- Continuous PFDA - tertile
##---------------------- Uadjusted


lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFDA_Aug21_q", "~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_unlm_q.csv")

##---------------------- Adjusted

lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFDA_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_adlm_q.csv")


#----------------------- Continuous PFOA - tertile
##---------------------- Uadjusted


lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFOA_Aug21_q", "~/Projects/BioMe/proteome/input/exwas/exwas_PFOA_inflam_unlm_q.csv")

##---------------------- Adjusted

lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFOA_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/exwas_PFOA_inflam_adlm_q.csv")


#----------------------- Continuous PFOS - tertile
##---------------------- Uadjusted


lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFOS_Aug21_q", "~/Projects/BioMe/proteome/input/exwas/exwas_PFOS_inflam_unlm_q.csv")

##---------------------- Adjusted

lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFOS_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/exwas_PFOS_inflam_adlm_q.csv")



#----------------------- Continuous PFHpA - tertile
##---------------------- Uadjusted


lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFHpA_Aug21_q", "~/Projects/BioMe/proteome/input/exwas/CHECK.csv")

##---------------------- Adjusted

lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFHpA_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/exwas_PFHpA_inflam_adlm_q.csv")



#----------------------- Continuous PFHxS - tertile
##---------------------- Uadjusted


lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFHxS_Aug21_q", "~/Projects/BioMe/proteome/input/exwas/exwas_PFHxS_inflam_unlm_q.csv")

##---------------------- Adjusted

lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFHxS_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/exwas_PFHxS_inflam_adlm_q.csv")




#----------------------- Continuous PFNA - tertile
##---------------------- Uadjusted


lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFNA_Aug21_q", "~/Projects/BioMe/proteome/input/exwas/exwas_PFNA_inflam_unlm_q.csv")

##---------------------- Adjusted

lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFNA_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/exwas_PFNA_inflam_adlm_q.csv")





#----------------------- Continuous PFHpS - tertile
##---------------------- Uadjusted


lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFHpS_Aug21_q", "~/Projects/BioMe/proteome/input/exwas/exwas_PFHpS_inflam_unlm_q.csv")

##---------------------- Adjusted

lm_fit_info(protein_in_inflammation, BioMe_proteome_PFAS_wide, BioMe_proteome_PFAS_long, 
            "~ PFHpS_Aug21_q + self_reported_race + age_at_enrollment + sex + bmi_at_enrollment + c_date_enrl", "~/Projects/BioMe/proteome/input/exwas/exwas_PFHpS_inflam_adlm_q.csv")

















