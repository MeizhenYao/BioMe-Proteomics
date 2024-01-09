library(car)
library(readr)
library(lattice)
library(nlme)
library(ggplot2)
library(GGally)
library(nnet)
library(foreign)
library(biotools)
library(glmmML)
library(MASS)
library(lme4)
library(multcomp)
library(dplyr)
library(knitr)
library(xtable)
library(kableExtra)
library(DT)
library(glmnet)
library(corrplot)
library(ggpubr)
library(lmerTest)
library("merTools")
library(reshape2)
library(ggplot2)
library(GGally)
library(mgcv)
library(gplots)
library(tidyr)
library(bkmr)
library(factoextra) 
library(spatstat)
library(Hmisc)
library(gtsummary)
library(blme)
library(grpreg)
library(robustHD)
library(gWQS)
library(gridExtra)
library(ggcorrplot)
library(BWQS)
library(qwraps2)
library(MatchIt)
library(data.table)
library(mice)
library(ggrepel)

##########################################################################################

setwd("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/Documents/Projects/BioMe/metabolome")


#################################################################################

## Epi data
data_all <- read.csv("merged_pfas_epi_liver_data.csv")
data_all <- data_all[,!(colnames(data_all) %in% c("NAFLD__status"))]
data_all$SAMPLEID <- as.character(data_all$SAMPLEID) 


ff1 <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/Documents/Projects/BioMe/metabolome/ff1_matching.csv")
dim(ff1)

ff1_prob <- ff1[!is.na(ff1$td2_case_all),]
dim(ff1_prob)

ff1_prob <- ff1_prob[ff1_prob$self_reported_race == "African American" |ff1_prob$self_reported_race == "European American" 
                     |ff1_prob$self_reported_race == "Hispanic", ]
dim(ff1_prob)

ff1_prob <- ff1_prob[ff1_prob$gender == "Female"| ff1_prob$gender == "Male", ]
dim(ff1_prob)

ff1_prob <- ff1_prob[ff1_prob$td2_case_incident == 0| ff1_prob$td2_case_incident == 2, ]
dim(ff1_prob)

ff1_prob$participant <- ff1_prob$masked_mrn %in% data_all$masked_mrn + 0


# estimation of denominator of ip weights

denom.fit <- glm(participant ~ as.character(self_reported_race) + as.factor(td2_case_all) + age_at_enrollment + gender, 
                 family = binomial(), data = ff1_prob)
summary(denom.fit)

pd.qsmk <- predict(denom.fit, type = "response")

# estimation of numerator of ip weights
numer.fit <- glm(participant~1, family = binomial(), data = ff1_prob)
summary(numer.fit)

pn.qsmk <- predict(numer.fit, type = "response")

ff1_prob$sw <- ifelse(ff1_prob$participant == 0, ((1-pn.qsmk)/(1-pd.qsmk)),
                      (pn.qsmk/pd.qsmk))

summary(ff1_prob$sw)

rrt <- ff1_prob[ff1_prob$masked_mrn %in% data_all$masked_mrn, c("sw","masked_mrn")]

# Final Inverse probability weights adjusted for case-control design

data_all$sw <- rrt[order(rrt$masked_mrn,data_all$masked_mrn),"sw"]



##################################################################################
## negative

mapfile <- read.delim("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/Documents/Projects/BioMe/metabolome/HRE0013_BioMe_AllModes_Mapfile_14May21.txt")
refmet_neg <- read.delim("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/Documents/Projects/BioMe/metabolome/RefMet_C18-HILICneg_Targets_7.5ppm_10percent-Studies_14May21.txt")


min(refmet_neg$Max_number_studies_reported) # All metabolites were reported in at least 14 
refmet_neg <- refmet_neg[,colnames(refmet_neg) %in% c(colnames(refmet_neg)[1:17],mapfile$Sample.ID[mapfile$Sample_Class == "Study_Sample"])]
# rownames(refmet_neg) <- refmet_neg$refmet_name_mode_id


# 1. metabolites with at least 75% non-zero values
y <- apply(refmet_neg[,c(18:376)], 1, function(x){sum(as.numeric(x) > 0)/length(as.numeric(x))})
refmet_neg <- refmet_neg[which(y>0.75),] 

# 2. keep unique metabolites
d_neg <- refmet_neg[,c(3,18:376)]
uniq <- unique(d_neg$refmet_name)
ccv <- NA_character_
for(i in 1:length(uniq)){
  ccv <- c(ccv, (names(which.max(apply(d_neg[which(d_neg$refmet_name == uniq[i]), -c(1)], 1, function(x){sd(as.numeric(x), na.rm = T)/mean(as.numeric(x), na.rm = T)})))))
}
ccv <- ccv[-1]
refmet_neg <- refmet_neg[ccv,]
refmet_neg$met <- paste0("neg",seq(1,dim(refmet_neg)[1]))
rownames(refmet_neg) <- paste0("neg",seq(1,dim(refmet_neg)[1]))

# 3. transpose dataset to make sure columns are metabolites, and rows are subject
t_refmet_neg <- data.table::transpose(refmet_neg)

# get row and colnames in order
colnames(t_refmet_neg) <- rownames(refmet_neg)
rownames(t_refmet_neg) <- colnames(refmet_neg)
t_refmet_neg <- t_refmet_neg[!(rownames(t_refmet_neg) %in% c(rownames(t_refmet_neg)[1:17],"met")),]

t_refmet_neg$SAMPLEID <- rownames(t_refmet_neg)
for(i in 1:length(t_refmet_neg$SAMPLEID)){
  if(length(strsplit(t_refmet_neg$SAMPLEID,"AB")[i][[1]]) == 2){
    t_refmet_neg$SAMPLEID[i] = (strsplit(t_refmet_neg$SAMPLEID,"AB")[i][[1]])[2]
  }
  else if (length(strsplit(t_refmet_neg$SAMPLEID,"X")[i][[1]]) == 2){
    t_refmet_neg$SAMPLEID[i] = (strsplit(t_refmet_neg$SAMPLEID,"X")[i][[1]])[2]
  }
  
}

t_refmet_neg <- apply(t_refmet_neg, 2, as.numeric)
t_refmet_neg <- as.data.frame(t_refmet_neg)
t_refmet_neg$SAMPLEID <- as.character(t_refmet_neg$SAMPLEID)


# 3. Change Metabolites into log2 + scale
chunk_18 <- apply(t_refmet_neg[,!(colnames(t_refmet_neg) %in% c("SAMPLEID"))], 2, function(x) scale(log(x + 1, base = 2)))
chunk_18 <- as.data.frame(chunk_18)
chunk_18 <- as.data.frame(scale(chunk_18))
chunk_18$SAMPLEID <-  t_refmet_neg$SAMPLEID
t_refmet_neg <- as.data.frame(chunk_18)
metabolites_refmet_neg <- colnames(t_refmet_neg)[colnames(t_refmet_neg) != "SAMPLEID"]

data_all__neg <- merge(data_all,t_refmet_neg, by = "SAMPLEID", all.x	= T, all.y = T)

# 4. Merge with more updated New PFAS - August 2021

aug_pfas <- read.csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/Documents/Projects/BioMe/metabolome/New PFAS data Aug 21.csv")
aug_pfas <- aug_pfas[(aug_pfas$Sample_Class %in% c("Study_Sample")),]

for(i in 1:length(aug_pfas$SAMPLEID)){
  if(length(strsplit(aug_pfas$SAMPLEID,"AB")[i][[1]]) == 2){
    aug_pfas$SAMPLEID[i] = (strsplit(aug_pfas$SAMPLEID,"AB")[i][[1]])[2]
  }
  else if (length(strsplit(aug_pfas$SAMPLEID,"X")[i][[1]]) == 2){
    aug_pfas$SAMPLEID[i] = (strsplit(aug_pfas$SAMPLEID,"X")[i][[1]])[2]
  }
  
}

data_all__neg <- merge(data_all__neg,aug_pfas, by = "SAMPLEID")


df <- data_all__neg[,colnames(data_all__neg) %in% c("PRS_score","age_at_enrollment","self_reported_race",
                                                    "smoking_at_enrollment","gender","bmi_at_enrollment","glucose_enrl","age_at_t2d",
                                                    "PFDA_Aug21","PFHpA_Aug21","PFHxS_Aug21","PFHpS_Aug21","PFNA_Aug21","PFOA_Aug21",
                                                    "PFOS_Aug21","status","selection_prop", "month_yr_enrl")]
init = mice(df, maxit=0) 
meth = init$method
predM = init$predictorMatrix
# predM[, c("redcap_survey_identifier")]=0

set.seed(1234)
imputed = mice(df, method=meth, predictorMatrix=predM, m=10)
imputed <- mice::complete(imputed, action = 10)
data_all__neg$glucose_enrl_imputed <- imputed$glucose_enrl
data_all__neg$bmi_at_enrollment_imputed <- imputed$bmi_at_enrollment
dim(data_all__neg)

data_pred <- data_all__neg[,colnames(data_all__neg) %in% c("SAMPLEID","DID","PRS_score","age_at_enrollment","self_reported_race",
                                                           "smoking_at_enrollment","gender","bmi_at_enrollment_imputed","glucose_enrl_imputed","age_at_t2d",
                                                           "PFDA_Aug21","PFHpA_Aug21","PFHxS_Aug21","PFHpS_Aug21","PFNA_Aug21","PFOA_Aug21",
                                                           "PFOS_Aug21","status","selection_prop","month_yr_enrl","sw",metabolites_refmet_neg)]
data_pred$status <- ifelse(data_pred$status == "case", 1, 0)
data_pred$date_enrl <- rep(NA_real_, nrow(data_pred))
for(i in 1:nrow(data_pred)){
  
  x <- anytime::anydate(paste((strsplit(data_pred$month_yr_enrl[i],"-")[[1]][2]), " 1,", 2000 + as.numeric(strsplit(data_pred$month_yr_enrl[i],"-")[[1]][1])))
  mydates <- as.Date(c("2011-01-01"))
  data_pred$date_enrl[i] <- as.numeric((x - mydates[1])/365 )
  
}



##############################
d_logis_dummy <- as.data.frame(data_pred[,!(colnames(data_pred) %in% c("month_yr_enrl"))])
new_quantile <- function(x, cuts ){
  
  qi <- unique(quantile(x, probs = seq(0, 1, by = 1/cuts), type = 2, na.rm = TRUE))
  
  if(length(qi) == 1){ 
    qi = c(-Inf, qi)
  } else{ 
    qi[1] <- -Inf
    qi[length(qi)] <- Inf
  }
  
  x[which(!is.na(x))] = cut(x[!is.na(x)], breaks = qi, labels = FALSE, include.lowest = TRUE)
  
  return(x)
  
}





# 5. Change PFAS into quartile


qqt <- as.data.frame(apply(d_logis_dummy[,c("PFDA_Aug21","PFHpA_Aug21","PFHxS_Aug21","PFHpS_Aug21","PFNA_Aug21","PFOA_Aug21",
                                            "PFOS_Aug21")], 2, function(x) new_quantile(x, cuts = 4)))
colnames(qqt) <- paste0(c("PFDA_Aug21","PFHpA_Aug21","PFHxS_Aug21","PFHpS_Aug21","PFNA_Aug21","PFOA_Aug21",
                          "PFOS_Aug21"), "_q")
d_logis_dummy <- cbind(d_logis_dummy, qqt)




d_logis_dummy$c_date_enrl <- ifelse(d_logis_dummy$date_enrl > 0, 1,0)



# 6. only restricted to participants who have both proteome and metabolome data
proteome<- fread("BioMe_proteome_PFAS_wide_imputed.txt")

PRO_DID<- proteome$DID

metabolome_restrict<- d_logis_dummy %>% 
                      filter(DID %in% PRO_DID)


write.table(refmet_neg, "~/Projects/BioMe/proteome/input/analysis_sample/met/BioMe_metabolome_neg_PFAS_annotation.txt", row.names = FALSE)
write.table(metabolome_restrict, "~/Projects/BioMe/proteome/input/analysis_sample/met/BioMe_metabolome_neg_PFAS_wide.txt", row.names = FALSE)


