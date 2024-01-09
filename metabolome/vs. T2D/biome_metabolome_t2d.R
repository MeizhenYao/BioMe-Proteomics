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
library(sandwich)
# 
# cores=detectCores()
# cl <- makeCluster(10) 
# registerDoParallel(cl)
# 
# start.time <- Sys.time()


########### negative

##------------------------------------------- import data
BioMe_proteome_PFAS_wide <- fread("~/Projects/BioMe/proteome/input/analysis_sample/met/BioMe_metabolome_neg_PFAS_wide.txt")
annotation <- fread("~/Projects/BioMe/proteome/input/analysis_sample/met/BioMe_metabolome_neg_PFAS_annotation.txt")


##------------------------------------------- prepare data


analysis_data<- BioMe_proteome_PFAS_wide %>% 
  dplyr::select(starts_with("neg"), ends_with("_q"), status, self_reported_race, gender, age_at_enrollment, smoking_at_enrollment, c_date_enrl, sw)


analysis_data_dummy<- as.data.frame(dummify(analysis_data))



##------------------------------------------- bwqs model fitting

protein<- colnames(analysis_data_dummy %>% dplyr::select(starts_with("neg")))


data_proteins<- analysis_data_dummy %>% 
  dplyr::select(starts_with("neg"))


proteome_t2d_model <- data.frame(Value = NA_real_, Std.Error = NA_real_, t.value = NA_real_ , p.value = NA_real_)

bwqs_pfas_weight<- data.frame(w1 = NA_real_, w2 = NA_real_, w3 = NA_real_,
                              w4 = NA_real_, w5 = NA_real_, w6 = NA_real_)

for(i in 1:length(protein)){
  
  s_lm <- (glm(status ~ data_proteins[,i]
               + self_reported_race.African.American + self_reported_race.European.American 
               + age_at_enrollment
               + gender.Female
               + c_date_enrl , data = analysis_data_dummy, family=binomial))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  proteome_t2d_model <- rbind(proteome_t2d_model, as.numeric(r.est[2,c(1,2,3,4)]))
}

proteome_t2d_model <- proteome_t2d_model[-1,]
proteome_t2d_model$t.value <- as.numeric(proteome_t2d_model$t.value)
proteome_t2d_model$p.value <- as.numeric(proteome_t2d_model$p.value)


proteome_t2d_model$met<- protein

q <- qvalue::qvalue(as.numeric(proteome_t2d_model$p.value), lambda=0)
proteome_t2d_model$q.value <-  q$qvalues


proteome_t2d_model<- proteome_t2d_model %>% 
  left_join(annotation[,c("met", "refmet_name", "super_class", "main_class", "sub_class")], by="met") 



write.table(proteome_t2d_model, "~/Projects/BioMe/proteome/input/mwas/metabolome_vs_t2d_neg.txt", row.names = FALSE)

# proteome_t2d_model<- fread("~/Projects/BioMe/proteome/input/pwas/proteome_vs_t2d_scale.txt")
# 
# end.time <- Sys.time()
# (time.taken <- end.time - start.time)
# 
# stopCluster(cl)


###### extract data for pathway analysis

prot<-  proteome_t2d_model %>%
  arrange(p.value) %>%
  filter(q.value < 0.1)

prot_pos<-  proteome_t2d_model %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)

prot_neg<-  proteome_t2d_model %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)


prot1<-  proteome_t2d_model %>%
  arrange(p.value) %>%
  filter(p.value < 0.05)


write.csv(prot, "~/Projects/BioMe/proteome/input/mwas/IPA/protein_sig_neg.csv", row.names = FALSE)

write.csv(prot1, "~/Projects/BioMe/proteome/input/mwas/IPA/protein_sig_p_value_neg.csv", row.names = FALSE)



########### positive

##------------------------------------------- import data
BioMe_proteome_PFAS_wide <- fread("~/Projects/BioMe/proteome/input/analysis_sample/met/BioMe_metabolome_pos_PFAS_wide.txt")
annotation <- fread("~/Projects/BioMe/proteome/input/analysis_sample/met/BioMe_metabolome_pos_PFAS_annotation.txt")


##------------------------------------------- prepare data


analysis_data<- BioMe_proteome_PFAS_wide %>% 
  dplyr::select(starts_with("pos"), ends_with("_q"), status, self_reported_race, gender, age_at_enrollment, smoking_at_enrollment, c_date_enrl, sw)


analysis_data_dummy<- as.data.frame(dummify(analysis_data))



##------------------------------------------- bwqs model fitting

protein<- colnames(analysis_data_dummy %>% dplyr::select(starts_with("pos")))


data_proteins<- analysis_data_dummy %>% 
  dplyr::select(starts_with("pos"))


proteome_t2d_model <- data.frame(Value = NA_real_, Std.Error = NA_real_, t.value = NA_real_ , p.value = NA_real_)

bwqs_pfas_weight<- data.frame(w1 = NA_real_, w2 = NA_real_, w3 = NA_real_,
                              w4 = NA_real_, w5 = NA_real_, w6 = NA_real_)

for(i in 1:length(protein)){
  
  s_lm <- (glm(status ~ data_proteins[,i]
               + self_reported_race.African.American + self_reported_race.European.American 
               + age_at_enrollment
               + gender.Female
               + c_date_enrl , data = analysis_data_dummy, family=binomial))
  
  cov.m1 <- vcovHC(s_lm, type = "HC3")
  
  std.err <- sqrt(diag(cov.m1))
  
  r.est <- cbind(
    Estimate = coef(s_lm)
    , "Robust SE" = std.err
    , z = (coef(s_lm)/std.err)
    , "Pr(>|z|) "= 2 * pnorm(abs(coef(s_lm)/std.err), lower.tail = FALSE))
  
  
  proteome_t2d_model <- rbind(proteome_t2d_model, as.numeric(r.est[2,c(1,2,3,4)]))
}

proteome_t2d_model <- proteome_t2d_model[-1,]
proteome_t2d_model$t.value <- as.numeric(proteome_t2d_model$t.value)
proteome_t2d_model$p.value <- as.numeric(proteome_t2d_model$p.value)


proteome_t2d_model$met<- protein

q <- qvalue::qvalue(as.numeric(proteome_t2d_model$p.value), lambda=0)
proteome_t2d_model$q.value <-  q$qvalues


proteome_t2d_model<- proteome_t2d_model %>% 
  left_join(annotation[,c("met", "refmet_name", "super_class", "main_class", "sub_class")], by="met") 



write.table(proteome_t2d_model, "~/Projects/BioMe/proteome/input/mwas/metabolome_vs_t2d_pos.txt", row.names = FALSE)

# proteome_t2d_model<- fread("~/Projects/BioMe/proteome/input/pwas/proteome_vs_t2d_scale.txt")
# 
# end.time <- Sys.time()
# (time.taken <- end.time - start.time)
# 
# stopCluster(cl)


###### extract data for pathway analysis

prot<-  proteome_t2d_model %>%
  arrange(p.value) %>%
  filter(q.value < 0.1)

prot_pos<-  proteome_t2d_model %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)

prot_neg<-  proteome_t2d_model %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)


prot1<-  proteome_t2d_model %>%
  arrange(p.value) %>%
  filter(p.value < 0.05)


write.csv(prot, "~/Projects/BioMe/proteome/input/mwas/IPA/protein_sig_pos.csv", row.names = FALSE)

write.csv(prot1, "~/Projects/BioMe/proteome/input/mwas/IPA/protein_sig_p_value_pos.csv", row.names = FALSE)





