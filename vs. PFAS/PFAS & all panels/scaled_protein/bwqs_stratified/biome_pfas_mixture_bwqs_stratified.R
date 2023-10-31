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

## case
#####################
## effect estimate ##
#####################

##------------------------------------------- import data
proteome_vs_pfas_bwqs_cardio <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/proteome_vs_pfas_bwqs_cardio_case.txt")
proteome_vs_pfas_bwqs_cardio2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/proteome_vs_pfas_bwqs_cardio2_case.txt")
proteome_vs_pfas_bwqs_inflammation <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/proteome_vs_pfas_bwqs_inflammation_case.txt")
proteome_vs_pfas_bwqs_inflammation2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/proteome_vs_pfas_bwqs_inflammation2_case.txt")
proteome_vs_pfas_bwqs_neuro <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/proteome_vs_pfas_bwqs_neuro_case.txt")
proteome_vs_pfas_bwqs_neuro2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/proteome_vs_pfas_bwqs_neuro2_case.txt")
proteome_vs_pfas_bwqs_oncology <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/proteome_vs_pfas_bwqs_oncology_case.txt")
proteome_vs_pfas_bwqs_oncology2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/proteome_vs_pfas_bwqs_oncology2_case.txt")


protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")


##------------------------------------------- data processing
proteome_vs_pfas_bwqs<- rbind(proteome_vs_pfas_bwqs_cardio,
                              proteome_vs_pfas_bwqs_cardio2,
                              proteome_vs_pfas_bwqs_inflammation,
                              proteome_vs_pfas_bwqs_inflammation2,
                              proteome_vs_pfas_bwqs_neuro,
                              proteome_vs_pfas_bwqs_neuro2,
                              proteome_vs_pfas_bwqs_oncology,
                              proteome_vs_pfas_bwqs_oncology2)


se = (proteome_vs_pfas_bwqs$upper - proteome_vs_pfas_bwqs$mean)/1.96
pval <- pnorm(abs(proteome_vs_pfas_bwqs$mean/se), lower.tail = F) 
proteome_vs_pfas_bwqs$p.value <- pval
q <- qvalue::qvalue(pval, lambda = 0)
proteome_vs_pfas_bwqs$q.value <-  q$qvalues


proteome_vs_pfas_bwqs<- proteome_vs_pfas_bwqs %>% 
                        left_join(protein_in_panel[,c("OlinkID", "Protein_name", "UniProt", "Gene_name")], by="OlinkID") 


prot<- proteome_vs_pfas_bwqs %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)




write.csv(prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/protein_sig.csv", row.names = FALSE)
write.table(proteome_vs_pfas_bwqs, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/proteome_vs_pfas_bwqs_case.txt", row.names = FALSE)


#############
## weights ##
#############

##------------------------------------------- import data
bwqs_pfas_weight_cardio <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/bwqs_pfas_weight_cardio_case.txt")
bwqs_pfas_weight_cardio2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/bwqs_pfas_weight_cardio2_case.txt")
bwqs_pfas_weight_inflammation <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/bwqs_pfas_weight_inflammation_case.txt")
bwqs_pfas_weight_inflammation2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/bwqs_pfas_weight_inflammation2_case.txt")
bwqs_pfas_weight_neuro <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/bwqs_pfas_weight_neuro_case.txt")
bwqs_pfas_weight_neuro2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/bwqs_pfas_weight_neuro2_case.txt")
bwqs_pfas_weight_oncology <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/bwqs_pfas_weight_oncology_case.txt")
bwqs_pfas_weight_oncology2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/bwqs_pfas_weight_oncology2_case.txt")

protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")


##------------------------------------------- data processing
bwqs_pfas_weight<- data.frame(rbind(bwqs_pfas_weight_cardio,
                                    bwqs_pfas_weight_cardio2,
                                    bwqs_pfas_weight_inflammation,
                                    bwqs_pfas_weight_inflammation2,
                                    bwqs_pfas_weight_neuro,
                                    bwqs_pfas_weight_neuro2,
                                    bwqs_pfas_weight_oncology,
                                    bwqs_pfas_weight_oncology2))


old_name<- colnames(bwqs_pfas_weight)
new_name<- c("PFDA", "PFHxS", "PFHpS", "PFNA", "PFOA", "PFOS", "OlinkID")



setnames(bwqs_pfas_weight, old = old_name, new = new_name)


bwqs_pfas_weight<- bwqs_pfas_weight %>% 
                   left_join(protein_in_panel[,c("OlinkID", "Protein_name", "UniProt", "Gene_name")], by="OlinkID")


write.table(bwqs_pfas_weight, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/bwqs_pfas_weight_case.txt", row.names = FALSE)




## control
#####################
## effect estimate ##
#####################

##------------------------------------------- import data
proteome_vs_pfas_bwqs_cardio <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/proteome_vs_pfas_bwqs_cardio_control.txt")
proteome_vs_pfas_bwqs_cardio2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/proteome_vs_pfas_bwqs_cardio2_control.txt")
proteome_vs_pfas_bwqs_inflammation <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/proteome_vs_pfas_bwqs_inflammation_control.txt")
proteome_vs_pfas_bwqs_inflammation2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/proteome_vs_pfas_bwqs_inflammation2_control.txt")
proteome_vs_pfas_bwqs_neuro <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/proteome_vs_pfas_bwqs_neuro_control.txt")
proteome_vs_pfas_bwqs_neuro2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/proteome_vs_pfas_bwqs_neuro2_control.txt")
proteome_vs_pfas_bwqs_oncology <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/proteome_vs_pfas_bwqs_oncology_control.txt")
proteome_vs_pfas_bwqs_oncology2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/proteome_vs_pfas_bwqs_oncology2_control.txt")


protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")


##------------------------------------------- data processing
proteome_vs_pfas_bwqs<- rbind(proteome_vs_pfas_bwqs_cardio,
                              proteome_vs_pfas_bwqs_cardio2,
                              proteome_vs_pfas_bwqs_inflammation,
                              proteome_vs_pfas_bwqs_inflammation2,
                              proteome_vs_pfas_bwqs_neuro,
                              proteome_vs_pfas_bwqs_neuro2,
                              proteome_vs_pfas_bwqs_oncology,
                              proteome_vs_pfas_bwqs_oncology2)


se = (proteome_vs_pfas_bwqs$upper - proteome_vs_pfas_bwqs$mean)/1.96
pval <- pnorm(abs(proteome_vs_pfas_bwqs$mean/se), lower.tail = F) 
proteome_vs_pfas_bwqs$p.value <- pval
q <- qvalue::qvalue(pval, lambda = 0)
proteome_vs_pfas_bwqs$q.value <-  q$qvalues


proteome_vs_pfas_bwqs<- proteome_vs_pfas_bwqs %>% 
  left_join(protein_in_panel[,c("OlinkID", "Protein_name", "UniProt", "Gene_name")], by="OlinkID") 



prot<- proteome_vs_pfas_bwqs %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)




write.csv(prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/protein_sig.csv", row.names = FALSE)
write.table(proteome_vs_pfas_bwqs, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/proteome_vs_pfas_bwqs_control.txt", row.names = FALSE)


#############
## weights ##
#############

##------------------------------------------- import data
bwqs_pfas_weight_cardio <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/bwqs_pfas_weight_cardio_control.txt")
bwqs_pfas_weight_cardio2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/bwqs_pfas_weight_cardio2_control.txt")
bwqs_pfas_weight_inflammation <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/bwqs_pfas_weight_inflammation_control.txt")
bwqs_pfas_weight_inflammation2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/bwqs_pfas_weight_inflammation2_control.txt")
bwqs_pfas_weight_neuro <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/bwqs_pfas_weight_neuro_control.txt")
bwqs_pfas_weight_neuro2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/bwqs_pfas_weight_neuro2_control.txt")
bwqs_pfas_weight_oncology <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/bwqs_pfas_weight_oncology_control.txt")
bwqs_pfas_weight_oncology2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/bwqs_pfas_weight_oncology2_control.txt")

protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")


##------------------------------------------- data processing
bwqs_pfas_weight<- data.frame(rbind(bwqs_pfas_weight_cardio,
                                    bwqs_pfas_weight_cardio2,
                                    bwqs_pfas_weight_inflammation,
                                    bwqs_pfas_weight_inflammation2,
                                    bwqs_pfas_weight_neuro,
                                    bwqs_pfas_weight_neuro2,
                                    bwqs_pfas_weight_oncology,
                                    bwqs_pfas_weight_oncology2))


old_name<- colnames(bwqs_pfas_weight)
new_name<- c("PFDA", "PFHxS", "PFHpS", "PFNA", "PFOA", "PFOS", "OlinkID")



setnames(bwqs_pfas_weight, old = old_name, new = new_name)


bwqs_pfas_weight<- bwqs_pfas_weight %>% 
  left_join(protein_in_panel[,c("OlinkID", "Protein_name", "UniProt", "Gene_name")], by="OlinkID")


write.table(bwqs_pfas_weight, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/bwqs_pfas_weight_control.txt", row.names = FALSE)



