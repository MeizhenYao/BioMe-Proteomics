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


#####################
## effect estimate ##
#####################

##------------------------------------------- import data
proteome_vs_pfas_bwqs_cardio <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/proteome_vs_pfas_bwqs_cardio.txt")
proteome_vs_pfas_bwqs_cardio2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/proteome_vs_pfas_bwqs_cardio2.txt")
proteome_vs_pfas_bwqs_inflammation <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/proteome_vs_pfas_bwqs_inflammation.txt")
proteome_vs_pfas_bwqs_inflammation2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/proteome_vs_pfas_bwqs_inflammation2.txt")
proteome_vs_pfas_bwqs_neuro <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/proteome_vs_pfas_bwqs_neuro.txt")
proteome_vs_pfas_bwqs_neuro2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/proteome_vs_pfas_bwqs_neuro2.txt")
proteome_vs_pfas_bwqs_oncology <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/proteome_vs_pfas_bwqs_oncology.txt")
proteome_vs_pfas_bwqs_oncology2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/proteome_vs_pfas_bwqs_oncology2.txt")


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


# write.table(proteome_vs_pfas_bwqs, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/proteome_vs_pfas_bwqs.txt", row.names = FALSE)


#############
## weights ##
#############

##------------------------------------------- import data
bwqs_pfas_weight_cardio <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/bwqs_pfas_weight_cardio.txt")
bwqs_pfas_weight_cardio2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/bwqs_pfas_weight_cardio2.txt")
bwqs_pfas_weight_inflammation <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/bwqs_pfas_weight_inflammation.txt")
bwqs_pfas_weight_inflammation2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/bwqs_pfas_weight_inflammation2.txt")
bwqs_pfas_weight_neuro <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/bwqs_pfas_weight_neuro.txt")
bwqs_pfas_weight_neuro2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/bwqs_pfas_weight_neuro2.txt")
bwqs_pfas_weight_oncology <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/bwqs_pfas_weight_oncology.txt")
bwqs_pfas_weight_oncology2 <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/bwqs_pfas_weight_oncology2.txt")

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


# write.table(bwqs_pfas_weight, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/bwqs_pfas_weight.txt", row.names = FALSE)






