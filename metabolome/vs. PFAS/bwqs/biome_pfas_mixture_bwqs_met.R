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
met_vs_pfas_bwqs_neg <- fread("~/Projects/BioMe/proteome/input/exwas/metabolome/met_vs_pfas_bwqs_neg.txt")
met_vs_pfas_bwqs_pos <- fread("~/Projects/BioMe/proteome/input/exwas/metabolome/met_vs_pfas_bwqs_pos.txt")

annotation_neg <- fread("~/Projects/BioMe/proteome/input/analysis_sample/met/BioMe_metabolome_neg_PFAS_annotation.txt")
annotation_pos <- fread("~/Projects/BioMe/proteome/input/analysis_sample/met/BioMe_metabolome_pos_PFAS_annotation.txt")



##------------------------------------------- negative
proteome_vs_pfas_bwqs<- met_vs_pfas_bwqs_neg


se = (proteome_vs_pfas_bwqs$upper - proteome_vs_pfas_bwqs$mean)/1.96
pval <- pnorm(abs(proteome_vs_pfas_bwqs$mean/se), lower.tail = F) 
proteome_vs_pfas_bwqs$p.value <- pval
q <- qvalue::qvalue(pval, lambda = 0)
proteome_vs_pfas_bwqs$q.value <-  q$qvalues


proteome_vs_pfas_bwqs_neg<- proteome_vs_pfas_bwqs %>% 
  left_join(annotation_neg[,c("met", "refmet_name", "super_class", "main_class", "sub_class", "pubchem_cid")], by="met") %>% 
  rowwise() %>% 
  mutate(Pubchem = strsplit(pubchem_cid, split=";")[[1]][1])

# 
prot<- proteome_vs_pfas_bwqs_neg %>%
       arrange(p.value) %>%
       filter(q.value < 0.1)
# 
# prot_neg<-  proteome_vs_pfas_bwqs_neg %>%
#             arrange(p.value) %>%
#             filter(q.value < 0.05 & mean < 0)
# 
# prot_pos<-  proteome_vs_pfas_bwqs_neg %>%
#             arrange(p.value) %>%
#             filter(q.value < 0.05 & mean > 0)

prot_pvalue<- proteome_vs_pfas_bwqs_neg %>%
              arrange(p.value) %>%
              filter(p.value < 0.05)

write.csv(prot_pvalue, "~/Projects/BioMe/proteome/input/exwas/metabolome/IPA/protein_sig_pvalue_neg.csv", row.names = FALSE)
write.csv(prot, "~/Projects/BioMe/proteome/input/exwas/metabolome/IPA/protein_sig_neg.csv", row.names = FALSE)
write.csv(proteome_vs_pfas_bwqs_neg, "~/Projects/BioMe/proteome/input/exwas/metabolome/IPA/metabolites_pvalue_neg.csv", row.names = FALSE)

write.table(proteome_vs_pfas_bwqs_neg, "~/Projects/BioMe/proteome/input/exwas/metabolome/met_vs_pfas_bwqs_plot_neg.txt", row.names = FALSE)



##------------------------------------------- positive
proteome_vs_pfas_bwqs<- met_vs_pfas_bwqs_pos


se = (proteome_vs_pfas_bwqs$upper - proteome_vs_pfas_bwqs$mean)/1.96
pval <- pnorm(abs(proteome_vs_pfas_bwqs$mean/se), lower.tail = F) 
proteome_vs_pfas_bwqs$p.value <- pval
q <- qvalue::qvalue(pval, lambda = 0)
proteome_vs_pfas_bwqs$q.value <-  q$qvalues


proteome_vs_pfas_bwqs_pos<- proteome_vs_pfas_bwqs %>% 
  left_join(annotation_pos[,c("met", "refmet_name", "super_class", "main_class", "sub_class", "pubchem_cid")], by="met") %>% 
  rowwise() %>% 
  mutate(Pubchem = strsplit(pubchem_cid, split=";")[[1]][1]) 

# 
prot<- proteome_vs_pfas_bwqs_pos %>%
  arrange(p.value) %>%
  filter(q.value < 0.1)
# 
# prot_neg<-  proteome_vs_pfas_bwqs_pos %>%
#   arrange(p.value) %>%
#   filter(q.value < 0.05 & mean < 0)
# 
# prot_pos<-  proteome_vs_pfas_bwqs_pos %>%
#   arrange(p.value) %>%
#   filter(q.value < 0.05 & mean > 0)

prot_pvalue<- proteome_vs_pfas_bwqs_pos %>%
  arrange(p.value) %>%
  filter(p.value < 0.05)

write.csv(prot_pvalue, "~/Projects/BioMe/proteome/input/exwas/metabolome/IPA/protein_sig_pvalue_pos.csv", row.names = FALSE)
write.csv(prot, "~/Projects/BioMe/proteome/input/exwas/metabolome/IPA/protein_sig_pos.csv", row.names = FALSE)
write.csv(proteome_vs_pfas_bwqs_pos, "~/Projects/BioMe/proteome/input/exwas/metabolome/IPA/metabolites_pvalue_pos.csv", row.names = FALSE)

write.table(proteome_vs_pfas_bwqs_pos, "~/Projects/BioMe/proteome/input/exwas/metabolome/met_vs_pfas_bwqs_plot_pos.txt", row.names = FALSE)





#############
## weights ##
#############

##------------------------------------------- import data
met_vs_pfas_bwqs_neg <- fread("~/Projects/BioMe/proteome/input/exwas/metabolome/met_pfas_weight_neg.txt")
met_vs_pfas_bwqs_pos <- fread("~/Projects/BioMe/proteome/input/exwas/metabolome/met_pfas_weight_pos.txt")

annotation_neg <- fread("~/Projects/BioMe/proteome/input/analysis_sample/met/BioMe_metabolome_neg_PFAS_annotation.txt")
annotation_pos <- fread("~/Projects/BioMe/proteome/input/analysis_sample/met/BioMe_metabolome_pos_PFAS_annotation.txt")


##------------------------------------------- negative
bwqs_pfas_weight<- met_vs_pfas_bwqs_neg


old_name<- colnames(bwqs_pfas_weight)
new_name<- c("PFDA", "PFHxS", "PFHpS", "PFNA", "PFOA", "PFOS", "met")



setnames(bwqs_pfas_weight, old = old_name, new = new_name)


bwqs_pfas_weight_neg<- bwqs_pfas_weight %>% 
                       left_join(annotation_neg[,c("met", "refmet_name", "super_class", "main_class", "sub_class")], by="met")  


write.table(bwqs_pfas_weight_neg, "~/Projects/BioMe/proteome/input/exwas/metabolome/met_pfas_weight_plot_neg.txt", row.names = FALSE)


##------------------------------------------- positive
bwqs_pfas_weight<- met_vs_pfas_bwqs_pos


old_name<- colnames(bwqs_pfas_weight)
new_name<- c("PFDA", "PFHxS", "PFHpS", "PFNA", "PFOA", "PFOS", "met")



setnames(bwqs_pfas_weight, old = old_name, new = new_name)


bwqs_pfas_weight_pos<- bwqs_pfas_weight %>% 
  left_join(annotation_pos[,c("met", "refmet_name", "super_class", "main_class", "sub_class")], by="met")  


write.table(bwqs_pfas_weight_pos, "~/Projects/BioMe/proteome/input/exwas/metabolome/met_pfas_weight_plot_pos.txt", row.names = FALSE)




