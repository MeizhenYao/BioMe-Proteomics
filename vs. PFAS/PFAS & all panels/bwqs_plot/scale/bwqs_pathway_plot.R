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



##------------------------------------------- import data
## whole
PFDA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFDA_allpanel_adlm_q.csv")
PFOA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOA_allpanel_adlm_q.csv")
PFOS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOS_allpanel_adlm_q.csv")
PFHxS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHxS_allpanel_adlm_q.csv")
PFNA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFNA_allpanel_adlm_q.csv")
PFHpS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpS_allpanel_adlm_q.csv")


## PFDA
PFDA_prot<- PFDA_all_adlm_results %>%
            arrange(p.value) %>%
            filter(q.value < 0.05)

PFDA_prot_neg<- PFDA_all_adlm_results %>%
                arrange(p.value) %>%
                filter(q.value < 0.05 & Value < 0)

PFDA_prot_pos<- PFDA_all_adlm_results %>%
                arrange(p.value) %>%
                filter(q.value < 0.05 & Value > 0)


write.csv(PFDA_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFDA_protein_sig.csv", row.names = FALSE)
write.csv(PFDA_prot_neg, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFDA_protein_sig_neg.csv", row.names = FALSE)
write.csv(PFDA_prot_pos, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFDA_protein_sig_pos.csv", row.names = FALSE)


## PFOA
PFOA_prot<- PFOA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFOA_prot_neg<- PFOA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFOA_prot_pos<- PFOA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFOA_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFOA_protein_sig.csv", row.names = FALSE)
write.csv(PFOA_prot_neg, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFOA_protein_sig_neg.csv", row.names = FALSE)
write.csv(PFOA_prot_pos, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFOA_protein_sig_pos.csv", row.names = FALSE)


## PFOS
PFOS_prot<- PFOS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFOS_prot_neg<- PFOS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFOS_prot_pos<- PFOS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFOS_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFOS_protein_sig.csv", row.names = FALSE)
write.csv(PFOS_prot_neg, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFOS_protein_sig_neg.csv", row.names = FALSE)
write.csv(PFOS_prot_pos, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFOS_protein_sig_pos.csv", row.names = FALSE)


## PFHxS
PFHxS_prot<- PFHxS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFHxS_prot_neg<- PFHxS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFHxS_prot_pos<- PFHxS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFHxS_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFHxS_protein_sig.csv", row.names = FALSE)
write.csv(PFHxS_prot_neg, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFHxS_protein_sig_neg.csv", row.names = FALSE)
write.csv(PFHxS_prot_pos, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFHxS_protein_sig_pos.csv", row.names = FALSE)



## PFNA
PFNA_prot<- PFNA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFNA_prot_neg<- PFNA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFNA_prot_pos<- PFNA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFNA_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFNA_protein_sig.csv", row.names = FALSE)
write.csv(PFNA_prot_neg, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFNA_protein_sig_neg.csv", row.names = FALSE)
write.csv(PFNA_prot_pos, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFNA_protein_sig_pos.csv", row.names = FALSE)


## PFHpS
PFHpS_prot<- PFHpS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFHpS_prot_neg<- PFHpS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFHpS_prot_pos<- PFHpS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFHpS_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFHpS_protein_sig.csv", row.names = FALSE)
write.csv(PFHpS_prot_neg, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFHpS_protein_sig_neg.csv", row.names = FALSE)
write.csv(PFHpS_prot_pos, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFHpS_protein_sig_pos.csv", row.names = FALSE)

##--------------------- cases
PFDA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFDA_allpanel_adlm_q_case.csv")
PFOA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOA_allpanel_adlm_q_case.csv")
PFOS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOS_allpanel_adlm_q_case.csv")
PFHxS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHxS_allpanel_adlm_q_case.csv")
PFNA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFNA_allpanel_adlm_q_case.csv")
PFHpS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpS_allpanel_adlm_q_case.csv")


## PFDA
PFDA_prot<- PFDA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFDA_prot_neg<- PFDA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFDA_prot_pos<- PFDA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFDA_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFDA_protein_sig_case.csv", row.names = FALSE)


## PFOA
PFOA_prot<- PFOA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFOA_prot_neg<- PFOA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFOA_prot_pos<- PFOA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFOA_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFOA_protein_sig_case.csv", row.names = FALSE)


## PFOS
PFOS_prot<- PFOS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFOS_prot_neg<- PFOS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFOS_prot_pos<- PFOS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFOS_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFOS_protein_sig_case.csv", row.names = FALSE)


## PFHxS
PFHxS_prot<- PFHxS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFHxS_prot_neg<- PFHxS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFHxS_prot_pos<- PFHxS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFHxS_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFHxS_protein_sig_case.csv", row.names = FALSE)



## PFNA
PFNA_prot<- PFNA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFNA_prot_neg<- PFNA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFNA_prot_pos<- PFNA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFNA_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFNA_protein_sig_case.csv", row.names = FALSE)


## PFHpS
PFHpS_prot<- PFHpS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFHpS_prot_neg<- PFHpS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFHpS_prot_pos<- PFHpS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFHpS_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFHpS_protein_sig_case.csv", row.names = FALSE)


##--------------------- control
PFDA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFDA_allpanel_adlm_q_control.csv")
PFOA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOA_allpanel_adlm_q_control.csv")
PFOS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOS_allpanel_adlm_q_control.csv")
PFHxS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHxS_allpanel_adlm_q_control.csv")
PFNA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFNA_allpanel_adlm_q_control.csv")
PFHpS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpS_allpanel_adlm_q_control.csv")


## PFDA
PFDA_prot<- PFDA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFDA_prot_neg<- PFDA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFDA_prot_pos<- PFDA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFDA_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFDA_protein_sig_control.csv", row.names = FALSE)


## PFOA
PFOA_prot<- PFOA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFOA_prot_neg<- PFOA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFOA_prot_pos<- PFOA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFOA_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFOA_protein_sig_control.csv", row.names = FALSE)


## PFOS
PFOS_prot<- PFOS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFOS_prot_neg<- PFOS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFOS_prot_pos<- PFOS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFOS_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFOS_protein_sig_control.csv", row.names = FALSE)


## PFHxS
PFHxS_prot<- PFHxS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFHxS_prot_neg<- PFHxS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFHxS_prot_pos<- PFHxS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFHxS_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFHxS_protein_sig_control.csv", row.names = FALSE)



## PFNA
PFNA_prot<- PFNA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFNA_prot_neg<- PFNA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFNA_prot_pos<- PFNA_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFNA_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFNA_protein_sig_control.csv", row.names = FALSE)


## PFHpS
PFHpS_prot<- PFHpS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05)

PFHpS_prot_neg<- PFHpS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value < 0)

PFHpS_prot_pos<- PFHpS_all_adlm_results %>%
  arrange(p.value) %>%
  filter(q.value < 0.05 & Value > 0)


write.csv(PFHpS_prot, "~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/pathway/PFHpS_protein_sig_control.csv", row.names = FALSE)

#######################
#######################
#######################
#######################
#######################
#######################
######## only bwqs
##------------------------------------------- import data
bwqs_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/scale/bwqs.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(cohort = "Total sample")
bwqs_sig_case <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/scale/bwqs_case.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(cohort = "Cases")
bwqs_sig_control <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/scale/bwqs_control.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(cohort = "Controls")


sig_pathway<- rbind(bwqs_sig,
                    bwqs_sig_case,
                    bwqs_sig_control)

sig_pathway$cohort<- factor(sig_pathway$cohort,
                          levels = c("Total sample", "Cases", "Controls"))

sig_pathway_plot<-  ggplot(sig_pathway, aes(x = cohort,  y = Pathway.name)) + 
  geom_point(aes(size = -log(Entities.pValue , base = 10)))+
  scale_x_discrete(limits=c("Total sample", "Cases", "Controls"))+
  scale_size_continuous(range = c(4,12)) +
  labs(x = NULL) + 
  theme_bw() +
  theme( 
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    axis.text.x= element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title=element_text(size=18,face="bold"))

sig_pathway_plot                  




jpeg("~/Projects/BioMe/proteome/output/reactome-pathway/bwqs_scale_pathwayplot.jpeg",
     units="in", width=30, height=14, res=500)

sig_pathway_plot

dev.off()










##------------------------------------------- import data
bwqs_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/scale/bwqs.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFAS Mixture")
PFDA_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/scale/PFDA.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFDA")
PFOS_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/scale/PFOS.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFOS")
PFHxS_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/scale/PFHxS.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFHxS")
PFNA_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/scale/PFNA.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFNA")
PFHpS_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/scale/PFHpS.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFHpS")


sig_pathway<- rbind(bwqs_sig,
                    PFDA_sig,
                    PFOS_sig,
                    PFHxS_sig,
                    PFNA_sig,
                    PFHpS_sig)

sig_pathway$PFAS<- factor(sig_pathway$PFAS,
                          levels = c("PFAS Mixture", "PFDA", "PFHpS", "PFHxS", "PFNA", "PFOA", "PFOS"))

sig_pathway_plot<-  ggplot(sig_pathway, aes(x = PFAS,  y = Pathway.name)) + 
                    geom_point(aes(size = -log(Entities.pValue , base = 10)))+
                    scale_x_discrete(limits=c("PFAS Mixture", "PFDA", "PFHpS", "PFHxS", "PFNA", "PFOA", "PFOS"))+
                    scale_size_continuous(range = c(4,12)) +
                    labs(x = NULL) + 
                    theme_bw() +
                    theme( 
                      legend.title = element_text(size = 14, face = "bold"),
                      legend.text = element_text(size = 14, face = "bold"),
                      axis.text.x= element_text(size = 18, face = "bold"),
                      axis.text.y = element_text(size = 18, face = "bold"),
                      axis.title=element_text(size=18,face="bold"))
                  
sig_pathway_plot                  
                               
             
             
             
jpeg("~/Projects/BioMe/proteome/output/reactome-pathway/pathwayplot.jpeg",
     units="in", width=30, height=14, res=500)

sig_pathway_plot

dev.off()

             
             
## negative directions
bwqs_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/bwqs/bwqs_allpanels_neg.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFAS Mixture")
PFDA_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/PFDA_neg.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFDA")
PFOS_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/PFOS_neg.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFOS")
PFHxS_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/PFHxS_neg.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFHxS")
PFNA_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/PFNA_neg.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFNA")
PFHpS_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/PFHpS_neg.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFHpS")


sig_pathway<- rbind(bwqs_sig,
                    PFDA_sig,
                    PFOS_sig,
                    PFHxS_sig,
                    PFNA_sig,
                    PFHpS_sig)

sig_pathway$PFAS<- factor(sig_pathway$PFAS,
                          levels = c("PFAS Mixture", "PFDA", "PFHpS", "PFHxS", "PFNA", "PFOA", "PFOS"))

sig_pathway_plot<-  ggplot(sig_pathway, aes(x = PFAS,  y = Pathway.name)) + 
  geom_point(aes(size = -log(Entities.pValue , base = 10)))+
  scale_x_discrete(limits=c("PFAS Mixture", "PFDA", "PFHpS", "PFHxS", "PFNA", "PFOA", "PFOS"))+
  scale_size_continuous(range = c(4,12)) +
  labs(x = NULL) + 
  theme_bw() +
  theme( 
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    axis.text.x= element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title=element_text(size=18,face="bold"))

sig_pathway_plot                  




jpeg("~/Projects/BioMe/proteome/output/reactome-pathway/pathwayplot_neg.jpeg",
     units="in", width=30, height=14, res=500)

sig_pathway_plot

dev.off()
             
             
             
## negative directions
bwqs_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/bwqs/bwqs_allpanels_pos.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFAS Mixture")
PFDA_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/PFDA_pos.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFDA")
PFOS_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/PFOS_pos.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFOS")
PFHxS_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/PFHxS_pos.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFHxS")
PFNA_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/PFNA_pos.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFNA")
PFHpS_sig <- read.csv("~/Projects/BioMe/proteome/output/reactome-pathway/PFHpS_pos.csv") %>% filter(Entities.FDR < 0.05) %>% mutate(PFAS = "PFHpS")


sig_pathway<- rbind(bwqs_sig,
                    PFDA_sig,
                    PFOS_sig,
                    PFHxS_sig,
                    PFNA_sig,
                    PFHpS_sig)

sig_pathway$PFAS<- factor(sig_pathway$PFAS,
                          levels = c("PFAS Mixture", "PFDA", "PFHpS", "PFHxS", "PFNA", "PFOA", "PFOS"))

sig_pathway_plot<-  ggplot(sig_pathway, aes(x = PFAS,  y = Pathway.name)) + 
  geom_point(aes(size = -log(Entities.pValue , base = 10)))+
  scale_x_discrete(limits=c("PFAS Mixture", "PFDA", "PFHpS", "PFHxS", "PFNA", "PFOA", "PFOS"))+
  scale_size_continuous(range = c(2,7)) +
  labs(x = NULL) + 
  theme_bw() +
  theme( 
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "bold"),
    axis.text.x= element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18, face = "bold"),
    axis.title=element_text(size=18,face="bold"))

sig_pathway_plot                  




jpeg("~/Projects/BioMe/proteome/output/reactome-pathway/pathwayplot_pos.jpeg",
     units="in", width=30, height=14, res=500)

sig_pathway_plot

dev.off()



