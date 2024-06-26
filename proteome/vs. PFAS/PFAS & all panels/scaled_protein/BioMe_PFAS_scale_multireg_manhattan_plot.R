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
library(ggrepel)



##------------------------------------------- import data

PFDA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/individual/exwas_PFAS_allpanel_adlm_q.csv") %>% filter (PFAS == "PFDA_Aug21_q") 
PFOA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/individual/exwas_PFAS_allpanel_adlm_q.csv") %>% filter (PFAS == "PFOA_Aug21_q")
PFOS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/individual/exwas_PFAS_allpanel_adlm_q.csv") %>% filter (PFAS == "PFOS_Aug21_q")
PFHxS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/individual/exwas_PFAS_allpanel_adlm_q.csv") %>% filter (PFAS == "PFHxS_Aug21_q")
PFNA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/individual/exwas_PFAS_allpanel_adlm_q.csv") %>% filter (PFAS == "PFNA_Aug21_q")
PFHpS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/individual/exwas_PFAS_allpanel_adlm_q.csv") %>% filter (PFAS == "PFHpS_Aug21_q")



length((PFHpS_all_adlm_results %>% filter(fdr<0.1 & Value > 0))$p.value)
length((PFHpS_all_adlm_results %>% filter(fdr<0.1 & Value < 0))$p.value)



PFDA<-  PFDA_all_adlm_results%>%
        mutate(sig = paste(Protein_name, "(", PFAS, ")"))%>%
        arrange(p.value) %>% 
        slice_head(n = 3)  %>%
        pull(sig)

PFOA<-  PFOA_all_adlm_results%>%
  mutate(sig = paste(Protein_name, "(", PFAS, ")"))%>%
  arrange(p.value) %>% 
  slice_head(n = 1)  %>%
  pull(sig)

PFOS<-  PFOS_all_adlm_results%>%
  mutate(sig = paste(Protein_name, "(", PFAS, ")"))%>%
  arrange(p.value) %>% 
  slice_head(n = 3)  %>%
  pull(sig)

PFHxS<-  PFHxS_all_adlm_results%>%
  mutate(sig = paste(Protein_name, "(", PFAS, ")"))%>%
  arrange(p.value) %>% 
  slice_head(n = 3)  %>%
  pull(sig)

PFNA<-  PFNA_all_adlm_results%>%
  mutate(sig = paste(Protein_name, "(", PFAS, ")"))%>%
  arrange(p.value) %>% 
  slice_head(n = 3)  %>%
  pull(sig)

PFHpS<-  PFHpS_all_adlm_results%>%
  mutate(sig = paste(Protein_name, "(", PFAS, ")"))%>%
  arrange(p.value) %>% 
  slice_head(n = 3)  %>%
  pull(sig)


## combine together
individual_pfas<- rbind(PFDA_all_adlm_results,
                        PFOA_all_adlm_results,
                        PFOS_all_adlm_results,
                        PFHxS_all_adlm_results,
                        PFNA_all_adlm_results,
                        PFHpS_all_adlm_results)%>%
                  mutate(sig = paste(Protein_name, "(", PFAS, ")")) 

# 
# top_15_name<- individual_pfas %>%
#               arrange(p.value) %>% 
#               slice_head(n = 20)  %>%
#               pull(sig)


# pval <- individual_pfas$p.value
# q <- qvalue::qvalue(pval, lambda = 0)
# individual_pfas$fdr.whole <-  q$qvalues
# 
# check<- individual_pfas %>%  
#         filter(fdr.whole < 0.05) %>% 
#         arrange(desc(p.value))%>% 
#         slice_head(n = 1)  %>%
#         pull(p.value)
           
individual_pfas$Association <- "FDR > 0.05 & p.value > 0.05"
individual_pfas$Association[individual_pfas$p.value < 0.05] <- "FDR > 0.05 & p.value < 0.05"
individual_pfas$Association[individual_pfas$Value > 0 & individual_pfas$fdr < 0.05] <- "FDR < 0.05 (Positive)"
individual_pfas$Association[individual_pfas$Value < 0 & individual_pfas$fdr < 0.05] <- "FDR < 0.05 (Negative)"


length((individual_pfas %>% filter(fdr<0.05 & Value > 0))$p.value)
length((individual_pfas %>% filter(fdr<0.05 & Value < 0))$p.value)


manhplot <- ggplot(individual_pfas, aes(x = PFAS, y = -log10(p.value),
                                        color = Association)) +
            geom_hline(yintercept = -log(0.05, base = 10), linewidth = 1.2,color = "grey40", linetype = "solid") + 
            geom_label_repel(data = subset(individual_pfas, 
                                           sig %in% c(PFDA, PFOA, PFOS, PFHxS, PFHpS, PFNA)),
                             size = 6,
                             aes(label = Protein_name), position = position_jitter(width = 0.4, seed = 3215),
                             max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                             segment.color = 'transparent') + 
            geom_point(aes(size = -log10(p.value)), position = position_jitter(width = 0.4, seed = 3215)) +
            scale_color_manual(values=c("FDR > 0.05 & p.value > 0.05" = "grey",
                                        "FDR > 0.05 & p.value < 0.05" = "black",
                                        "FDR < 0.05 (Negative)" = "blue", 
                                        "FDR < 0.05 (Positive)" = "red")) +
            scale_size_continuous(range = c(0.5,3)) +
            labs(x = NULL) + 
            theme_classic() +
            theme( 
              legend.title = element_text(size = 24, face = "bold"),
              legend.text = element_text(size = 24, face = "bold"),
              plot.title = element_text(size = 24, face = "bold"),
              axis.text.x= element_text(size = 24, face = "bold"),
              axis.text.y = element_text(size = 24, face = "bold"),
              axis.title=element_text(size=24,face="bold"))

manhplot



jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/manhat/manhplot_0.05.jpeg",
     units="in", width=26, height=16, res=500)

manhplot

dev.off()






