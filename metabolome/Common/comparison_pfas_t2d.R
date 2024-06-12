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
library(ggrepel)


##------------------------------------------- q.value
proteome_vs_pfas_bwqs <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/scale/proteome_vs_pfas_bwqs.txt")
proteome_t2d_model <- fread("~/Projects/BioMe/proteome/input/pwas/proteome_vs_t2d_scale.txt")


proteome_pfas<- (proteome_vs_pfas_bwqs %>% filter(q.value < 0.1))$Protein_name
proteome_t2d<- (proteome_t2d_model %>% filter(q.value < 0.1))$Protein_name
common<- intersect(proteome_pfas, proteome_t2d)



proteome_pfas_value<- proteome_vs_pfas_bwqs %>% 
                      filter(Protein_name %in% common) %>% 
                      pull(mean)


proteome_t2d_value<-  proteome_t2d_model %>% 
                      filter(Protein_name %in% common) %>% 
                      pull(Value)


check<- data.frame(cbind(common,round(proteome_pfas_value,2), round(proteome_t2d_value,2)))
colnames(check)[1]<- "protein"
colnames(check)[2]<- "proteome_pfas_value"
colnames(check)[3]<- "proteome_t2d_value"



##------------------------------------------- p.value

proteome_pfas<- (proteome_vs_pfas_bwqs %>% filter(p.value < 0.05))$Protein_name
proteome_t2d<- (proteome_t2d_model %>% filter(p.value < 0.05))$Protein_name
common<- intersect(proteome_pfas, proteome_t2d)



proteome_pfas_value<- proteome_vs_pfas_bwqs %>% 
  filter(Protein_name %in% common) %>% 
  pull(mean)


proteome_t2d_value<-  proteome_t2d_model %>% 
  filter(Protein_name %in% common) %>% 
  pull(Value)


check<- data.frame(cbind(common,round(proteome_pfas_value,2), round(proteome_t2d_value,2)))
colnames(check)[1]<- "protein"
colnames(check)[2]<- "proteome_pfas_value"
colnames(check)[3]<- "proteome_t2d_value"



##------------------------------------------- dataset

proteome_pfas<- (proteome_vs_pfas_bwqs %>% filter(p.value < 0.05))$Protein_name
proteome_t2d<- (proteome_t2d_model %>% filter(p.value < 0.05))$Protein_name
common<- intersect(proteome_pfas, proteome_t2d)



proteome_pfas_value<- proteome_vs_pfas_bwqs %>% 
  filter(Protein_name %in% common) %>% 
  rename_with( ~ paste0("pfas_", .x))


proteome_t2d_value<-  proteome_t2d_model %>% 
  filter(Protein_name %in% common) %>% 
  rename_with( ~ paste0("t2d_", .x))


pfas_t2d_value<- data.frame(proteome_pfas_value,
                            proteome_t2d_value)


write.csv(pfas_t2d_value, "~/Projects/BioMe/proteome/input/comparison/pfas_t2d_pvalue.csv", row.names = FALSE)








