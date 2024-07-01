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
proteome_vs_pfas_bwqs <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/nonimputed/proteome_vs_pfas_bwqs.txt")
proteome_t2d_model <- fread("~/Projects/BioMe/proteome/input/pwas/nonimputed/proteome_vs_t2d_scale.txt")
protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")



proteome_pfas<- (proteome_vs_pfas_bwqs %>% filter(q.value < 0.05))$Protein_name
proteome_t2d<- (proteome_t2d_model %>% filter(q.value < 0.05))$Protein_name
common<- intersect(proteome_pfas, proteome_t2d)



proteome_pfas_value<- proteome_vs_pfas_bwqs %>% 
                      filter(Protein_name %in% common) %>% 
                      pull(mean)


proteome_t2d_value<-  proteome_t2d_model %>% 
                      filter(Protein_name %in% common) %>% 
                      pull(Value)


check_mixture<- data.frame(cbind(common,round(proteome_pfas_value,2), round(proteome_t2d_value,2)))
colnames(check_mixture)[1]<- "Protein_name"
colnames(check_mixture)[2]<- "proteome_pfas_value"
colnames(check_mixture)[3]<- "proteome_t2d_value"

check_mixture_samedir<- check_mixture %>% 
                        filter((proteome_pfas_value < 0 & proteome_t2d_value < 0) | (proteome_pfas_value > 0 & proteome_t2d_value > 0))%>% 
                        inner_join(assay_list, by = "Protein_name")




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


write.csv(pfas_t2d_value, "~/Projects/BioMe/proteome/input/comparison/nonimputed/pfas_t2d_pvalue.csv", row.names = FALSE)



##------------------------------------------- PFAS: Holm-Bonferroni theorem or FDR
# cutoff<- 0.05/length(proteome_vs_pfas_bwqs$p.value)


proteome_vs_pfas_bwqs_sig<- proteome_vs_pfas_bwqs %>% 
                            filter(q.value < 0.05) %>% 
                            inner_join(protein_in_panel[, c("Protein_name", "Panel")], by = "Protein_name") %>% 
                            select(Panel, Protein_name, Gene_name, mean, p.value, q.value)



proteome_vs_pfas_bwqs_sig



write_xlsx(proteome_vs_pfas_bwqs_sig, "~/Projects/BioMe/proteome/output/PFAS vs. all panels/bwqs/scale/proteome_vs_pfas_bwqs_sig.xlsx")






