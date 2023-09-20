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


##------------------------------------------- import data
BioMe_proteome_PFAS_wide <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_wide_imputed.txt")
protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")
BioMe_proteome_PFAS_long <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_long.txt")


protein_in_allpanels<- protein_in_panel$OlinkID

# ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_allpanels)), use = 'complete.obs')
# et <- eigen(ght)
# cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.0003990423

#################################
###### volcano plot #############
#################################

## whole samples
#-------------------------------------------- Continuous PFDA - quartile
##------------------------------------------- adjusted
PFDA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFDA_allpanel_adlm_q.csv")
d_lm_pfas_plot <- PFDA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"




top_5_name <- d_lm_pfas_plot %>%
              arrange(p.value) %>% 
              slice_head(n = 6)  %>%
              pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFDA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFDA_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFOA - tertile
##------------------------------------------- adjusted
PFOA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOA_allpanel_adlm_q.csv")
d_lm_pfas_plot <- PFOA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"


d_lm_pfas_plot$Association<- factor(d_lm_pfas_plot$Association,
                                    levels = c("Negative", "Null", "Positive"))


top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFOA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFOA_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFOS - tertile
##------------------------------------------- adjusted
PFOS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOS_allpanel_adlm_q.csv")
d_lm_pfas_plot <- PFOS_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"




top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFOS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))



volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFOS_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFHxS - tertile
##------------------------------------------- adjusted
PFHxS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHxS_allpanel_adlm_q.csv")
d_lm_pfas_plot <- PFHxS_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"



top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHxS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFHxS_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()


#-------------------------------------------- Continuous PFNA - tertile
##------------------------------------------- adjusted
PFNA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFNA_allpanel_adlm_q.csv")
d_lm_pfas_plot <- PFNA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"



top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFNA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFNA_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFHpS - tertile
##------------------------------------------- adjusted
PFHpS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpS_allpanel_adlm_q.csv")
d_lm_pfas_plot <- PFHpS_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"



top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHpS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFHpS_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()





#-------------------------------------------- Continuous PFHpA - binary
##------------------------------------------- adjusted
PFHpA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpA_allpanel_adlm_bi.csv")
d_lm_pfas_plot <- PFHpA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"


top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHpA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFHpA_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



## cases samples
#-------------------------------------------- Continuous PFDA - quartile
##------------------------------------------- adjusted
PFDA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFDA_allpanel_adlm_q_case.csv")
d_lm_pfas_plot <- PFDA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"




top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFDA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFDA_all_adlm_q_case.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFOA - tertile
##------------------------------------------- adjusted
PFOA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOA_allpanel_adlm_q_case.csv")
d_lm_pfas_plot <- PFOA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"


d_lm_pfas_plot$Association<- factor(d_lm_pfas_plot$Association,
                                    levels = c("Negative", "Null", "Positive"))


top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFOA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFOA_all_adlm_q_case.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFOS - tertile
##------------------------------------------- adjusted
PFOS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOS_allpanel_adlm_q_case.csv")
d_lm_pfas_plot <- PFOS_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"




top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFOS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))



volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFOS_all_adlm_q_case.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFHxS - tertile
##------------------------------------------- adjusted
PFHxS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHxS_allpanel_adlm_q_case.csv")
d_lm_pfas_plot <- PFHxS_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"



top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHxS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFHxS_all_adlm_q_case.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()


#-------------------------------------------- Continuous PFNA - tertile
##------------------------------------------- adjusted
PFNA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFNA_allpanel_adlm_q_case.csv")
d_lm_pfas_plot <- PFNA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"



top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFNA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFNA_all_adlm_q_case.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFHpS - tertile
##------------------------------------------- adjusted
PFHpS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpS_allpanel_adlm_q_case.csv")
d_lm_pfas_plot <- PFHpS_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"



top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHpS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFHpS_all_adlm_q_case.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()





#-------------------------------------------- Continuous PFHpA - binary
##------------------------------------------- adjusted
PFHpA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpA_allpanel_adlm_bi_case.csv")
d_lm_pfas_plot <- PFHpA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"


top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHpA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-1, 1)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFHpA_all_adlm_q_case.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()





## control samples
#-------------------------------------------- Continuous PFOS - tertile
##------------------------------------------- adjusted
PFOS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOS_allpanel_adlm_q_control.csv")
d_lm_pfas_plot <- PFOS_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"




top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFOS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))



volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFOS_all_adlm_q_control.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFOA - tertile
##------------------------------------------- adjusted
PFOA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFOA_allpanel_adlm_q_control.csv")
d_lm_pfas_plot <- PFOA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"


d_lm_pfas_plot$Association<- factor(d_lm_pfas_plot$Association,
                                    levels = c("Negative", "Null", "Positive"))


top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFOA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFOA_all_adlm_q_control.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()

#-------------------------------------------- Continuous PFNA - tertile
##------------------------------------------- adjusted
PFNA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFNA_allpanel_adlm_q_control.csv")
d_lm_pfas_plot <- PFNA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"



top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFNA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFNA_all_adlm_q_control.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()




#-------------------------------------------- Continuous PFDA - quartile
##------------------------------------------- adjusted
PFDA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFDA_allpanel_adlm_q_control.csv")
d_lm_pfas_plot <- PFDA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"




top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFDA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFDA_all_adlm_q_control.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFHxS - tertile
##------------------------------------------- adjusted
PFHxS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHxS_allpanel_adlm_q_control.csv")
d_lm_pfas_plot <- PFHxS_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"



top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHxS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFHxS_all_adlm_q_control.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()




#-------------------------------------------- Continuous PFHpS - tertile
##------------------------------------------- adjusted
PFHpS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpS_allpanel_adlm_q_control.csv")
d_lm_pfas_plot <- PFHpS_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"



top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHpS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))



jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFHpS_all_adlm_q_control.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()





#-------------------------------------------- Continuous PFHpA - binary
##------------------------------------------- adjusted
PFHpA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/exwas_PFHpA_allpanel_adlm_bi_control.csv")
d_lm_pfas_plot <- PFHpA_all_adlm_results
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"


top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHpA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-1, 1)+
          theme_bw() +
          scale_color_manual(values=c("Negative" = "blue", 
                                      "Null" = "black", 
                                      "Positive" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/imputed/adjusted/PFHpA_all_adlm_q_control.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()






