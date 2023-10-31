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
proteome_vs_pfas_bwqs <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/proteome_vs_pfas_bwqs.txt")
bwqs_pfas_weight <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/bwqs_pfas_weight.txt")
proteome_vs_pfas_bwqs_case <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case/proteome_vs_pfas_bwqs_case.txt")
bwqs_pfas_weight_case <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case/bwqs_pfas_weight_case.txt")
proteome_vs_pfas_bwqs_control <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control/proteome_vs_pfas_bwqs_control.txt")
bwqs_pfas_weight_control <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control/bwqs_pfas_weight_control.txt")
protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")


protein_in_allpanels<- protein_in_panel$OlinkID




#################################
###### volcano plot #############
#################################

## whole samples

d_lm_pfas_plot <- proteome_vs_pfas_bwqs
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$mean > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$mean < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"




top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 15)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=mean, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "BWQS: PFAS mixture vs. all panels Proteomics") +
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


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/bwqs/bwqs_all_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()


## case samples

d_lm_pfas_plot <- proteome_vs_pfas_bwqs_case
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$mean > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$mean < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"




top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=mean, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "BWQS: PFAS mixture vs. all panels Proteomics") +
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

length((d_lm_pfas_plot %>% filter(p.value<0.05))$p.value)

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/bwqs/bwqs_all_q_case.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()






## case samples

d_lm_pfas_plot <- proteome_vs_pfas_bwqs_control
cutoff <- 0.05
cut_label<- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$mean > 0 & d_lm_pfas_plot$q.value < cut_label] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$mean < 0 & d_lm_pfas_plot$q.value < cut_label] <- "Negative"




top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=mean, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          labs(x = "Beta Coefficients", title = "BWQS: PFAS mixture vs. all panels Proteomics") +
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

length((d_lm_pfas_plot %>% filter(p.value<0.05))$p.value)
length((d_lm_pfas_plot %>% filter(q.value<0.05))$p.value)
length((d_lm_pfas_plot %>% filter(q.value<0.1))$p.value)
length((d_lm_pfas_plot %>% filter(q.value<0.2))$p.value)

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/bwqs/bwqs_all_q_control.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()







