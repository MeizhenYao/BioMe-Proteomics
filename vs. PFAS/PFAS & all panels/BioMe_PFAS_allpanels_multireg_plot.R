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
BioMe_proteome_PFAS_wide <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_wide.txt")
protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")
BioMe_proteome_PFAS_long <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_long.txt")


protein_in_allpanels<- protein_in_panel$protein

ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_allpanels)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.0003990423

#################################
###### volcano plot #############
#################################
##------------------------------------------- T2D vs. all panels
PFDA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_T2D_all_adlm_cont.csv")
d_lm_pfas_plot <- PFDA_all_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_allpanels)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: T2D vs. Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/T2D_all_adlm_cont.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Binary PFDA
##------------------------------------------- unadjusted
PFDA_all_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_unlm.csv")
d_lm_pfas_plot <- PFDA_all_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"


d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
              arrange(p.value) %>% 
              slice_head(n = 5)  %>%
              pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFDA (Higher vs. Lower) vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/PFDA_all_unlm.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFDA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_adlm.csv")
d_lm_pfas_plot <- PFDA_all_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
              arrange(p.value) %>% 
              slice_head(n = 5)  %>%
              pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFDA (Higher vs. Lower) vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/PFDA_all_adlm.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()


#-------------------------------------------- Continuous PFDA
##------------------------------------------- unadjusted
PFDA_all_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_unlm_cont.csv")
d_lm_pfas_plot <- PFDA_all_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFDA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/PFDA_all_unlm_cont.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFDA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_adlm_cont.csv")
d_lm_pfas_plot <- PFDA_all_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
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
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/PFDA_all_adlm_cont.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()




#-------------------------------------------- Continuous PFDA - tertile
##------------------------------------------- unadjusted
PFDA_all_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_unlm_q.csv")
d_lm_pfas_plot <- PFDA_all_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFDA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/unadjusted/PFDA_all_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFDA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFDA_all_adlm_q.csv")
d_lm_pfas_plot <- PFDA_all_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
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
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/adjusted/PFDA_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFOA - tertile
##------------------------------------------- unadjusted
PFOA_all_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFOA_all_unlm_q.csv")
d_lm_pfas_plot <- PFOA_all_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFOA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/unadjusted/PFOA_all_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFOA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFOA_all_adlm_q.csv")
d_lm_pfas_plot <- PFOA_all_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
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
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/adjusted/PFOA_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFOS - tertile
##------------------------------------------- unadjusted
PFOS_all_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFOS_all_unlm_q.csv")
d_lm_pfas_plot <- PFOS_all_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFOS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/unadjusted/PFOS_all_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFOS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFOS_all_adlm_q.csv")
d_lm_pfas_plot <- PFOS_all_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"


d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
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
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/adjusted/PFOS_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()


#-------------------------------------------- Continuous PFHpA - tertile
##------------------------------------------- unadjusted
PFHpA_all_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHpA_all_unlm_q.csv")
d_lm_pfas_plot <- PFHpA_all_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"


d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFHpA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/unadjusted/PFHpA_all_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFHpA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHpA_all_adlm_q.csv")
d_lm_pfas_plot <- PFHpA_all_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
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
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/adjusted/PFHpA_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()

#-------------------------------------------- Continuous PFHxS - tertile
##------------------------------------------- unadjusted
PFHxS_all_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHxS_all_unlm_q.csv")
d_lm_pfas_plot <- PFHxS_all_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFHxS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/unadjusted/PFHxS_all_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFHxS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHxS_all_adlm_q.csv")
d_lm_pfas_plot <- PFHxS_all_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"


d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
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
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/adjusted/PFHxS_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()


#-------------------------------------------- Continuous PFNA - tertile
##------------------------------------------- unadjusted
PFNA_all_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFNA_all_unlm_q.csv")
d_lm_pfas_plot <- PFNA_all_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFNA vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/unadjusted/PFNA_all_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFNA_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFNA_all_adlm_q.csv")
d_lm_pfas_plot <- PFNA_all_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
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
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/adjusted/PFNA_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFHpS - tertile
##------------------------------------------- unadjusted
PFHpS_all_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHpS_all_unlm_q.csv")
d_lm_pfas_plot <- PFHpS_all_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFHpS vs. all panels Proteomics") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 8,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/unadjusted/PFHpS_all_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFHpS_all_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/all panels/exwas_PFHpS_all_adlm_q.csv")
d_lm_pfas_plot <- PFHpS_all_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  slice_head(n = 6)  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
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
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/multireg/tertile/adjusted/PFHpS_all_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()









