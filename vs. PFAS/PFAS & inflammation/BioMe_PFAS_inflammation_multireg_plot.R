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


BioMe_proteome_PFAS_wide <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_wide_norm_imputed_v2.txt")
protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")


protein_in_inflammation<- (protein_in_panel %>% filter(Panel == "Inflammation"))$OlinkID



#################################
###### volcano plot #############
#################################
##------------------------------------------- T2D vs. inflammation

PFDA_inflam_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_T2D_inflam_adlm_cont.csv")


d_lm_pfas_plot <- PFDA_inflam_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
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
          labs(x = "Beta Coefficients", title = "Adjusted regression: T2D vs. Inflammation Proteome") +
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

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/T2D_inflammation_adlm_cont.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Binary PFDA
##------------------------------------------- unadjusted
PFDA_inflam_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_unlm.csv")
d_lm_pfas_plot <- PFDA_inflam_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

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
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFDA (Higher vs. Lower) vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/PFDA_inflammation_unlm.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFDA_inflam_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_adlm.csv")
d_lm_pfas_plot <- PFDA_inflam_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

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
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFDA (Higher vs. Lower) vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/PFDA_inflammation_adlm.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()


#-------------------------------------------- Continuous PFDA
##------------------------------------------- unadjusted
PFDA_inflam_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_unlm_cont.csv")
d_lm_pfas_plot <- PFDA_inflam_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

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
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFDA vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/PFDA_inflammation_unlm_cont.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFDA_inflam_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_adlm_cont.csv")
d_lm_pfas_plot <- PFDA_inflam_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
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
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFDA vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/PFDA_inflammation_adlm_cont.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()




#-------------------------------------------- Continuous PFDA - tertile
##------------------------------------------- unadjusted
PFDA_inflam_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_unlm_q.csv")
d_lm_pfas_plot <- PFDA_inflam_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

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
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFDA vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/unadjusted/PFDA_inflammation_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFDA_inflam_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFDA_inflam_adlm_q_imputed1.csv")
d_lm_pfas_plot <- PFDA_inflam_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986
# [1] 0.004541536
# [1] 0.004537802

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
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFDA vs. Inflammation Proteome") +
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

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/adjusted/PFDA_inflammation_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFOA - tertile
##------------------------------------------- unadjusted
PFOA_inflam_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFOA_inflam_unlm_q.csv")
d_lm_pfas_plot <- PFOA_inflam_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

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
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFOA vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/unadjusted/PFOA_inflammation_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFOA_inflam_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFOA_inflam_adlm_q.csv")
d_lm_pfas_plot <- PFOA_inflam_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
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
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFOA vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/adjusted/PFOA_inflammation_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFOS - tertile
##------------------------------------------- unadjusted
PFOS_inflam_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFOS_inflam_unlm_q.csv")
d_lm_pfas_plot <- PFOS_inflam_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

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
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFOS vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/unadjusted/PFOS_inflammation_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFOS_inflam_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFOS_inflam_adlm_q.csv")
d_lm_pfas_plot <- PFOS_inflam_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
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
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFOS vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/adjusted/PFOS_inflammation_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()


#-------------------------------------------- Continuous PFHpA - tertile
##------------------------------------------- unadjusted
PFHpA_inflam_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFHpA_inflam_unlm_q.csv")
d_lm_pfas_plot <- PFHpA_inflam_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

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
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFHpA vs. Inflammation Proteome") +
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

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/unadjusted/PFHpA_inflammation_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFHpA_inflam_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFHpA_inflam_adlm_q.csv")
d_lm_pfas_plot <- PFHpA_inflam_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

d_lm_pfas_plot$delabel <- NA
d_lm_pfas_plot$delabel[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label ] <- d_lm_pfas_plot$Protein_name[d_lm_pfas_plot$Association != "Null" & d_lm_pfas_plot$p.value < cut_label]

top_5_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>% 
  filter(Protein_name=="Tumor necrosis factor")  %>%
  pull(Protein_name)


vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1) + 
          # geom_hline(yintercept= -log(cut_label, base = 10), color = "green", size = 1) + 
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHpA vs. Inflammation Proteome") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         Protein_name %in% top_5_name),
                           aes(label = Protein_name),
                           size = 12,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-0.7, 0.7)+
          theme_bw() + 
          scale_color_manual(values=c("blue", "black", "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "none",
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.title.x = element_text(size = 20, face = "bold"),
                                    axis.title.y = element_text(size = 20, face = "bold"),
                                    axis.text.x= element_text(size = 18, face = "bold"),
                                    axis.text.y = element_text(size = 18, face = "bold"))






jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/adjusted/PFHpA_inflammation_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()

#-------------------------------------------- Continuous PFHxS - tertile
##------------------------------------------- unadjusted
PFHxS_inflam_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFHxS_inflam_unlm_q.csv")
d_lm_pfas_plot <- PFHxS_inflam_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

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
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFHxS vs. Inflammation Proteome") +
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

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/unadjusted/PFHxS_inflammation_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFHxS_inflam_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFHxS_inflam_adlm_q.csv")
d_lm_pfas_plot <- PFHxS_inflam_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
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
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHxS vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/adjusted/PFHxS_inflammation_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()


#-------------------------------------------- Continuous PFNA - tertile
##------------------------------------------- unadjusted
PFNA_inflam_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFNA_inflam_unlm_q.csv")
d_lm_pfas_plot <- PFNA_inflam_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

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
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFNA vs. Inflammation Proteome") +
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

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/unadjusted/PFNA_inflammation_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFNA_inflam_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFNA_inflam_adlm_q.csv")
d_lm_pfas_plot <- PFNA_inflam_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
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
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFNA vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/adjusted/PFNA_inflammation_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#-------------------------------------------- Continuous PFHpS - tertile
##------------------------------------------- unadjusted
PFHpS_inflam_unlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFHpS_inflam_unlm_q.csv")
d_lm_pfas_plot <- PFHpS_inflam_unlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
et <- eigen(ght)
cut_label<- 1/ (sum((et$values>1 + 0)* (et$values - 1)))
# [1] 0.004458986

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
          labs(x = "Beta Coefficients", title = "Unadjusted regression: PFHpS vs. Inflammation Proteome") +
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

jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/unadjusted/PFHpS_inflammation_unlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



##------------------------------------------- adjusted
PFHpS_inflam_adlm_results <- read.csv("~/Projects/BioMe/proteome/input/exwas/exwas_PFHpS_inflam_adlm_q.csv")
d_lm_pfas_plot <- PFHpS_inflam_adlm_results
cutoff <- 0.05
d_lm_pfas_plot$Association <- "Null"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$p.value < cutoff] <- "Positive"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$p.value < cutoff] <- "Negative"



ght <- cor(BioMe_proteome_PFAS_wide %>% select(all_of(protein_in_inflammation)), use = 'complete.obs')
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
          labs(x = "Beta Coefficients", title = "Adjusted regression: PFHpS vs. Inflammation Proteome") +
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
jpeg("~/Projects/BioMe/proteome/output/PFAS vs. inflammation/multireg/tertile/adjusted/PFHpS_inflammation_adlm_q.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()









