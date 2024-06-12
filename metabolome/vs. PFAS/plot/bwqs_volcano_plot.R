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
proteome_vs_pfas_bwqs_neg <- fread("~/Projects/BioMe/proteome/input/exwas/metabolome/met_vs_pfas_bwqs_plot_neg.txt")
bwqs_pfas_weight_neg <- fread("~/Projects/BioMe/proteome/input/exwas/metabolome/met_pfas_weight_plot_neg.txt")

proteome_vs_pfas_bwqs_pos <- fread("~/Projects/BioMe/proteome/input/exwas/metabolome/met_vs_pfas_bwqs_plot_pos.txt")
bwqs_pfas_weight_pos <- fread("~/Projects/BioMe/proteome/input/exwas/metabolome/met_pfas_weight_plot_pos.txt")



#################################
###### volcano plot #############
#################################

## negative

d_lm_pfas_plot <- proteome_vs_pfas_bwqs_neg
cutoff <- max(d_lm_pfas_plot[d_lm_pfas_plot$q.value<0.1]$p.value)
cut_label<- 0.05
d_lm_pfas_plot$Association <- "FDR > 0.1 & p.value > 0.05"
d_lm_pfas_plot$Association[d_lm_pfas_plot$p.value < cut_label] <- "FDR > 0.1 & p.value < 0.05"
d_lm_pfas_plot$Association[d_lm_pfas_plot$mean > 0 & d_lm_pfas_plot$p.value <= cutoff] <- "FDR < 0.1 (Positive)"
d_lm_pfas_plot$Association[d_lm_pfas_plot$mean < 0 & d_lm_pfas_plot$p.value <= cutoff] <- "FDR < 0.1 (Negative)"


top_5_name <- d_lm_pfas_plot%>%
  filter(mean < 0) %>% 
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(refmet_name)


top_3_pos<- d_lm_pfas_plot %>%
  filter(mean > 0) %>% 
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(refmet_name)






vol <- (ggplot(d_lm_pfas_plot, aes(x=mean, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cut_label, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1, linetype = "dashed") + 
          labs(x = "Beta Coefficients", title = "") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         refmet_name %in% c(top_5_name, top_3_pos)),
                           aes(label = refmet_name),
                           size = 9,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 60),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-1, 1)+
          theme_bw() + 
          scale_color_manual(values=c("FDR > 0.1 & p.value > 0.05" = "grey",
                                      "FDR > 0.1 & p.value < 0.05" = "black",
                                      "FDR < 0.1 (Negative)" = "blue", 
                                      "FDR < 0.1 (Positive)" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "bottom",
                                    legend.text = element_text(size = 22, face = "bold"),
                                    legend.title = element_text(size = 22, face = "bold"),
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 22, face = "bold"),
                                    axis.text.y = element_text(size = 22, face = "bold"),
                                    axis.title=element_text(size=22,face="bold"))

volcano_pos_pfas_met

length((d_lm_pfas_plot %>% filter(q.value<0.1 & mean > 0))$p.value)
length((d_lm_pfas_plot %>% filter(q.value<0.1 & mean < 0))$p.value)



jpeg("~/Projects/BioMe/proteome/output/PFAS vs. metabolome/pfas_metabolome_neg.jpeg",
     units="in", width=20, height=15, res=500)

volcano_pos_pfas_met

dev.off()


## positive

d_lm_pfas_plot <- proteome_vs_pfas_bwqs_pos
cutoff <- max(d_lm_pfas_plot[d_lm_pfas_plot$q.value<0.1]$p.value)
cut_label<- 0.05
d_lm_pfas_plot$Association <- "FDR > 0.1 & p.value > 0.05"
d_lm_pfas_plot$Association[d_lm_pfas_plot$p.value < cut_label] <- "FDR > 0.1 & p.value < 0.05"
d_lm_pfas_plot$Association[d_lm_pfas_plot$mean > 0 & d_lm_pfas_plot$p.value <= cutoff] <- "FDR < 0.1 (Positive)"
d_lm_pfas_plot$Association[d_lm_pfas_plot$mean < 0 & d_lm_pfas_plot$p.value <= cutoff] <- "FDR < 0.1 (Negative)"


top_5_name <- d_lm_pfas_plot%>%
  filter(mean < 0) %>% 
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(refmet_name)


top_3_pos<- d_lm_pfas_plot %>%
  filter(mean > 0) %>% 
  arrange(p.value) %>% 
  slice_head(n = 5)  %>%
  pull(refmet_name)






vol <- (ggplot(d_lm_pfas_plot, aes(x=mean, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cut_label, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1, linetype = "dashed") + 
          labs(x = "Beta Coefficients", title = "") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         refmet_name %in% c(top_5_name, top_3_pos)),
                           aes(label = refmet_name),
                           size = 9,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 60),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-1, 1)+
          theme_bw() + 
          scale_color_manual(values=c("FDR > 0.1 & p.value > 0.05" = "grey",
                                      "FDR > 0.1 & p.value < 0.05" = "black",
                                      "FDR < 0.1 (Negative)" = "blue", 
                                      "FDR < 0.1 (Positive)" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "bottom",
                                    legend.text = element_text(size = 22, face = "bold"),
                                    legend.title = element_text(size = 22, face = "bold"),
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 22, face = "bold"),
                                    axis.text.y = element_text(size = 22, face = "bold"),
                                    axis.title=element_text(size=22,face="bold"))

volcano_pos_pfas_met

length((d_lm_pfas_plot %>% filter(q.value<0.1 & mean > 0))$p.value)
length((d_lm_pfas_plot %>% filter(q.value<0.1 & mean < 0))$p.value)



jpeg("~/Projects/BioMe/proteome/output/PFAS vs. metabolome/pfas_metabolome_pos.jpeg",
     units="in", width=20, height=15, res=500)

volcano_pos_pfas_met

dev.off()










