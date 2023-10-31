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
proteome_t2d_model <- fread("~/Projects/BioMe/proteome/input/pwas/proteome_vs_t2d_scale.txt")

protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")


protein_in_allpanels<- protein_in_panel$OlinkID



#################################
###### volcano plot #############
#################################

## whole samples

d_lm_pfas_plot <- proteome_t2d_model
cutoff <- max(d_lm_pfas_plot[d_lm_pfas_plot$q.value<0.05]$p.value)
cut_label<- 0.05
d_lm_pfas_plot$Association <- "FDR > 0.05 & p.value > 0.05"
d_lm_pfas_plot$Association[d_lm_pfas_plot$p.value < cut_label] <- "FDR > 0.05 & p.value < 0.05"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value > 0 & d_lm_pfas_plot$q.value < cut_label] <- "FDR < 0.05 (Positive)"
d_lm_pfas_plot$Association[d_lm_pfas_plot$Value < 0 & d_lm_pfas_plot$q.value < cut_label] <- "FDR < 0.05 (Negative)"




sig_name <- d_lm_pfas_plot %>%
  arrange(p.value) %>%
  filter(q.value < 0.05) %>%
  pull(OlinkID)






vol <- (ggplot(d_lm_pfas_plot, aes(x=Value, y=-log10(p.value), col=Association)) +# Show all points
          geom_point(size=2) +
          geom_hline(yintercept= -log(cut_label, base = 10), color = "black", size = 1) + 
          geom_hline(yintercept= -log(cutoff, base = 10), color = "black", size = 1, linetype = "dashed") + 
          labs(x = "Beta Coefficients", title = "") +
          geom_label_repel(data = subset(d_lm_pfas_plot, 
                                         OlinkID %in% c(sig_name)),
                           aes(label = Protein_name),
                           size = 9,
                           box.padding = unit(0.5, "lines"),
                           max.overlaps = getOption("ggrepel.max.overlaps", default = 60),
                           force = 2, force_pull = 2, show.legend = FALSE) + 
          xlim(-1, 1)+
          theme_bw() + 
          scale_color_manual(values=c("FDR > 0.05 & p.value > 0.05" = "grey",
                                      "FDR > 0.05 & p.value < 0.05" = "black",
                                      "FDR < 0.05 (Negative)" = "blue", 
                                      "FDR < 0.05 (Positive)" = "red")))


volcano_pos_pfas_met <- vol + theme(legend.position = "bottom",
                                    legend.text = element_text(size = 22, face = "bold"),
                                    legend.title = element_text(size = 22, face = "bold"),
                                    plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 22, face = "bold"),
                                    axis.text.y = element_text(size = 22, face = "bold"),
                                    axis.title=element_text(size=22,face="bold"))

volcano_pos_pfas_met

length((d_lm_pfas_plot %>% filter(q.value<0.05 & Value > 0))$p.value)
length((d_lm_pfas_plot %>% filter(q.value<0.05 & Value < 0))$p.value)



jpeg("~/Projects/BioMe/proteome/output/Proteome vs. T2D/proteome_t2d.jpeg",
     units="in", width=20, height=15, res=500)

volcano_pos_pfas_met

dev.off()

