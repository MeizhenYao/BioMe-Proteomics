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
proteome_vs_pfas_bwqs <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/scale/proteome_vs_pfas_bwqs.txt")
bwqs_pfas_weight <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/scale/bwqs_pfas_weight.txt")
proteome_vs_pfas_bwqs_case <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/proteome_vs_pfas_bwqs_case.txt")
bwqs_pfas_weight_case <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/case_scale/bwqs_pfas_weight_case.txt")
proteome_vs_pfas_bwqs_control <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/proteome_vs_pfas_bwqs_control.txt")
bwqs_pfas_weight_control <- fread("~/Projects/BioMe/proteome/input/exwas/all panels/batch_imputed/bwqs/control_scale/bwqs_pfas_weight_control.txt")
protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")



protein_in_allpanels<- protein_in_panel$OlinkID

bwqs_pfas_weight_long<- bwqs_pfas_weight %>% 
                        pivot_longer(
                          cols = c("PFDA", "PFHxS", "PFHpS", "PFNA", "PFOA", "PFOS"),
                          names_to = "PFAS",
                          values_to = "weight"
                        )

bwqs_pfas_weight_case_long<-  bwqs_pfas_weight_case %>% 
                              pivot_longer(
                                cols = c("PFDA", "PFHxS", "PFHpS", "PFNA", "PFOA", "PFOS"),
                                names_to = "PFAS",
                                values_to = "weight"
                              )

bwqs_pfas_weight_control_long<- bwqs_pfas_weight_control %>% 
                                pivot_longer(
                                  cols = c("PFDA", "PFHxS", "PFHpS", "PFNA", "PFOA", "PFOS"),
                                  names_to = "PFAS",
                                  values_to = "weight"
                                )

#################################
###### heatmap plot #############
#################################

## whole samples

vol <- (ggplot(bwqs_pfas_weight_long, aes(x=PFAS, y=Gene_name)) +
          geom_tile(aes(fill = weight)) +
          scale_fill_distiller(type = "seq",
                               palette = 8,
                               direction = 1,
                               values = NULL,
                               space = "Lab",
                               na.value = "grey50",
                               guide = "colourbar",
                               aesthetics = "fill") +
          labs(title = "",
               y = "Proteins",
               x = ""))


volcano_pos_pfas_met <- vol + theme(plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 22, face = "bold"),
                                    axis.text.y = element_blank(),
                                    legend.text = element_text(size = 22, face = "bold"),
                                    legend.title = element_text(size = 22, face = "bold"),
                                    axis.title=element_text(size=22,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/bwqs/scale/bwqs_all_weight_heatmap.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()

## case samples

vol <- (ggplot(bwqs_pfas_weight_case_long, aes(x=PFAS, y=Gene_name)) +
          geom_tile(aes(fill = weight)) +
          scale_fill_distiller(type = "seq",
                               palette = 8,
                               direction = 1,
                               values = NULL,
                               space = "Lab",
                               na.value = "grey50",
                               guide = "colourbar",
                               aesthetics = "fill") +
          labs(title = "Estimated Weight from BWQS",
               y = "Protein",
               x = "PFAS"))


volcano_pos_pfas_met <- vol + theme(plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_blank(),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/bwqs/bwqs_all_weight_heatmap_case.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()

## control samples

vol <- (ggplot(bwqs_pfas_weight_control_long, aes(x=PFAS, y=Gene_name)) +
          geom_tile(aes(fill = weight)) +
          scale_fill_distiller(type = "seq",
                               palette = 8,
                               direction = 1,
                               values = NULL,
                               space = "Lab",
                               na.value = "grey50",
                               guide = "colourbar",
                               aesthetics = "fill") +
          labs(title = "Estimated Weight from BWQS",
               y = "Protein",
               x = "PFAS"))


volcano_pos_pfas_met <- vol + theme(plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_blank(),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/bwqs/bwqs_all_weight_heatmap_control.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



#################################
###### boxplot #############
#################################

## whole samples

vol <- (ggplot(bwqs_pfas_weight_long, aes(x=PFAS, y=weight)) +
          geom_boxplot(fill = "dark green") +
          labs(title = "",
               y = "Estimated Weight from BWQS",
               x = "PFAS"))


volcano_pos_pfas_met <- vol + theme(plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/bwqs/bwqs_all_weight_boxplot.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()



## CASE samples

vol <- (ggplot(bwqs_pfas_weight_case_long, aes(x=PFAS, y=weight)) +
          geom_boxplot(fill = "dark green") +
          labs(title = "",
               y = "Estimated Weight from BWQS",
               x = "PFAS"))


volcano_pos_pfas_met <- vol + theme(plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/bwqs/bwqs_all_weight_boxplot_case.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()


## Control samples

vol <- (ggplot(bwqs_pfas_weight_control_long, aes(x=PFAS, y=weight)) +
          geom_boxplot(fill = "dark green") +
          labs(title = "",
               y = "Estimated Weight from BWQS",
               x = "PFAS"))


volcano_pos_pfas_met <- vol + theme(plot.title = element_text(size = 24, face = "bold"),
                                    axis.text.x= element_text(size = 14, face = "bold"),
                                    axis.text.y = element_text(size = 14, face = "bold"),
                                    axis.title=element_text(size=14,face="bold"))


jpeg("~/Projects/BioMe/proteome/output/PFAS vs. all panels/bwqs/bwqs_all_weight_boxplot_control.jpeg",
     units="in", width=16, height=12, res=500)

volcano_pos_pfas_met

dev.off()






