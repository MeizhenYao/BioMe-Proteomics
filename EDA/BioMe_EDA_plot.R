
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
library(naniar)
library(reactable)
library(sandwich)
library(ggrepel)
library(mice)


##------------------------------------------- import data
BioMe_proteome_PFAS_wide <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_wide.txt")
protein_in_panel <- fread("~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt")
BioMe_proteome_PFAS_long <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_long.txt")


protein_in_allpanels<- protein_in_panel$protein







###################################### PFAS
#--- Table
plasma_PFAS<- c('PFOS_Aug21',	'PFOA_Aug21',	'PFNA_Aug21',	"PFDA_Aug21", "PFHpS_Aug21", "PFHxS_Aug21", "PFHpA_Aug21")
LOD_value<- c(0.01,	0.01, 0.01,	0.01, 0.01, 0.01, 0.01)


statistics<- c('LOD value',	'% > LOD',	'Min',	'Percentile 25', 'Median', 'Percentile 75', 'Max', 'Missing')




plasma_PFAS_baseline_info<- data.frame(matrix(ncol = 7, nrow = 7))

for (i in 1:7){
  BioMe_proteome_PFAS_wide_plasma_PFAS<- BioMe_proteome_PFAS_wide[ ,..plasma_PFAS]
  
  plasma_PFAS_baseline_info[i, ]<- data.frame (as.character((BioMe_proteome_PFAS_wide_plasma_PFAS %>% 
                                                 filter(BioMe_proteome_PFAS_wide_plasma_PFAS[,..i]>0.01 | BioMe_proteome_PFAS_wide_plasma_PFAS[,..i]==0.01) %>% 
                                                 dplyr::summarise(percent=percent(n()/nrow(BioMe_proteome_PFAS_wide_plasma_PFAS %>% na.omit()),1)))$percent), # %>LOD
                                              format(round(min(BioMe_proteome_PFAS_wide_plasma_PFAS[ ,..i], na.rm=TRUE),3),nsmall = 3), #MIN
                                              format( round(quantile(BioMe_proteome_PFAS_wide_plasma_PFAS[ ,..i],c(.25),na.rm=TRUE),3),nsmall = 3), #P25
                                              format(round(quantile(BioMe_proteome_PFAS_wide_plasma_PFAS[ ,..i],c(.50),na.rm=TRUE),3),nsmall = 3), # MEDIAN
                                              format(round(quantile(BioMe_proteome_PFAS_wide_plasma_PFAS[ ,..i],c(.75),na.rm=TRUE),3),nsmall = 3), #P75
                                              format(round(max(BioMe_proteome_PFAS_wide_plasma_PFAS[ ,..i], na.rm=TRUE),3),nsmall = 3), # MAX
                                              (BioMe_proteome_PFAS_wide_plasma_PFAS %>% 
                                                 filter(is.na(BioMe_proteome_PFAS_wide_plasma_PFAS[ ,..i])==TRUE) %>% 
                                                 dplyr::summarise(n=n()))$n) # MISSING NUMBER
}

plasma_PFAS_baseline_info_name<- cbind(statistics,
                                      as.data.frame(t(data.frame(LOD_value,
                                                                 plasma_PFAS_baseline_info))))


plasma_PFAS_baseline_info_name


plasma_POP_info_name_table<- flextable(plasma_PFAS_baseline_info_name) %>% 
  set_header_labels(V1='PFOS_Aug21',	V2='PFOA_Aug21',	V3='PFNA_Aug21',	V4="PFDA_Aug21", V5="PFHpS_Aug21", V6="PFHxS_Aug21", V7="PFHpA_Aug21") %>% 
  theme_box() 

plasma_POP_info_name_table


#--- corrplot
corr_fun<- function(data){
  
  corrplot(data, 
           method="color" ,
           col = colorRampPalette(c("steelblue", "white", "darkred"))(100),cl.lim=c(0,1),
           type="upper",
           # order="hclust" ,
           tl.pos = 'tp',
           tl.srt=30,
           tl.col = "black",
  ) 
  corrplot(data, 
           method="number", 
           type="lower", 
           # order="hclust" ,
           col = 'black', 
           tl.pos = 'n',
           cl.pos = 'n',
           add=TRUE
  ) 
}


# specify metal mixture
mixture_PFAS = c('PFOS_Aug21',	'PFOA_Aug21',	'PFNA_Aug21',	"PFDA_Aug21", "PFHpS_Aug21", "PFHxS_Aug21", "PFHpA_Aug21")
PFAS_name = c('PFOS',	'PFOA',	'PFNA',	"PFDA", "PFHpS", "PFHxS", "PFHpA")

# mixture_PFAS

PFAS_corr = cor(data.frame(BioMe_proteome_PFAS_wide[,..mixture_PFAS]), method = "spearman")
colnames(PFAS_corr)<- PFAS_name[1:7]
rownames(PFAS_corr)<- PFAS_name[1:7]

jpeg("~/Projects/BioMe/proteome/output/EDA/corr_pfas.jpeg",
     units="in", width=8, height=6, res=500)


corr_fun(PFAS_corr) 

dev.off()


###################################### Proteome
## missing pattern
vis_miss((BioMe_proteome_PFAS_wide %>% select(starts_with("OID"))), warn_large_data = FALSE)

## whole distribution
protein_distribution<-  ggplot(BioMe_proteome_PFAS_long) +
                        geom_density(aes(x = NPX, fill=OlinkID), alpha=0.3) +
                        theme(legend.position = "none")



jpeg("~/Projects/BioMe/proteome/output/EDA/distribution.jpeg",
     units="in", width=20, height=16, res=500)

protein_distribution

dev.off()

### protein number for each sample
BioMe_proteome<- BioMe_proteome_PFAS_wide %>% 
                 select(starts_with("OID"))

summary(rowSums(is.na(BioMe_proteome)==FALSE))


number<- data.frame(number = rowSums(is.na(BioMe_proteome)==FALSE))

number_plot1<- ggplot(number, aes(x=number)) + 
                geom_histogram(fill = "#2E5FA1")+
                labs(y ="Frequency",x = "Protein number",
                     title = "Available Protein # for each Sample") +
                theme(plot.title = element_text(size = 16, face="bold"),
                      axis.text = element_text(size = 12, face="bold"),
                      axis.title = element_text(size = 14, face="bold", vjust=0.5),
                      axis.title.x = element_text(size = 14, face="bold", vjust=-0.8),
                      legend.title = element_text(size = 14, face="bold"),
                      legend.text = element_text(size = 12, face="bold")) 

jpeg("~/Projects/BioMe/proteome/output/EDA/number_plot.jpeg",
     units="in", width=14, height=10, res=500)

number_plot1

dev.off()


### sample size for each proteins
count<- BioMe_proteome_PFAS_long %>% 
        group_by(OlinkID) %>% 
        dplyr::summarise(count = n())

summary(count$count)

sample_size_plot1<- ggplot(count, aes(x=count)) + 
                    geom_histogram(fill = "#2E5FA1")+
                    labs(y ="Frequency",x = "Sample Size",
                         title = "Sample Size for Proteins") +
                    theme(plot.title = element_text(size = 16, face="bold"),
                          axis.text = element_text(size = 12, face="bold"),
                          axis.title = element_text(size = 14, face="bold", vjust=0.5),
                          axis.title.x = element_text(size = 14, face="bold", vjust=-0.8),
                          legend.title = element_text(size = 14, face="bold"),
                          legend.text = element_text(size = 12, face="bold")) 

jpeg("~/Projects/BioMe/proteome/output/EDA/sample_size_plot.jpeg",
     units="in", width=14, height=10, res=500)

sample_size_plot1

dev.off()


### median
Median<-  BioMe_proteome_PFAS_long %>% 
                            group_by(OlinkID) %>% 
                            dplyr::summarise(median = median(NPX, na.rm = TRUE))

summary(Median$median)

  
median_plot1<- ggplot(Median, aes(x=reorder(OlinkID, median), y=median)) + 
               geom_bar(stat = "identity", color = "#2E5FA1")+
               labs(y ="Median",x = "Proteins",
                    title = "Distribution of NPX median from all proteins") +
               theme(axis.text.x=element_blank())

jpeg("~/Projects/BioMe/proteome/output/EDA/median_plot.jpeg",
     units="in", width=30, height=10, res=500)

median_plot1

dev.off()


### standard deviation
SD<-  BioMe_proteome_PFAS_long %>% 
      group_by(OlinkID) %>% 
      dplyr::summarise(sd = sd(NPX, na.rm = TRUE))

summary(SD$sd)

sd_plot1<-  ggplot(SD, aes(x=reorder(OlinkID, sd), y=sd)) + 
            geom_bar(stat = "identity", color = "#2E5FA1")+
            labs(y ="SD",x = "Proteins",
                 title = "Distribution of NPX standard deviation from all proteins") +
            theme(axis.text.x=element_blank())


jpeg("~/Projects/BioMe/proteome/output/EDA/sd_plot.jpeg",
     units="in", width=30, height=10, res=500)

sd_plot1

dev.off()


