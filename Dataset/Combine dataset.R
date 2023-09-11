library(tidyverse)
library(ggpubr)
library(data.table)
library(readr)
library(readxl)



#------------------------------------------- import data
BioMe_proteome <- fread("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/Documents/Projects/BioMe/proteome/input/Q-00771_Loos_NPX_2022-05-03/Q-00771_Loos_NPX_2022-05-03.txt")
assay_list<- read_excel("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/Documents/Projects/BioMe/proteome/input/assay-list-olink-explore-3072.xlsx", sheet = "data")
BioMe_proteome_remove<- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_warn_remove.txt")
BioMe_customer <- read_excel("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/Documents/Projects/BioMe/proteome/input/Q-00771_Customer_sample_manifest_-_consolidated.xlsx")

merged<- read_csv("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/Documents/Projects/BioMe/proteome/input/merged_pfas_epi_liver_data.csv")
PFAS<- read_excel("C:/Users/yaom03/OneDrive - The Mount Sinai Hospital/Documents/Projects/BioMe/proteome/input/New_PFAS_data_Aug_21.xlsx")

PFAS<- PFAS %>%
       filter(Sample_Class == "Study_Sample")


colnames(merged)[101]<- "SampleID"
colnames(PFAS)[4]<- "SampleID"
colnames(BioMe_customer)[2]<- "SampleID"
colnames(BioMe_customer)[7]<- "DID"

# PFAS$SampleID<- gsub("AB", "", PFAS$SampleID)

#------------------------------------------- remove controls
BioMe_proteome_nocontrol<- BioMe_proteome %>% 
                           filter(!grepl("CONTROL", SampleID))


#------------------------------------------- combine data with annotation
str(BioMe_proteome_nocontrol)
str(assay_list)

colnames(assay_list)[1]<- "UniProt"
colnames(assay_list)[2]<- "Protein_name"
colnames(assay_list)[3]<- "Gene_name"
colnames(assay_list)[4]<- "Panel"


BioMe_proteome_protname<- BioMe_proteome_nocontrol %>% 
                          left_join(assay_list, by = c("Panel", "UniProt"))


# colSums(is.na(BioMe_proteome_protname))


#------------------------------------------- assign WARN and LOD indicator column
BioMe_proteome_protname<- BioMe_proteome_protname %>% 
                          mutate(Warning = if_else(QC_Warning == "PASS" & Assay_Warning == "PASS", "PASS", "WARN"),
                                 LOD_ind = if_else(NPX >= LOD, 1, 0)) %>% 
                          inner_join(BioMe_customer[, c("SampleID", "DID")], by="SampleID")




#------------------------------------------- calculate WARN and LOD percent for each proteins (cutoff 10% and 25% for warning and LOD)
BioMe_proteome_QC_data<-  BioMe_proteome_protname %>% 
                          dplyr::select(Panel, Protein_name, Gene_name, OlinkID, Warning) %>% 
                          group_by(Panel, Protein_name, Gene_name, OlinkID, Warning) %>% 
                          summarise(Warning_num=n())%>%
                          mutate(Warning_sign = case_when(Warning=="PASS" & Warning_num==352 ~ 0, 
                                                          Warning=="PASS" & Warning_num!=352 ~ NA,
                                                          Warning=="WARN" ~ Warning_num)) %>% 
                          filter(is.na(Warning_sign) == FALSE) %>% 
                          mutate(Warning_percent = Warning_sign/352) %>% 
                          dplyr::select(Panel, Protein_name, Gene_name, OlinkID, Warning_percent) %>% 
                          as.data.frame() 




BioMe_proteome_LOD_data<- BioMe_proteome_protname %>% 
                          dplyr::select(Panel, Protein_name, Gene_name, OlinkID, LOD_ind) %>% 
                          group_by(Panel, Protein_name, Gene_name, OlinkID, LOD_ind) %>% 
                          summarise(LOD_num=n())%>%
                          mutate(LOD_sign = case_when(LOD_ind==0 & LOD_num==352 ~ 0, 
                                                      LOD_ind==0 & LOD_num!=352 ~ NA,
                                                      LOD_ind==1 ~ LOD_num)) %>% 
                          filter(is.na(LOD_sign) == FALSE) %>% 
                          mutate(Percentage_above_LOD = LOD_sign/352) %>% 
                          dplyr::select(Panel, Protein_name, Gene_name, OlinkID, Percentage_above_LOD) %>% 
                          as.data.frame() 

BioMe_proteome_check_data<- BioMe_proteome_QC_data %>% 
                            inner_join(BioMe_proteome_LOD_data, by=c('Panel', 'Protein_name', 'Gene_name', 'OlinkID'))



BioMe_proteome_remove<- BioMe_proteome_check_data %>% 
                        filter(Warning_percent > 0.10 | Percentage_above_LOD < 0.25 | 
                                 OlinkID %in% c("OID20101", "OID20911", "OID21276",
                                                "OID30230", "OID30563", "OID31050",
                                                "OID20153", "OID20631", "OID20997", 
                                                "OID30212", "OID30547", "OID31416",
                                                "OID30225", "OID30558", "OID31046", 
                                                "OID20074", "OID20473", "OID20848"))

# write.table(BioMe_proteome_remove, "~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_warn_remove.txt", row.names = FALSE)



#------------------------------------------- remove unqualified proteins & assign NA to warning and remove & replace LOD values 
colnames(BioMe_proteome_protname)[12]<- "NPX_raw"
removed_protein<- BioMe_proteome_remove$OlinkID
BioMe_proteome_protname_clean<- BioMe_proteome_protname %>% 
                                filter(!OlinkID %in% removed_protein) %>% 
                                mutate(NPX_LOD=if_else(LOD_ind==0, LOD/sqrt(2), NPX_raw),
                                       NPX=if_else(Warning=="WARN", NA, NPX_LOD)) %>% 
                                filter(is.na(NPX)==FALSE)

BioMe_proteome_protname_removed<- BioMe_proteome_protname %>% 
                                filter(OlinkID %in% removed_protein) %>% 
                                mutate(NPX_LOD=if_else(LOD_ind==0, LOD/sqrt(2), NPX_raw),
                                       NPX=if_else(Warning=="WARN", NA, NPX_LOD)) %>% 
                                filter(is.na(NPX)==FALSE)





#------------------------------------------- transpose to wide format & combine with DID
BioMe_proteome_wide<- BioMe_proteome_protname_clean %>% 
                      dplyr::select(SampleID, OlinkID, NPX) %>% 
                      pivot_wider(
                        names_from = OlinkID,
                        values_from = NPX
                      ) %>% 
                      inner_join(BioMe_customer[, c("SampleID", "DID")], by="SampleID")

BioMe_proteome_wide_removed<- BioMe_proteome_protname_removed %>% 
                              dplyr::select(SampleID, OlinkID, NPX) %>% 
                              pivot_wider(
                                names_from = OlinkID,
                                values_from = NPX
                              ) %>% 
                              inner_join(BioMe_customer[, c("SampleID", "DID")], by="SampleID")





length(intersect(BioMe_proteome_wide$SampleID, BioMe_customer$SampleID))
length(intersect(BioMe_customer$DID, merged$DID))




#------------------------------------------- Combine PFAS & EPI data 
merged<- subset(merged, select = -c(PFDeA, PFHxA, PFHxS, PFNA, PFOA, PFOS))
PFAS$SampleID<- as.numeric(gsub("AB", "", PFAS$SampleID))

length(intersect(merged$SampleID, PFAS$SampleID))

PFAS_epi_AUG21<- merged %>% 
                 inner_join(PFAS, by = "SampleID")


##################################################
# function to format variables in quantile
##################################################
new_quantile <- function(x, cuts ){
  
  qi <- unique(quantile(x, probs = seq(0, 1, by = 1/cuts), type = 2, na.rm = TRUE))
  
  if(length(qi) == 1){ 
    qi = c(-Inf, qi)
  } else{ 
    qi[1] <- -Inf
    qi[length(qi)] <- Inf
  }
  
  x[which(!is.na(x))] = cut(x[!is.na(x)], breaks = qi, labels = FALSE, include.lowest = TRUE)
  
  return(x)
  
}


#--------------------------- PFDA
## dichotomous PFDA 
PFAS_epi_AUG21$PFDA_Aug21_bi<- new_quantile(PFAS_epi_AUG21$PFDA_Aug21, cuts=2)


PFAS_epi_AUG21$PFDA_Aug21_bi<- factor(PFAS_epi_AUG21$PFDA_Aug21_bi,
                                      levels = c(1, 2),
                                      labels = c("Lower PFDA", "Higher PFDA"))


## Quantile PFDA - but continuous

PFAS_epi_AUG21$PFDA_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFDA_Aug21, cuts=4)


## continuous PFDA - LOD/SQRT(2) 
PFAS_epi_AUG21$PFDA_Aug21<- ifelse(PFAS_epi_AUG21$PFDA_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFDA_Aug21)



#--------------------------- PFOA
## dichotomous PFOA 
PFAS_epi_AUG21$PFOA_Aug21_bi<- new_quantile(PFAS_epi_AUG21$PFOA_Aug21, cuts=2)


PFAS_epi_AUG21$PFOA_Aug21_bi<- factor(PFAS_epi_AUG21$PFOA_Aug21_bi,
                                      levels = c(1, 2),
                                      labels = c("Lower PFOA", "Higher PFOA"))

## Tertile PFDA - but continuous
PFAS_epi_AUG21$PFOA_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFOA_Aug21, cuts=4)


## continuous PFOA - LOD/SQRT(2) 
PFAS_epi_AUG21$PFOA_Aug21<- ifelse(PFAS_epi_AUG21$PFOA_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFOA_Aug21)



#--------------------------- PFOS
## dichotomous PFOS 
PFAS_epi_AUG21$PFOS_Aug21_bi<- new_quantile(PFAS_epi_AUG21$PFOS_Aug21, cuts=2)


PFAS_epi_AUG21$PFOS_Aug21_bi<- factor(PFAS_epi_AUG21$PFOS_Aug21_bi,
                                      levels = c(1, 2),
                                      labels = c("Lower PFOS", "Higher PFOS"))


## Tertile PFDA - but continuous
PFAS_epi_AUG21$PFOS_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFOS_Aug21, cuts=4)


## continuous PFOS - LOD/SQRT(2) 
PFAS_epi_AUG21$PFOS_Aug21<- ifelse(PFAS_epi_AUG21$PFOS_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFOS_Aug21)



#--------------------------- PFHpA
## dichotomous PFHpA 
PFAS_epi_AUG21$PFHpA_Aug21_bi<- PFAS_epi_AUG21$PFHpA_Aug21

for(i in 1:nrow(PFAS_epi_AUG21)){
  if (is.na(PFAS_epi_AUG21$PFHpA_Aug21[i]) == TRUE){PFAS_epi_AUG21$PFHpA_Aug21_bi[i] <- NA}
  else if (PFAS_epi_AUG21$PFHpA_Aug21[i] < 0.01){PFAS_epi_AUG21$PFHpA_Aug21_bi[i] <- "Lower PFHpA"}
  else {PFAS_epi_AUG21$PFHpA_Aug21_bi[i] <- "Higher PFHpA"}
}

PFAS_epi_AUG21$PFHpA_Aug21_bi<- factor(PFAS_epi_AUG21$PFHpA_Aug21_bi,
                                      levels = c("Lower PFHpA", "Higher PFHpA"))


## Tertile PFHpA - but continuous
PFAS_epi_AUG21$PFHpA_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFHpA_Aug21, cuts=4)


## continuous PFHpA - LOD/SQRT(2) 
PFAS_epi_AUG21$PFHpA_Aug21<- ifelse(PFAS_epi_AUG21$PFHpA_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFHpA_Aug21)




#--------------------------- PFHxS
## dichotomous PFHxS 
PFAS_epi_AUG21$PFHxS_Aug21_bi<- new_quantile(PFAS_epi_AUG21$PFHxS_Aug21, cuts=2)


PFAS_epi_AUG21$PFHxS_Aug21_bi<- factor(PFAS_epi_AUG21$PFHxS_Aug21_bi,
                                      levels = c(1, 2),
                                      labels = c("Lower PFHxS", "Higher PFHxS"))


## Tertile PFHxS - but continuous
PFAS_epi_AUG21$PFHxS_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFHxS_Aug21, cuts=4)



## continuous PFHxS - LOD/SQRT(2) 
PFAS_epi_AUG21$PFHxS_Aug21<- ifelse(PFAS_epi_AUG21$PFHxS_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFHxS_Aug21)


#--------------------------- PFNA
## dichotomous PFNA 
PFAS_epi_AUG21$PFNA_Aug21_bi<- new_quantile(PFAS_epi_AUG21$PFNA_Aug21, cuts=2)


PFAS_epi_AUG21$PFNA_Aug21_bi<- factor(PFAS_epi_AUG21$PFNA_Aug21_bi,
                                       levels = c(1, 2),
                                       labels = c("Lower PFNA", "Higher PFNA"))


## Tertile PFNA - but continuous
PFAS_epi_AUG21$PFNA_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFNA_Aug21, cuts=4)


## continuous PFNA - LOD/SQRT(2) 
PFAS_epi_AUG21$PFNA_Aug21<- ifelse(PFAS_epi_AUG21$PFNA_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFNA_Aug21)


#--------------------------- PFHpS
## dichotomous PFHpS 
PFAS_epi_AUG21$PFHpS_Aug21_bi<- new_quantile(PFAS_epi_AUG21$PFHpS_Aug21, cuts=2)


PFAS_epi_AUG21$PFHpS_Aug21_bi<- factor(PFAS_epi_AUG21$PFHpS_Aug21_bi,
                                      levels = c(1, 2),
                                      labels = c("Lower PFHpS", "Higher PFHpS"))



## Tertile PFHpS - but continuous
PFAS_epi_AUG21$PFHpS_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFHpS_Aug21, cuts=4)



## continuous PFHpS - LOD/SQRT(2) 
PFAS_epi_AUG21$PFHpS_Aug21<- ifelse(PFAS_epi_AUG21$PFHpS_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFHpS_Aug21)



#------------------------------------------- Combine Proteomics & PFAS & EPI data 
BioMe_proteome_PFAS_wide<- BioMe_proteome_wide %>% 
                           inner_join(subset(PFAS_epi_AUG21, select = -c(SampleID)), by="DID")

BioMe_proteome_PFAS_long<- BioMe_proteome_protname_clean %>% 
                           inner_join(subset(PFAS_epi_AUG21, select = -c(SampleID)), by="DID")
                           




# 
# write.table(BioMe_proteome_PFAS_wide, "~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_wide.txt", row.names = FALSE)
# 
# write.table(BioMe_proteome_PFAS_long, "~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_long.txt", row.names = FALSE)
# 





#------------------------------------------- protein in each panel
protein_in_panel<-  BioMe_proteome_PFAS_long %>% 
                    group_by(Panel, OlinkID) %>% 
                    dplyr::summarise(count = n())


# write.table(protein_in_panel, "~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt", row.names = FALSE)










