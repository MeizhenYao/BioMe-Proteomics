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



#------------------------------------------- remove unqualified proteins & assign NA to warning and remove & replace LOD values 
colnames(BioMe_proteome_protname)[12]<- "NPX_raw"
removed_protein<- BioMe_proteome_remove$OlinkID
BioMe_proteome_protname_clean<- BioMe_proteome_protname %>% 
                                filter(!OlinkID %in% removed_protein) %>% 
                                mutate(NPX_LOD=if_else(LOD_ind==0, LOD/sqrt(2), NPX_raw),
                                       NPX=if_else(Warning=="WARN", NA, NPX_LOD)) %>% 
                                filter(is.na(NPX)==FALSE)





#------------------------------------------- transpose to wide format & combine with DID
BioMe_proteome_wide<- BioMe_proteome_protname_clean %>% 
                      select(SampleID, OlinkID, NPX) %>% 
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
                 left_join(PFAS, by = "SampleID")


## Function: Quantile
new_quantile <- function(x, cuts ){
  
  y <- x[x!= 0 & !is.na(x)]
  qi <- unique(quantile(y, probs = seq(0, 1, by = 1/cuts), na.rm = TRUE))
  
  if(length(qi) == 1){ 
    qi = c(-Inf, qi)
  } else{ 
    qi[1] <- -Inf
    qi[length(qi)] <- Inf
  }
  
  x[which(x!= 0 & !is.na(x))] = cut(x[x!= 0 & !is.na(x)], breaks = qi, labels = FALSE, include.lowest = TRUE)
  
  return(x)
  
}


#--------------------------- PFDA
## dichotomous PFDA 
PFAS_epi_AUG21$PFDA_Aug21_bi<- PFAS_epi_AUG21$PFDA_Aug21
PFDA_Aug21_q<- quantile(PFAS_epi_AUG21$PFDA_Aug21, probs = c(1/2), type = 2, na.rm = TRUE)
for(i in 1:nrow(PFAS_epi_AUG21)){
  if (is.na(PFAS_epi_AUG21$PFDA_Aug21[i]) == TRUE){PFAS_epi_AUG21$PFDA_Aug21_bi[i] <- NA}
  else if (PFAS_epi_AUG21$PFDA_Aug21[i] <= PFDA_Aug21_q[1]){PFAS_epi_AUG21$PFDA_Aug21_bi[i] <- "Lower PFDA"}
  else if (PFAS_epi_AUG21$PFDA_Aug21[i] > PFDA_Aug21_q[1]){PFAS_epi_AUG21$PFDA_Aug21_bi[i] <- "Higher PFDA"}
}

PFAS_epi_AUG21$PFDA_Aug21_bi<- factor(PFAS_epi_AUG21$PFDA_Aug21_bi,
                                      levels = c("Lower PFDA", "Higher PFDA"))


## Tertile PFDA - but continuous
PFAS_epi_AUG21$PFDA_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFDA_Aug21, cuts=3)


## continuous PFDA - LOD/SQRT(2) 
PFAS_epi_AUG21$PFDA_Aug21<- ifelse(PFAS_epi_AUG21$PFDA_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFDA_Aug21)



#--------------------------- PFOA
## dichotomous PFOA 
PFAS_epi_AUG21$PFOA_Aug21_bi<- PFAS_epi_AUG21$PFOA_Aug21
PFOA_Aug21_q<- quantile(PFAS_epi_AUG21$PFOA_Aug21, probs = c(1/2), type = 2, na.rm = TRUE)
for(i in 1:nrow(PFAS_epi_AUG21)){
  if (is.na(PFAS_epi_AUG21$PFOA_Aug21[i]) == TRUE){PFAS_epi_AUG21$PFOA_Aug21_bi[i] <- NA}
  else if (PFAS_epi_AUG21$PFOA_Aug21[i] <= PFOA_Aug21_q[1]){PFAS_epi_AUG21$PFOA_Aug21_bi[i] <- "Lower PFOA"}
  else if (PFAS_epi_AUG21$PFOA_Aug21[i] > PFOA_Aug21_q[1]){PFAS_epi_AUG21$PFOA_Aug21_bi[i] <- "Higher PFOA"}
}

PFAS_epi_AUG21$PFOA_Aug21_bi<- factor(PFAS_epi_AUG21$PFOA_Aug21_bi,
                                      levels = c("Lower PFOA", "Higher PFOA"))

## Tertile PFDA - but continuous
PFAS_epi_AUG21$PFOA_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFOA_Aug21, cuts=3)


## continuous PFOA - LOD/SQRT(2) 
PFAS_epi_AUG21$PFOA_Aug21<- ifelse(PFAS_epi_AUG21$PFOA_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFOA_Aug21)



#--------------------------- PFOS
## dichotomous PFOS 
PFAS_epi_AUG21$PFOS_Aug21_bi<- PFAS_epi_AUG21$PFOS_Aug21
PFOS_Aug21_q<- quantile(PFAS_epi_AUG21$PFOS_Aug21, probs = c(1/2), type = 2, na.rm = TRUE)
for(i in 1:nrow(PFAS_epi_AUG21)){
  if (is.na(PFAS_epi_AUG21$PFOS_Aug21[i]) == TRUE){PFAS_epi_AUG21$PFOS_Aug21_bi[i] <- NA}
  else if (PFAS_epi_AUG21$PFOS_Aug21[i] <= PFOS_Aug21_q[1]){PFAS_epi_AUG21$PFOS_Aug21_bi[i] <- "Lower PFOS"}
  else if (PFAS_epi_AUG21$PFOS_Aug21[i] > PFOS_Aug21_q[1]){PFAS_epi_AUG21$PFOS_Aug21_bi[i] <- "Higher PFOS"}
}

PFAS_epi_AUG21$PFOS_Aug21_bi<- factor(PFAS_epi_AUG21$PFOS_Aug21_bi,
                                      levels = c("Lower PFOS", "Higher PFOS"))


## Tertile PFDA - but continuous
PFAS_epi_AUG21$PFOS_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFOS_Aug21, cuts=3)


## continuous PFOS - LOD/SQRT(2) 
PFAS_epi_AUG21$PFOS_Aug21<- ifelse(PFAS_epi_AUG21$PFOS_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFOS_Aug21)



#--------------------------- PFHpA
## dichotomous PFHpA 
PFAS_epi_AUG21$PFHpA_Aug21_bi<- PFAS_epi_AUG21$PFHpA_Aug21
PFHpA_Aug21_q<- quantile(PFAS_epi_AUG21$PFHpA_Aug21, probs = c(1/2), type = 2, na.rm = TRUE)
for(i in 1:nrow(PFAS_epi_AUG21)){
  if (is.na(PFAS_epi_AUG21$PFHpA_Aug21[i]) == TRUE){PFAS_epi_AUG21$PFHpA_Aug21_bi[i] <- NA}
  else if (PFAS_epi_AUG21$PFHpA_Aug21[i] <= PFHpA_Aug21_q[1]){PFAS_epi_AUG21$PFHpA_Aug21_bi[i] <- "Lower PFHpA"}
  else if (PFAS_epi_AUG21$PFHpA_Aug21[i] > PFHpA_Aug21_q[1]){PFAS_epi_AUG21$PFHpA_Aug21_bi[i] <- "Higher PFHpA"}
}

PFAS_epi_AUG21$PFHpA_Aug21_bi<- factor(PFAS_epi_AUG21$PFHpA_Aug21_bi,
                                      levels = c("Lower PFHpA", "Higher PFHpA"))


## Tertile PFHpA - but continuous
PFAS_epi_AUG21$PFHpA_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFHpA_Aug21, cuts=3)


## continuous PFHpA - LOD/SQRT(2) 
PFAS_epi_AUG21$PFHpA_Aug21<- ifelse(PFAS_epi_AUG21$PFHpA_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFHpA_Aug21)




#--------------------------- PFHxS
## dichotomous PFHxS 
PFAS_epi_AUG21$PFHxS_Aug21_bi<- PFAS_epi_AUG21$PFHxS_Aug21
PFHxS_Aug21_q<- quantile(PFAS_epi_AUG21$PFHxS_Aug21, probs = c(1/2), type = 2, na.rm = TRUE)
for(i in 1:nrow(PFAS_epi_AUG21)){
  if (is.na(PFAS_epi_AUG21$PFHxS_Aug21[i]) == TRUE){PFAS_epi_AUG21$PFHxS_Aug21_bi[i] <- NA}
  else if (PFAS_epi_AUG21$PFHxS_Aug21[i] <= PFHxS_Aug21_q[1]){PFAS_epi_AUG21$PFHxS_Aug21_bi[i] <- "Lower PFHxS"}
  else if (PFAS_epi_AUG21$PFHxS_Aug21[i] > PFHxS_Aug21_q[1]){PFAS_epi_AUG21$PFHxS_Aug21_bi[i] <- "Higher PFHxS"}
}

PFAS_epi_AUG21$PFHxS_Aug21_bi<- factor(PFAS_epi_AUG21$PFHxS_Aug21_bi,
                                       levels = c("Lower PFHxS", "Higher PFHxS"))


## Tertile PFHxS - but continuous
PFAS_epi_AUG21$PFHxS_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFHxS_Aug21, cuts=3)



## continuous PFHxS - LOD/SQRT(2) 
PFAS_epi_AUG21$PFHxS_Aug21<- ifelse(PFAS_epi_AUG21$PFHxS_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFHxS_Aug21)


#--------------------------- PFNA
## dichotomous PFNA 
PFAS_epi_AUG21$PFNA_Aug21_bi<- PFAS_epi_AUG21$PFNA_Aug21
PFNA_Aug21_q<- quantile(PFAS_epi_AUG21$PFNA_Aug21, probs = c(1/2), type = 2, na.rm = TRUE)
for(i in 1:nrow(PFAS_epi_AUG21)){
  if (is.na(PFAS_epi_AUG21$PFNA_Aug21[i]) == TRUE){PFAS_epi_AUG21$PFNA_Aug21_bi[i] <- NA}
  else if (PFAS_epi_AUG21$PFNA_Aug21[i] <= PFNA_Aug21_q[1]){PFAS_epi_AUG21$PFNA_Aug21_bi[i] <- "Lower PFNA"}
  else if (PFAS_epi_AUG21$PFNA_Aug21[i] > PFNA_Aug21_q[1]){PFAS_epi_AUG21$PFNA_Aug21_bi[i] <- "Higher PFNA"}
}

PFAS_epi_AUG21$PFNA_Aug21_bi<- factor(PFAS_epi_AUG21$PFNA_Aug21_bi,
                                       levels = c("Lower PFNA", "Higher PFNA"))


## Tertile PFNA - but continuous
PFAS_epi_AUG21$PFNA_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFNA_Aug21, cuts=3)


## continuous PFNA - LOD/SQRT(2) 
PFAS_epi_AUG21$PFNA_Aug21<- ifelse(PFAS_epi_AUG21$PFNA_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFNA_Aug21)


#--------------------------- PFHpS
## dichotomous PFHpS 
PFAS_epi_AUG21$PFHpS_Aug21_bi<- PFAS_epi_AUG21$PFHpS_Aug21
PFHpS_Aug21_q<- quantile(PFAS_epi_AUG21$PFHpS_Aug21, probs = c(1/2), type = 2, na.rm = TRUE)
for(i in 1:nrow(PFAS_epi_AUG21)){
  if (is.na(PFAS_epi_AUG21$PFHpS_Aug21[i]) == TRUE){PFAS_epi_AUG21$PFHpS_Aug21_bi[i] <- NA}
  else if (PFAS_epi_AUG21$PFHpS_Aug21[i] <= PFHpS_Aug21_q[1]){PFAS_epi_AUG21$PFHpS_Aug21_bi[i] <- "Lower PFHpS"}
  else if (PFAS_epi_AUG21$PFHpS_Aug21[i] > PFHpS_Aug21_q[1]){PFAS_epi_AUG21$PFHpS_Aug21_bi[i] <- "Higher PFHpS"}
}

PFAS_epi_AUG21$PFHpS_Aug21_bi<- factor(PFAS_epi_AUG21$PFHpS_Aug21_bi,
                                      levels = c("Lower PFHpS", "Higher PFHpS"))


## Tertile PFHpS - but continuous
PFAS_epi_AUG21$PFHpS_Aug21_q<- new_quantile(PFAS_epi_AUG21$PFHpS_Aug21, cuts=3)



## continuous PFHpS - LOD/SQRT(2) 
PFAS_epi_AUG21$PFHpS_Aug21<- ifelse(PFAS_epi_AUG21$PFHpS_Aug21 < 0.01, 0.01/sqrt(2), PFAS_epi_AUG21$PFHpS_Aug21)



#------------------------------------------- Combine Proteomics & PFAS & EPI data 
BioMe_proteome_PFAS_wide<- BioMe_proteome_wide %>% 
                           inner_join(subset(PFAS_epi_AUG21, select = -c(SampleID)), by="DID")

BioMe_proteome_PFAS_long<- BioMe_proteome_protname_clean %>% 
                           inner_join(subset(PFAS_epi_AUG21, select = -c(SampleID)), by="DID")
                           





write.table(BioMe_proteome_PFAS_wide, "~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_wide.txt", row.names = FALSE)

write.table(BioMe_proteome_PFAS_long, "~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_long.txt", row.names = FALSE)






#------------------------------------------- protein in each panel

BioMe_proteome_PFAS_inflammation<-  BioMe_proteome_PFAS_long %>% 
                                    filter(Panel == "Inflammation")%>% 
                                    select(SampleID, OlinkID, NPX) %>% 
                                    pivot_wider(
                                      names_from = OlinkID,
                                      values_from = NPX
                                    )
  
  
inflammation_OLinkID<- data.frame(protein = colnames(BioMe_proteome_PFAS_inflammation[-1])) %>% 
                       mutate(Panel = "Inflammation")

BioMe_proteome_PFAS_inflammation2<- BioMe_proteome_PFAS_long %>% 
                                    filter(Panel == "Inflammation_II")%>% 
                                    select(SampleID, OlinkID, NPX) %>% 
                                    pivot_wider(
                                      names_from = OlinkID,
                                      values_from = NPX
                                    )

inflammation2_OLinkID<- data.frame(protein = colnames(BioMe_proteome_PFAS_inflammation2[-1])) %>% 
                        mutate(Panel = "Inflammation_II")

BioMe_proteome_PFAS_Cardiometabolic<- BioMe_proteome_PFAS_long %>% 
                                      filter(Panel == "Cardiometabolic")%>% 
                                      select(SampleID, OlinkID, NPX) %>% 
                                      pivot_wider(
                                        names_from = OlinkID,
                                        values_from = NPX
                                      )


cardiometabolic_OLinkID<- data.frame(protein = colnames(BioMe_proteome_PFAS_Cardiometabolic[-1])) %>% 
                          mutate(Panel = "Cardiometabolic")


BioMe_proteome_PFAS_Cardiometabolic2<-  BioMe_proteome_PFAS_long %>% 
                                        filter(Panel == "Cardiometabolic_II")%>% 
                                        select(SampleID, OlinkID, NPX) %>% 
                                        pivot_wider(
                                          names_from = OlinkID,
                                          values_from = NPX
                                        )


cardiometabolic2_OLinkID<- data.frame(protein = colnames(BioMe_proteome_PFAS_Cardiometabolic2[-1])) %>% 
                           mutate(Panel = "Cardiometabolic_II")



BioMe_proteome_PFAS_Neurology<- BioMe_proteome_PFAS_long %>% 
                                filter(Panel == "Neurology")%>% 
                                select(SampleID, OlinkID, NPX) %>% 
                                pivot_wider(
                                  names_from = OlinkID,
                                  values_from = NPX
                                )


neurology_OLinkID<- data.frame(protein = colnames(BioMe_proteome_PFAS_Neurology[-1])) %>% 
                    mutate(Panel = "Neurology")


BioMe_proteome_PFAS_Neurology2<-  BioMe_proteome_PFAS_long %>% 
                                  filter(Panel == "Neurology_II")%>% 
                                  select(SampleID, OlinkID, NPX) %>% 
                                  pivot_wider(
                                    names_from = OlinkID,
                                    values_from = NPX
                                  )


neurology2_OLinkID<- data.frame(protein = colnames(BioMe_proteome_PFAS_Neurology2[-1])) %>% 
                     mutate(Panel = "Neurology_II")

BioMe_proteome_PFAS_Oncology<-  BioMe_proteome_PFAS_long %>% 
                                filter(Panel == "Oncology")%>% 
                                select(SampleID, OlinkID, NPX) %>% 
                                pivot_wider(
                                  names_from = OlinkID,
                                  values_from = NPX
                                )


oncology_OLinkID<-  data.frame(protein = colnames(BioMe_proteome_PFAS_Oncology[-1])) %>% 
                    mutate(Panel = "Oncology")

BioMe_proteome_PFAS_Oncology2<- BioMe_proteome_PFAS_long %>% 
                                filter(Panel == "Oncology_II")%>% 
                                select(SampleID, OlinkID, NPX) %>% 
                                pivot_wider(
                                  names_from = OlinkID,
                                  values_from = NPX
                                )


oncology2_OLinkID<-  data.frame(protein = colnames(BioMe_proteome_PFAS_Oncology2[-1])) %>% 
                     mutate(Panel = "Oncology_II")


protein_in_panel<-  rbind(cardiometabolic_OLinkID,
                          cardiometabolic2_OLinkID,
                          inflammation_OLinkID,
                          inflammation2_OLinkID,
                          neurology_OLinkID,
                          neurology2_OLinkID,
                          oncology_OLinkID,
                          oncology2_OLinkID)



# write.table(protein_in_panel, "~/Projects/BioMe/proteome/input/analysis_sample/protein_in_panel.txt", row.names = FALSE)










