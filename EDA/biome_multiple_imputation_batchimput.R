library(foreach)
library(doParallel)
library(iterators)
library(parallel)
library(tidyverse)
library(ggpubr)
library(data.table)
library(readr)
library(readxl)
library(mice)
library(hdImpute)



BioMe_proteome_PFAS_wide <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_wide.txt")
BioMe_proteome_PFAS_long <- fread("~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_long.txt")



#------------------------------------------- Multiple imputation for covariates
covariates<-  BioMe_proteome_PFAS_wide %>% 
  dplyr::select(PFDA_Aug21, PFHpA_Aug21, PFHxS_Aug21, PFHpS_Aug21, PFNA_Aug21, PFOA_Aug21,
                PFOS_Aug21, self_reported_race, gender, age_at_enrollment, smoking_at_enrollment, bmi_at_enrollment)

covariates$smoking_at_enrollment<- as.factor(covariates$smoking_at_enrollment)                        
sapply(covariates, function(x) sum(is.na(x)))  # to check missingness


# Initialize the Imputation
init = mice(covariates, maxit=0) 
meth = init$method
predM = init$predictorMatrix

set.seed(23761)
imputed_covariates1 = mice(covariates, method=meth, predictorMatrix=predM, m=10, print =  FALSE)

# Use the second completed data set (you can choose any one of the 10-20 completed Datasets)
imputed_covariates <- complete(imputed_covariates1, action=10)


# Check if you actually imputed all the variables you wanted to impute
sapply(imputed_covariates, function(x) sum(is.na(x)))


# combine imputed variables in dataset
BioMe_proteome_PFAS_wide$smoking_at_enrollment_imputed<- imputed_covariates$smoking_at_enrollment
BioMe_proteome_PFAS_wide$bmi_at_enrollment_imputed<- imputed_covariates$bmi_at_enrollment




#------------------------------------------- Multiple imputation for proteins
proteins<-  BioMe_proteome_PFAS_wide %>% 
            dplyr::select(starts_with("OID"))



# stage 1: calculate correlation matrix and store as matrix/array
all_cor <- feature_cor(proteins)

# stage 2: flatten and rank the features by absolute correlation and store as df/tbl
flat_mat <- flatten_mat(all_cor) # can set return_mat = TRUE to print the flattened and ranked feature list

# stage 3: impute, join, and return
imputed_proteins <- impute_batches(data = proteins,
                                   features = flat_mat, 
                                   batch = 10,
                                   n_trees = 100,
                                   seed = 31524)


# stage 4: check if imputation is done for all variables

table(sapply(imputed_proteins, function(x) sum(is.na(x)))) 


# stage 5: evaluate imputation using mean absolute difference (MAD) scores

evaluation<-  mad(original = proteins,
              imputed = imputed_proteins,
              round = 3)


evaluation_plot<- evaluation %>%
                  ggplot(aes(x = mad)) +
                  geom_histogram(fill = "dark green") +
                  labs(x = "MAD Scores (%)", y = "Count of Variables", title = "Distribution of MAD Scores (batch=100)") +
                  theme_minimal() +
                  theme(legend.position = "none")

evaluation_plot


# jpeg("~/Projects/BioMe/proteome/output/EDA/eva3.jpeg",
#      units="in", width=8, height=6, res=500)
# 
# 
# evaluation_plot
# 
# dev.off()


# stage 6: combine with other data

BioMe_proteome_PFAS_wide_imputed<- cbind((BioMe_proteome_PFAS_wide %>% 
                                            dplyr::select(-starts_with("OID"))), imputed_proteins)



#------------------------------------------- transpose to long format

BioMe_proteome_PFAS_long_imputed<- BioMe_proteome_PFAS_wide_imputed %>% 
                                    pivot_longer(
                                      cols = starts_with("OID"),
                                      names_to = "OlinkID",
                                      values_to = "NPX"
                                    ) %>% 
                                   left_join(BioMe_proteome_PFAS_long[,c("OlinkID", "Panel", "UniProt", "Protein_name", "Gene_name")], by = "OlinkID")


write.table(BioMe_proteome_PFAS_wide_imputed, "~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_wide_imputed.txt", row.names = FALSE)
write.table(BioMe_proteome_PFAS_long_imputed, "~/Projects/BioMe/proteome/input/analysis_sample/BioMe_proteome_PFAS_long_imputed.txt", row.names = FALSE)

