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


cores=detectCores()
cl <- makeCluster(10) 
registerDoParallel(cl)

start.time <- Sys.time()

BioMe_proteome_PFAS_wide <- fread("/sc/arion/work/yaom03/biome_proteome/dataset/BioMe_proteome_PFAS_wide.txt")



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
imputed_covariates1 = mice(covariates, method=meth, predictorMatrix=predM, m=10)

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

sapply(proteins, function(x) sum(is.na(x)))  # to check missingness


# Initialize the Imputation
init = mice(proteins, maxit=0) 
meth = init$method
predM = init$predictorMatrix

set.seed(23761)
imputed_proteins1 = mice(proteins, method=meth, predictorMatrix=predM, m=10)

# Use the second completed data set (you can choose any one of the 10-20 completed Datasets)
imputed_proteins <- complete(imputed_proteins1, action=10)


# Check if you actually imputed all the variables you wanted to impute
sapply(imputed_proteins, function(x) sum(is.na(x)))

BioMe_proteome_PFAS_wide_imputed<- cbind((BioMe_proteome_PFAS_wide %>% 
                                            dplyr::select(-starts_with("OID"))), imputed_proteins)


write.table(BioMe_proteome_PFAS_wide_imputed, "/sc/arion/work/yaom03/biome_proteome/dataset/BioMe_proteome_PFAS_wide_imputed.txt", row.names = FALSE)


end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)