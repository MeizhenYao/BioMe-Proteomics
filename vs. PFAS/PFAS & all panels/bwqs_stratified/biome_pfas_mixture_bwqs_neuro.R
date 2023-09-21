library(car)
library(readr)
library(lattice)
library(nlme)
library(ggplot2)
library(GGally)
library(foreign)
library(MASS)
library(lme4)
library(multcomp)
library(dplyr)
library(knitr)
library(xtable)
library(kableExtra)
library(glmnet)
library(corrplot)
library(ggpubr)
library(lmerTest)
library("merTools")
library(reshape2)
library(gplots)
library(tidyr)
library(blme)
library(grpreg)
library(gridExtra)
library(ggcorrplot)
library(BWQS)
library(data.table)
library(mice)
library(spatstat)
library(foreach)
library(doParallel)
library(iterators)
library(parallel)
library(anytime)

cores=detectCores()
cl <- makeCluster(10) 
registerDoParallel(cl)

start.time <- Sys.time()


##------------------------------------------- import data
BioMe_proteome_PFAS_wide <- fread("/sc/arion/work/yaom03/biome_proteome/pfas_proteome/bwqs/BioMe_proteome_PFAS_wide_imputed.txt")
protein_in_panel <- fread("/sc/arion/work/yaom03/biome_proteome/pfas_proteome/bwqs/protein_in_panel.txt")
BioMe_proteome_PFAS_long <- fread("/sc/arion/work/yaom03/biome_proteome/pfas_proteome/bwqs/BioMe_proteome_PFAS_long.txt")



##------------------------------------------- prepare data
## date of blood draw
BioMe_proteome_PFAS_wide$date_enrl <- rep(NA_real_, nrow(BioMe_proteome_PFAS_wide))
BioMe_proteome_PFAS_wide$month_yr_enrl<- as.character(BioMe_proteome_PFAS_wide$month_yr_enrl)


for(i in 1:nrow(BioMe_proteome_PFAS_wide)){
  
  x <- anytime::anydate(paste((strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][2]), " 1,", 2000 + as.numeric(strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][1])))
  mydates <- as.Date(c("2011-01-01"))
  BioMe_proteome_PFAS_wide$date_enrl[i] <- as.numeric((x - mydates[1])/365 )
  
}

BioMe_proteome_PFAS_wide$c_date_enrl <- ifelse(BioMe_proteome_PFAS_wide$date_enrl > 0, 1,0)


## stratified by status
BioMe_proteome_PFAS_wide_case<- BioMe_proteome_PFAS_wide %>% 
  filter(td2_case_all == 1)

BioMe_proteome_PFAS_wide_control<- BioMe_proteome_PFAS_wide %>% 
  filter(td2_case_all == 0)


bwqs_data<- BioMe_proteome_PFAS_wide_control %>% 
  dplyr::select(starts_with("OID"), ends_with("_q"), self_reported_race, gender, age_at_enrollment, smoking_at_enrollment, c_date_enrl, ipw)


bwqs_data_dummy<- as.data.frame(dummify(bwqs_data))



##------------------------------------------- bwqs model fitting

protein<- (protein_in_panel %>% filter(Panel == "Neurology"))$OlinkID

bwqs_data_proteins<- bwqs_data_dummy %>% 
  dplyr::select(starts_with("OID"))

name_data_wqs <- c("PFDA_Aug21_q","PFHxS_Aug21_q","PFHpS_Aug21_q","PFNA_Aug21_q","PFOA_Aug21_q",
                   "PFOS_Aug21_q")

model_bwqs_gaussian_lasso <- "data {

int<lower=0> N;              // number of individual
int<lower=0> C1;             // number of element in the mix
int<lower=0> K;              // number of covariates
matrix[N,C1] XC1;            // matrix of first mix
matrix[N,K] KV;	             // matrix of covariates
vector[C1] DalpC1;           // vector of the Dirichlet coefficients for first mix
vector[N] sw;                // IPW weights
real y[N];                   // outcome gaussian variable
}

parameters {

real <lower=0> sigma;
real mu;                              // intercepts
real beta;                            // coeffs by group
vector[K] delta;                      // covariates coefficients
real<lower=0> lambda_squared;         // penalization factor
simplex[C1] WC1;                      // weights of first mix

}
transformed parameters {

vector[N] Xb;
Xb = mu + (XC1*WC1)*beta  + KV*delta;
}
model {

mu ~ normal(0, 10);
sigma ~ inv_gamma(0.01,0.01);
lambda_squared ~ gamma(2,0.5);
beta ~ normal(0,lambda_squared);
for(j in 1:K) delta[j] ~ normal(0,K);
WC1 ~ dirichlet(DalpC1);
for(n in 1:N){
  target +=  normal_lpdf(y[n]| Xb[n], sigma) * sw[n];
}
}

"
m_lasso_data_challenge <- rstan::stan_model(model_code =  model_bwqs_gaussian_lasso)

data = bwqs_data_dummy

bwqs_pfas_met_model <- data.frame(mean = NA_real_, se_mean = NA_real_, sd = NA_real_,
                                  lower = NA_real_, upper = NA_real_, n_eff = NA_real_,
                                  Rhat = NA_real_)

bwqs_pfas_weight<- data.frame(w1 = NA_real_, w2 = NA_real_, w3 = NA_real_,
                              w4 = NA_real_, w5 = NA_real_, w6 = NA_real_)


start.time <- Sys.time()

for(i in 1:length(protein)){
  ## specify parameter
  y_name  <- protein[i]
  formula = as.formula( ~ self_reported_race.African.American
                        + self_reported_race.European.American + age_at_enrollment
                        + smoking_at_enrollment.No + gender.Female
                        + c_date_enrl)
  
  KV_name <- all.vars(formula)
  mix_name_1 <- name_data_wqs
  
  X1 = bwqs_data_dummy[,name_data_wqs]
  
  data_reg <- list(
    
    N   = nrow(bwqs_data_dummy),
    C1  = length(mix_name_1),
    XC1 = cbind(X1),
    DalpC1 = rep(1, length(mix_name_1)),
    KV = data[,KV_name],
    K   = length(KV_name),
    sw = as.vector(data[,"ipw"]),
    y = as.vector(data[,y_name])
  )
  
  
  ## model fit
  fit_lasso <- rstan::sampling(m_lasso_data_challenge,
                               data = data_reg,
                               chains = 1,
                               iter = 2e3,
                               thin = 1,
                               refresh = 0, verbose = T,
                               control=list(max_treedepth = 20,
                                            adapt_delta = 0.999999999999999))
  
  
  ## model result summary
  sum_fit_lasso <- (summary(fit_lasso,
                            probs = c(0.025, 0.975))$summary)
  
  
  bwqs_pfas_met_model <-  rbind(bwqs_pfas_met_model,as.numeric(sum_fit_lasso["beta",])) 
  
  bwqs_pfas_weight <-  rbind(bwqs_pfas_weight,as.numeric(sum_fit_lasso[c("WC1[1]", "WC1[2]", "WC1[3]", "WC1[4]", "WC1[5]", "WC1[6]"),"mean"])) 
  
}


bwqs_pfas_met_model <- bwqs_pfas_met_model[-1,]
bwqs_pfas_weight <- bwqs_pfas_weight[-1,]


bwqs_pfas_met_model$OlinkID <- protein
bwqs_pfas_weight$OlinkID <- protein


write.table(bwqs_pfas_met_model, "/sc/arion/work/yaom03/biome_proteome/pfas_proteome/bwqs_control/proteome_vs_pfas_bwqs_neuro_control.txt", row.names = FALSE)
write.table(bwqs_pfas_weight, "/sc/arion/work/yaom03/biome_proteome/pfas_proteome/bwqs_control/bwqs_pfas_weight_neuro_control.txt", row.names = FALSE)


end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)












