library(car)
library(readr)
library(lattice)
library(nlme)
library(ggplot2)
library(GGally)
library(nnet)
library(foreign)
library(biotools)
library(glmmML)
library(MASS)
library(lme4)
library(multcomp)
library(dplyr)
library(knitr)
library(xtable)
library(kableExtra)
library(DT)
library(glmnet)
library(corrplot)
library(ggpubr)
library(lmerTest)
library("merTools")
library(reshape2)
library(ggplot2)
library(GGally)
library(mgcv)
library(gplots)
library(tidyr)
library(bkmr)
library(factoextra) 
library(spatstat)
library(Hmisc)
library(blme)
library(grpreg)
library(robustHD)
library(gWQS)
library(gridExtra)
library(ggcorrplot)
library(BWQS)
library(qwraps2)
library(MatchIt)
library(data.table)
library(mice)
library(ggrepel)
library(spatstat)

cores=detectCores()
cl <- makeCluster(10) 
registerDoParallel(cl)

start.time <- Sys.time()


##------------------------------------------- import data
BioMe_proteome_PFAS_wide <- fread("/sc/arion/work/yaom03/biome_proteome/dataset/BioMe_proteome_PFAS_wide.txt")



##------------------------------------------- prepare data
## date of blood draw
BioMe_proteome_PFAS_wide$date_enrl <- rep(NA_real_, nrow(BioMe_proteome_PFAS_wide))
BioMe_proteome_PFAS_wide$month_yr_enrl<- as.character(BioMe_proteome_PFAS_wide$month_yr_enrl)


for(i in 1:nrow(BioMe_proteome_PFAS_wide)){
  
  x <- anytime::anydate(paste((strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][2]), " 1,", 2000 + as.numeric(strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][1])))
  mydates <- as.Date(c("2011-01-01"))
  BioMe_proteome_PFAS_wide$date_enrl[i] <- as.numeric((x - mydates[1])/365 )
  BioMe_proteome_PFAS_wide$year_enrl[i] <- round(2000 + as.numeric(strsplit(BioMe_proteome_PFAS_wide$month_yr_enrl[i],"-")[[1]][1]), 0)
  
}

BioMe_proteome_PFAS_wide$c_date_enrl <- ifelse(BioMe_proteome_PFAS_wide$date_enrl > 0, 1,0)


bwqs_data<- BioMe_proteome_PFAS_wide %>% 
            select(starts_with("OID"), ends_with("_q"), self_reported_race, gender, age_at_enrollment, smoking_at_enrollment, c_date_enrl, ipw)


bwqs_data_dummy<- as.data.frame(dummify(bwqs_data))



##------------------------------------------- bwqs model fitting
bwqs_data_proteins<- bwqs_data_dummy %>% 
                     select(starts_with("OID"))

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


bwqs_pfas_met_model <- data.frame(mean = NA_real_, se_mean = NA_real_, sd = NA_real_,
                                lower = NA_real_, upper = NA_real_, n_eff = NA_real_,
                                Rhat = NA_real_)



start.time <- Sys.time()

for(i in 1:2612){
  ## specify parameter
  y_name  <- colnames(bwqs_data_proteins)[i]
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
  
  
  bwqs_pfas_met_model[i] <-  as.numeric(sum_fit_lasso[3,])
  
  
}

write.csv(bwqs_pfas_met_model, "/sc/arion/work/yaom03/biome_proteome/pfas_proteome/proteome_vs_pfas_bwqs.csv",
          row.names = F)

end.time <- Sys.time()
(time.taken <- end.time - start.time)

stopCluster(cl)








