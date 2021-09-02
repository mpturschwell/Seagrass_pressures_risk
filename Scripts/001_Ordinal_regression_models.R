# ---------------------------------------------------------------------------------------------------------------------------
# Anthropogenic pressures and life-history predict trajectories of seagrass meadow extent at a global scale
# June 2021
# Ordinal seagrass trajectories model combinations and plots
# Bayesian ordinal regression models with unequal variance - ONLY 2000-2010 
# Author: MP Turschwell
# https://github.com/mpturschwell
# ---------------------------------------------------------------------------------------------------------------------------

library(tidyverse)
library(dplyr)
library(brms)
library(patchwork)
library(ggplot2)
library(sf)
library(ggthemes)
library(bayestestR)
library(RColorBrewer)
library(tidybayes)

mytime <- format(Sys.time(), "%Y_%m_%d")

# import data
dat2000 <- read_csv("Data/Data_seagrass_trends_pressures_2000-2010.csv", col_types = cols())

# Fit null model first 
m1_null_no_RE <- brm(category ~ 1, 
               prior = c(prior(normal(0, 4), class = Intercept)),  # weakly informative 
               iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
               data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

# Fit null model with bioregion random effect 
m1_null <- brm(category ~ 1 + (1|bioregion), 
               prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
               prior(cauchy(0, 2), class = sd)),         # weakly informative
               iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
               data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))


# --------------------------------------------------------------------------------------------------
# Fit distributional model that also models variance by lh category 
# --------------------------------------------------------------------------------------------------

# 5km buffer
dist_mod_5k <- brm(
  formula = bf(category ~ mean_Dd_5k +mean_Dh_5k + mean_Np_5k +  mean_Oa_5k +  
                 mean_Ocp_5k + mean_Ship_5k + mean_Slr_5k +  mean_Sst_5k +   
                 Turbidity.Mean_5k + Turbidity.Coefficient.of.Variance_5k +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

# 5km buffer with interaction between LH and predictors 
dist_mod_5k_int <- brm(
  formula = bf(category ~ (mean_Dd_5k +mean_Dh_5k + mean_Np_5k +  mean_Oa_5k +  
                             mean_Ocp_5k + mean_Ship_5k + mean_Slr_5k +  mean_Sst_5k +   
                             Turbidity.Mean_5k + Turbidity.Coefficient.of.Variance_5k) * lh +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

# 10km buffer
dist_mod_10k <- brm(
  formula = bf(category ~ mean_Dd_10k +mean_Dh_10k + mean_Np_10k +  mean_Oa_10k +  
                 mean_Ocp_10k + mean_Ship_10k + mean_Slr_10k +  mean_Sst_10k +   
                 Turbidity.Mean_10k + Turbidity.Coefficient.of.Variance_10k +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

# 10km buffer with interaction between LH and predictors 
dist_mod_10k_int <- brm(
  formula = bf(category ~ (mean_Dd_10k +mean_Dh_10k + mean_Np_10k +  mean_Oa_10k +  
                             mean_Ocp_10k + mean_Ship_10k + mean_Slr_10k +  mean_Sst_10k +   
                             Turbidity.Mean_10k + Turbidity.Coefficient.of.Variance_10k) * lh +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

# 20km buffer
dist_mod_20k <- brm(
  formula = bf(category ~ mean_Dd_20k +mean_Dh_20k + mean_Np_20k +  mean_Oa_20k +  
                 mean_Ocp_20k + mean_Ship_20k + mean_Slr_20k +  mean_Sst_20k +   
                 Turbidity.Mean_20k + Turbidity.Coefficient.of.Variance_20k +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

#20km buffer with interaction between LH and predictors 
dist_mod_20k_int <- brm(
  formula = bf(category ~ (mean_Dd_20k +mean_Dh_20k + mean_Np_20k +  mean_Oa_20k +  
                             mean_Ocp_20k + mean_Ship_20k + mean_Slr_20k +  mean_Sst_20k +   
                             Turbidity.Mean_20k + Turbidity.Coefficient.of.Variance_20k) * lh +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

# 30km buffer
dist_mod_30k <- brm(
  formula = bf(category ~ mean_Dd_30k +mean_Dh_30k + mean_Np_30k +  mean_Oa_30k +  
                 mean_Ocp_30k + mean_Ship_30k + mean_Slr_30k +  mean_Sst_30k +   
                 Turbidity.Mean_30k + Turbidity.Coefficient.of.Variance_30k +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

#30km buffer with interaction between LH and predictors 
dist_mod_30k_int <- brm(
  formula = bf(category ~ (mean_Dd_30k +mean_Dh_30k + mean_Np_30k +  mean_Oa_30k +  
                             mean_Ocp_30k + mean_Ship_30k + mean_Slr_30k +  mean_Sst_30k +   
                             Turbidity.Mean_30k + Turbidity.Coefficient.of.Variance_30k) * lh +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

# 50km buffer
dist_mod_50k <- brm(
  formula = bf(category ~ mean_Dd_50k +mean_Dh_50k + mean_Np_50k +  mean_Oa_50k +  
                             mean_Ocp_50k + mean_Ship_50k + mean_Sst_50k+ mean_Slr_50k +     
                             Turbidity.Mean_50k + Turbidity.Coefficient.of.Variance_50k +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

#50km buffer with interaction between LH and predictors 
dist_mod_50k_int <- brm(
  formula = bf(category ~ (mean_Dd_50k +mean_Dh_50k + mean_Np_50k +  mean_Oa_50k +  
                             mean_Ocp_50k + mean_Ship_50k + mean_Slr_50k +  mean_Sst_50k +   
                             Turbidity.Mean_50k + Turbidity.Coefficient.of.Variance_50k) * lh +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))


# 100km buffer
dist_mod_100k <- brm(
  formula = bf(category ~ mean_Dd_100k +mean_Dh_100k + mean_Np_100k +  mean_Oa_100k +  
                 mean_Ocp_100k + mean_Ship_100k + mean_Slr_100k +  mean_Sst_100k +   
                 Turbidity.Mean_100k + Turbidity.Coefficient.of.Variance_100k +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

#100km buffer with interaction between LH and predictors 
dist_mod_100k_int <- brm(
  formula = bf(category ~ (mean_Dd_100k +mean_Dh_100k + mean_Np_100k +  mean_Oa_100k +  
                             mean_Ocp_100k + mean_Ship_100k + mean_Slr_100k +  mean_Sst_100k +   
                             Turbidity.Mean_100k + Turbidity.Coefficient.of.Variance_100k) * lh +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

# 200km buffer
dist_mod_200k <- brm(
  formula = bf(category ~ mean_Dd_200k +mean_Dh_200k + mean_Np_200k +  mean_Oa_200k +  
                 mean_Ocp_200k + mean_Ship_200k + mean_Slr_200k +  mean_Sst_200k +   
                 Turbidity.Mean_200k + Turbidity.Coefficient.of.Variance_200k +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

#200km buffer with interaction between LH and predictors 
dist_mod_200k_int <- brm(
  formula = bf(category ~ (mean_Dd_200k +mean_Dh_200k + mean_Np_200k +  mean_Oa_200k +  
                             mean_Ocp_200k + mean_Ship_200k + mean_Slr_200k +  mean_Sst_200k +   
                             Turbidity.Mean_200k + Turbidity.Coefficient.of.Variance_200k) * lh +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))


# based on best above model, remove collinear predictors bases on Grech et al. 2012 vulnerability weightings
# removing Sea Level Rise 
dist_mod_100k_noSLR <- brm(
  formula = bf(category ~ mean_Dd_100k + mean_Dh_100k + mean_Np_100k +  mean_Oa_100k + mean_Ocp_100k+  
                 mean_Ship_100k +  mean_Sst_100k +    
                 Turbidity.Mean_100k + Turbidity.Coefficient.of.Variance_100k +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

# removing Sea Level Rise but including interaction between LH and predictors
dist_mod_100k_noSLR_int <- brm(
  formula = bf(category ~ (mean_Dd_100k + mean_Dh_100k + mean_Np_100k +  mean_Oa_100k + mean_Ocp_100k+  
                             mean_Ship_100k +  mean_Sst_100k +    
                             Turbidity.Mean_100k + Turbidity.Coefficient.of.Variance_100k) *lh +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

# removing Sea Level Rise & Population Density 
dist_mod_100k_noSLR_DH <- brm(
  formula = bf(category ~ mean_Dd_100k + mean_Np_100k +  mean_Oa_100k + mean_Ocp_100k+  
                 mean_Ship_100k +  mean_Sst_100k +    
                 Turbidity.Mean_100k + Turbidity.Coefficient.of.Variance_100k +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))

# removing Sea Level Rise & Population Density including interaction between LH and predictors 
dist_mod_100k_noSLR_DH_int <- brm(
  formula = bf(category ~ (mean_Dd_100k + mean_Np_100k +  mean_Oa_100k + mean_Ocp_100k+  
                             mean_Ship_100k +  mean_Sst_100k +    
                             Turbidity.Mean_100k + Turbidity.Coefficient.of.Variance_100k)*lh +
                 (1 | bioregion)) +
    lf(disc ~ 0+lh, cmc = FALSE), 
  prior = c(prior(normal(0, 4), class = Intercept),  # weakly informative 
            prior(normal(0, 2), class = b),          # weakly informative
            prior(exponential(1), class = sd)),         # weakly informative
  iter = 4000, warmup = 1000, cores = 2, chains = 4, thin = 1,
  data = dat2000, family = cumulative("probit"), control = list(max_treedepth = 15, adapt_delta = 0.99))


# add criterion for model comparison 
m1_null          <- add_criterion(m1_null,  c("loo", "waic"))
m1_null_no_RE    <- add_criterion(m1_null_no_RE,      c("loo", "waic"))
dist_mod_5k_int  <- add_criterion(dist_mod_5k_int,  c("loo", "waic"))
dist_mod_5k      <- add_criterion(dist_mod_5k,      c("loo", "waic"))
dist_mod_10k_int <- add_criterion(dist_mod_10k_int, c("loo", "waic"))
dist_mod_10k     <- add_criterion(dist_mod_10k,     c("loo", "waic"))
dist_mod_20k_int <- add_criterion(dist_mod_20k_int, c("loo", "waic"))
dist_mod_20k     <- add_criterion(dist_mod_20k,     c("loo", "waic"))
dist_mod_30k_int <- add_criterion(dist_mod_30k_int, c("loo", "waic"))
dist_mod_30k     <- add_criterion(dist_mod_30k,     c("loo", "waic"))
dist_mod_50k_int <- add_criterion(dist_mod_50k_int, c("loo", "waic"))
dist_mod_50k     <- add_criterion(dist_mod_50k,     c("loo", "waic"))
dist_mod_100k_int<- add_criterion(dist_mod_100k_int, c("loo", "waic"))
dist_mod_100k    <- add_criterion(dist_mod_100k,     c("loo", "waic"))
dist_mod_200k_int<- add_criterion(dist_mod_200k_int, c("loo", "waic"))
dist_mod_200k    <- add_criterion(dist_mod_200k,     c("loo", "waic"))
dist_mod_100k_noSLR    <- add_criterion(dist_mod_100k_noSLR,     c("loo", "waic"))
dist_mod_100k_noSLR_DH    <- add_criterion(dist_mod_100k_noSLR_DH,     c("loo", "waic"))
dist_mod_100k_noSLR_int    <- add_criterion(dist_mod_100k_noSLR_int,     c("loo", "waic"))
dist_mod_100k_noSLR_DH_int    <- add_criterion(dist_mod_100k_noSLR_DH_int,     c("loo", "waic"))

# Save all models 
save(dist_mod_5k_int ,
     dist_mod_5k     ,
     dist_mod_10k_int,
     dist_mod_10k    ,
     dist_mod_20k_int,
     dist_mod_20k    ,
     dist_mod_30k_int,
     dist_mod_30k    ,
     dist_mod_50k_int,
     dist_mod_50k    ,
     dist_mod_100k_int,
     dist_mod_100k    ,
     dist_mod_100k_noSLR,
     dist_mod_100k_noSLR_int,
     dist_mod_100k_noSLR_DH,
     dist_mod_100k_noSLR_DH_int,
     dist_mod_200k_int,
     dist_mod_200k    ,
     m1_null,
     m1_null_no_RE,
     dat, dat2000,
     file = "Data/All_model_combinations.rda")

##################################################################################################################################
load(file = "Data/All_model_combinations.rda")

# Compare model information criterion
print(loo_compare(
          dist_mod_5k_int ,
          dist_mod_10k_int,
          dist_mod_20k_int,
          dist_mod_30k_int,
          dist_mod_50k_int,
          dist_mod_100k_int,
          dist_mod_100k_noSLR_int,
          dist_mod_100k_noSLR_DH_int,
          dist_mod_200k_int,
          dist_mod_5k     ,
          dist_mod_10k    ,
          dist_mod_20k    ,
          dist_mod_30k    ,
          dist_mod_50k,
          dist_mod_100k    ,
          dist_mod_100k_noSLR,
          dist_mod_100k_noSLR_DH,
          dist_mod_200k),
      criterion = "loo")

# only save best model
save(dist_mod_100k_noSLR_DH,file = "Data/brms_model.rda")
#
# Load and look at summary of best model
load(file = "Data/brms_model.rda")
mod_summary <- summary(dist_mod_100k_noSLR_DH)

# Check model residuals
resid <- dat2000 %>%
  dplyr::select(-geom) %>% 
  add_predicted_draws(dist_mod_100k_noSLR_DH) %>%
  mutate(.prediction = ordered(levels(category)[.prediction], levels = levels(category))) %>%
  summarise(
    p_lower = mean(.prediction < category),
    p_upper = mean(.prediction <= category),
    p_residual = runif(1, p_lower, p_upper),
    z_residual = qnorm(p_residual)
  ) %>%
  ggplot(aes(x = .row, y = z_residual)) +
  geom_point()

resid 

residuals <- resid$data

# residuals by LH category 
ggplot(residuals, aes(x = z_residual))+
  geom_density()+
  facet_wrap(~lh)

# residuals by bioregion 
ggplot(residuals, aes(x = z_residual, y = bioregion))+
  geom_boxplot()

#----------------------------------------------------------------------------------------#
# Calculate Standard Deviation of each LH category - remembering colonising is fixed at 1
#----------------------------------------------------------------------------------------#
# disc is inverse of SD: s = 1/disc 
# standard deviations s, we transformed log(disc) to s using: s = 1/ exp(log(disc))
Opportunistic_LH <- 1/exp(mod_summary$fixed[c(13,0)])
Mixed_LH <- 1/exp(mod_summary$fixed[c(14,0)])
Persistent <- 1/exp(mod_summary$fixed[c(15,0)])

