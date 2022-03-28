library(PKPDsim)
library(PKPDposterior)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pkvoriconazolefriberg)
library(posterior)

# set Stan path
Sys.setenv(STAN_PATH = "/home/dominic@insight-rx.com/Torsten")

# Compile or reload model
# this creates a binary in the installed package folder
# at ~/R/x86_64-pc-linux-gnu-library/4.1/PKPDposterior/models/
mod <- load_model(
  "pk_voriconazole_friberg", 
  force = T
)

## define init values (use population values): 
prior <- get_init(
  "pkvoriconazolefriberg", 
  map = list("V1" = "V"), 
  drop = c("T50", "TDM_INIT", "F1", "TLAG", "BCF")
)

## Define regimen, covariates, and TDM data
regimen <- PKPDsim::new_regimen(
  amt = 1000, 
  n = 4, 
  times = c(0, 12, 24, 36), 
  type = 'infusion',
  t_inf = 2
)
covariates <- list(
  WT = PKPDsim::new_covariate(value = 70, unit = "kg")
)
tdm_data <- data.frame(
  t = c(3, 23), 
  dv = c(10, 4.2)
)

## Create combined dataset for Torsten/Stan to read:
data <- prepare_data(
  regimen, 
  covariates, 
  tdm_data,
  dose_cmt = 2
)

## Sample from posterior
post <- get_mcmc_posterior(
  mod = mod,
  data = data,
  init = prior,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.95
)

## Show info about the posterior draws
post

## Plot parameter distributions
plot_params(post)

## Simulate from posterior and prior:
covariates$AGE <- new_covariate(25)
covariates$CYP2C19unknown <- new_covariate(0)
covariates$CYP2C19a1a2 <- new_covariate(0)
covariates$CYP2C19a2a2 <- new_covariate(0)
covariates$CYP2C19a2a3 <- new_covariate(0)
covariates$CYP2C19a1a3 <- new_covariate(0)
covariates$CYP2C19a3a3 <- new_covariate(0)
pred_post <- sim_from_draws(
  post, 
  model = pkvoriconazolefriberg::model(),
  map = list("V" = "V1"), 
  parameters = list(T50 = 2.41, TDM_INIT = 0, F1 = 0, TLAG = 0, BCF = 0),
  regimen = regimen,
  covariates = covariates,
  n = 200,
  summarize = TRUE
)
pred_prior <- sim_from_draws(
  post, 
  model = pkvoriconazolefriberg::model(), 
  map = list("V" = "V1"), 
  parameters = list(T50 = 2.41, TDM_INIT = 0, F1 = 0, TLAG = 0, BCF = 0),
  regimen = regimen,
  covariates = covariates,
  prior = TRUE,
  n = 200,
  summarize = TRUE
)

## Plot confidence interval posterior and prior predictions
plot_predictions(pred_prior, obs = tdm_data)
plot_predictions(pred_post, obs = tdm_data)

## Plot posterior AUC distribution
pred_post_full <- sim_from_draws( # don't summarize
  post, 
  map = list("V" = "V1"), 
  parameters = list(T50 = 2.41, TDM_INIT = 0, F1 = 0, TLAG = 0, BCF = 0),
  model = pkvoriconazolefriberg::model(), 
  regimen = regimen,
  covariates = covariates,
  n = 200
)
pred_post_full %>%
  filter(t %in% c(36, 48)) %>%
  filter(comp == 3) %>%
  group_by(id) %>%
  tidyr::pivot_wider(names_from = t, values_from = y) %>%
  mutate(auc24 = 2 * (`48` - `36`)) %>%
  ggplot() + 
  aes(x = auc24) +
  geom_histogram()
