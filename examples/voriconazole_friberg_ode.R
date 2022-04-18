## Voriconazole PK model (Friberg et al)
## Implemented using ODE system

library(PKPDsim)
library(PKPDposterior)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pkvoriconazolefriberg)
library(posterior)
cmdstanr::set_cmdstan_path(
  path = file.path(Sys.getenv("STAN_PATH"), "cmdstan")
)

mapping <- list("V1" = "V")

# Compile or reload model
# this creates a binary in the installed package folder
# at ~/R/x86_64-pc-linux-gnu-library/4.1/PKPDposterior/models/
mod <- load_model(
  "pk_voriconazole_friberg",
  force = T,
  verbose = T
)

## define init values (use population values): 
prior <- prior_from_PKPDsim_model(
  "pkvoriconazolefriberg", 
  map = mapping, 
  drop = c("T50", "TDM_INIT", "F1", "TLAG", "BCF")
)

## Define regimen, covariates, and TDM data
regimen <- PKPDsim::new_regimen(
  amt = 500, 
  n = 7, 
  times = c(0, 12, 24, 36, 48, 60, 72), 
  type = 'infusion',
  t_inf = 2
)
covariates <- list(
  WT = PKPDsim::new_covariate(value = 70, unit = "kg")
)
tdm_data <- data.frame(
  t = c(59, 70.5), 
  dv = c(1.6, 1.5)
)

## Create combined dataset for Torsten/Stan to read:
data <- new_stan_data(
  regimen, 
  covariates, 
  tdm_data,
  dose_cmt = 2,
  parameters = prior,
  iiv = list(
    CL = 0.5,
    Q = 0.42,
    V1 = 0.136,
    V2 = 0.77,
    KA = 0.9,
    KM = 1.0,
    VMAX1 = 0.5
  ),
  ruv = list(
    prop = 0.15,
    add = 1.6
  )
)

## Sample from posterior
post <- get_mcmc_posterior(
  mod = mod,
  data = data,
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
  map = mapping, 
  parameters = list(T50 = 2.41, TDM_INIT = 0, F1 = 0, TLAG = 0, BCF = 0),
  regimen = regimen,
  covariates = covariates,
  n = 200,
  summarize = TRUE
)
pred_prior <- sim_from_draws(
  post,
  model = pkvoriconazolefriberg::model(), 
  map = mapping, 
  parameters = list(T50 = 2.41, TDM_INIT = 0, F1 = 0, TLAG = 0, BCF = 0),
  regimen = regimen,
  covariates = covariates,
  prior = TRUE,
  n = 200,
  summarize = TRUE
)

## Plot confidence interval posterior and prior predictions
plot_predictions(pred_prior, obs = tdm_data) +
  geom_hline(yintercept = 1.15)
plot_predictions(pred_post, obs = tdm_data) +
  geom_hline(yintercept = 1.15)
