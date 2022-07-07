## Busulfan PK model in children (Shukla et al)
## Implemented using analytic equation
## Stan model generated using new_stan_model()

library(PKPDsim)
library(PKPDposterior)
library(pkbusulfanucsf)

## Busulfan model & parameters
parameters <- list(
  CL = 3.96,
  V = 10.8,
  MAT_MAG = 0.451, 
  K_MAT = 1.37, 
  TH_REGI = -0.2, 
  TH_DAY = -0.135
)
iiv <- list(CL = 0.24, V = 0.17)
ruv <- list(prop = 0.106, add = 22.2)

model <- new_stan_model(
  parameters = parameters,
  fixed = c("MAT_MAG", "K_MAT", "TH_REGI", "TH_DAY"),
  variable_definitions = list(
    "BMI" = "WT[j]/(HT[j]*HT[j]/10000)",
    "FFM" = "(SEX[j] * (0.88 + ((1-0.88)/(1+(AGE[j]/13.4)^-12.7))) * ((9270 * WT[j])/(6680 + (216 * BMI)))) + ((1-SEX[j])*(1.11 + ((1-1.11)/(1+(AGE[j]/7.1)^-1.1))) * ((9270 * WT[j])/(8780 + (244 * BMI))))"
  ),
  parameter_definitions = list(
    "CL" = "CL * (FFM/12)^0.75 * (MAT_MAG + (1-MAT_MAG)*(1-exp(-AGE[j]*K_MAT))) * (1 + (-TH_REGI * REGI[j])) * (1 -TH_DAY*(time[j]>24))",
    "V"  = "V * FFM/12",
    "KA" = "0" 
  ),
  covariate_definitions = list(
    "SEX" = "int",
    "WT" = "real",
    "HT" = "real",
    "AGE" = "real",
    "REGI" = "int"
  ),
  solver = "pmx_solve_onecpt",
  scale = "((V * FFM/12)/1000)",
  verbose = T
)

model_file <- write_stan_model(model)

# Compile or reload model
mod <- load_model(
  model_file, 
  force = T,
  verbose = T
)

## mapping of parameter names between PKPDsim and Stan/Torsten
mapping <- list()

## Define regimen, covariates, and TDM data
regimen <- new_regimen(
  amt = 150, 
  n = 4, 
  times = c(0, 24, 48, 72), 
  type = 'infusion',
  t_inf = 3
)
covariates <- list(
  WT = new_covariate(value = 25, unit = "kg"),
  HT = new_covariate(value = 100, unit = "cm"),
  SEX = new_covariate(value = 0), # male
  AGE = new_covariate(value = 6, unit = "l/hr"),
  REGI = new_covariate(0)
)
tdm_data <- data.frame(
  t = c(3.25, 6), 
  dv = c(3600, 1200)
)

## allow to specify extra calcualtion lines
## `iiv` object requires same list elements as `parameters` object --> add ". Fix parameters"
data <- new_stan_data(
  regimen,
  covariates,
  tdm_data,
  dose_cmt = 2,
  parameters = parameters,
  fixed = c("MAT_MAG", "K_MAT", "TH_REGI", "TH_DAY"),
  iiv = iiv,
  ruv = ruv,
  ltbs = FALSE
)

## Sample from posterior
post <- get_mcmc_posterior(
  mod = mod,
  data = data,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.8
)

## Show info about the posterior draws
post

## Plot parameter distributions
plot_params(post)

## Validate model vs PKPDsim implementation
validate_stan_model(
  stan_model = mod,
  pkpdsim_model = pkbusulfanucsf::model(),
  max_abs_delta = 0.1,
  data = data,
  mapping = mapping
)

## simulate
pred_post <- sim_from_draws(
  post,
  model = pkbusulfanucsf::model(),
  map = mapping,
  regimen = regimen,
  covariates = covariates,
  n = 200,
  summarize = TRUE
)
pred_prior <- sim_from_draws(
  post,
  model = pkbusulfanucsf::model(),
  map = mapping,
  regimen = regimen,
  covariates = covariates,
  n = 200,
  prior = TRUE,
  summarize = TRUE
)

plot_predictions(pred_post, obs = tdm_data) +
  ggplot2::scale_y_log10()

plot_predictions(pred_prior, obs = tdm_data) +
  ggplot2::scale_y_log10()
