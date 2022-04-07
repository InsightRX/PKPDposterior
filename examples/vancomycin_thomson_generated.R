## Vancomycin PK model (Thomson et al)
## Implemented using analytic equation
## Stan model generated using new_stan_model()

library(PKPDposterior)

## Vancomycin model (Thomson et al)
parameters <- list(CL = 2.99, Q = 2.28, V2 = 0.732, V1 = 0.675)
iiv <- list(CL = 0.27, Q = 0.49, V1 = 0.15, V2 = 1.3)
ruv <- list(add = 1.6, prop = 0.15)

model <- new_stan_model(
  parameters = parameters,
  parameter_definitions = list(
    "CL" = "CL * (1.0 + 0.0154 * ((CRCL[j] * 16.6667) - 66.0))",
    "Q"  = "Q",
    "V1" = "V1 * WT[j]",
    "V2" = "V2 * WT[j]",
    "KA" = "0" 
  ),
  covariate_definitions = list(
    "CRCL" = "real", # CrCl in L/hr
    "WT" = "real"    # WT in kg
  ),
  solver = "pmx_solve_twocpt",
  scale = "(V1 * mean(WT))",
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
mapping <- list("V1" = "V")

## Define regimen, covariates, and TDM data
regimen <- PKPDsim::new_regimen(
  amt = 1500, 
  n = 4, 
  times = c(0, 12, 24, 36), 
  type = 'infusion',
  t_inf = 2
)
covariates <- list(
  WT = PKPDsim::new_covariate(value = 70, unit = "kg"),
  CRCL = PKPDsim::new_covariate(value = 5, unit = "l/hr")
)
tdm_data <- data.frame(
  t = c(2.5, 11.5), 
  dv = c(40, 14)
)

## Create combined dataset for Torsten/Stan to read:
data <- PKPDsim_to_stan_data(
  regimen,
  covariates, 
  tdm_data,
  dose_cmt = 2,
  parameters = parameters,
  iiv = iiv,
  ruv = list(
    prop = 0.15,
    add = 1.6
  ),
  ltbs = FALSE
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
