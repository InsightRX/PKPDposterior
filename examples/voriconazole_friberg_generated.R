## Voriconazole PK model (Friberg et al)
## Implemented using ODE system
## Stan model generated using new_stan_model()

library(PKPDposterior)

parameters <- list(KA = 1.19, VMAX1 = 114, 
                   KM = 1.15, CL = 6.16, Q = 15.5, 
                   V2 = 103, V1 = 79)
iiv <- list(CL = 0.5, Q = 0.42, V1 = 0.136, V2 = 0.77, KA = 0.9, KM = 1.0, VMAX1 = 0.5)
ruv <- list(prop = 0.3, add = 0.01)

ode <- c(
  "real k10 = CL / V1;",
  "real k12 = Q / V1;",
  "real k21 = Q / V2;",
  "real VMAXINH = exp(1.5) / (1 + exp(1.5));",
  "real Vmax = VMAX1 * (1 - VMAXINH * (t-1) / ((t-1) + (2.41 - 1)));",
  "dAdt[1] = -KA*A[1];",
  "dAdt[2] =  KA*A[1] - (k10 + k12)*A[2] + k21*A[3] - (Vmax * A[2]/V1)/(KM + A[2]/V1);",
  "dAdt[3] = k12*A[2] - k21*A[3];"
)
parameter_definitions <- list(
  "CL" = "CL * pow(mean(WT)/70, 0.75)",
  "Q" = "Q * 1.637 * pow(mean(WT)/70, 0.75)",
  "V1" = "V1 * mean(WT)/70",
  "V2" = "V2 * mean(WT)/70",
  "KA" = "KA",
  "KM" = "KM",
  "VMAX1" = "VMAX1 * pow(mean(WT)/70, 0.75)"
)

model <- new_stan_model(
  parameters = parameters,
  parameter_definitions = parameter_definitions,
  ode = ode,
  covariate_definitions = list(
    "WT" = "real"    # WT in kg
  ),
  solver = "pmx_solve_rk45",
  obs_cmt = 2,
  scale = "V1 * pow(mean(WT)/70, 0.75)",
  verbose = T
)
model_file <- write_stan_model(model)

# Compile or reload model
mod <- load_model(
  model_file, 
  force = T,
  verbose = T
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
data <- prepare_data(
  regimen, 
  covariates, 
  tdm_data,
  dose_cmt = 2,
  parameters = parameters,
  iiv = iiv,
  ruv = ruv
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
