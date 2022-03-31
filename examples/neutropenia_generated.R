## Neutropenia PK-PD example
## Implemented using ODE system for both PK and PD
## Stan model generated using new_stan_model()

library(PKPDsim)
library(PKPDposterior)
library(pkpdneutropeniatemplate1)
library(dplyr)
library(ggplot2)
cmdstanr::set_cmdstan_path(
  path = file.path(Sys.getenv("STAN_PATH"), "cmdstan")
)
parameters <- list(CL = 5, V = 50,
                   SLOPE = 0.1, MTT = 100, 
                   CIRC0 = 5, GAMMA = 0.2)
iiv <- list(
  CL = 0.2, V = 0.5, 
  SLOPE = 1.0, MTT = 0.2, 
  CIRC0 = 0.5, GAMMA = 0.2
)
ruv <- list(
  pk = list(add = 0.2),
  pd = list(add = 0.3)
)
parameter_definitions <- list(
  "CL" = "CL",
  "Q" = 0,
  "V" = "V",
  "V2" = 1,
  "KA" = 0,
  "MTT" = "MTT",
  "CIRC0" = "CIRC0",
  "GAMMA" = "GAMMA",
  "ALPHA" = "SLOPE"
)
ode <- "
  real ka  = 0;
  real k10 = CL / V;
  real k12 = 0;
  real k21 = 0;
  real ktr = 4 / MTT;
  
  real conc;
  real EDrug;
  real transit1;
  real transit2;
  real transit3;
  real circ;
  real prol;
  
  dAdt[1] = -KA * A[1];
  dAdt[2] =  KA * A[1] - (k10 + k12) * A[2] + k21 * A[3];
  dAdt[3] = k12 * A[2] - k21 * A[3];
  conc = A[2] / V;
  
  EDrug = ALPHA * conc; // slope model, not Emax
  prol = A[4] + CIRC0;
  transit1 = A[5] + CIRC0;
  transit2 = A[6] + CIRC0;
  transit3 = A[7] + CIRC0;
  circ = fmax(machine_precision(), A[8] + CIRC0); // Device for implementing a modeled 
  
  // initial condition
  dAdt[4] = ktr * prol * ((1 - EDrug) * ((CIRC0 / circ)^GAMMA) - 1);
  dAdt[5] = ktr * (prol - transit1);
  dAdt[6] = ktr * (transit1 - transit2);
  dAdt[7] = ktr * (transit2 - transit3);
  dAdt[8] = ktr * (transit3 - circ);
"

model <- new_stan_model(
  parameters = parameters,
  parameter_definitions = parameter_definitions,
  ode = ode,
  covariate_definitions = NULL,
  solver = 'pmx_solve_rk45',
  obs_types = c("pk", "pd"),
  custom_ipred = list(
    "pk" = "A[2, ] ./ V;",
    "pd" = "A[8, ] + theta[7];"
  ),
  verbose = T
)
model_file <- write_stan_model(model)
mod <- load_model(
  model_file = model_file,
  force = T,
  verbose = T
)

## Define regimen, covariates, and TDM data
regimen <- new_regimen(
  amt = 1500, 
  n = 3, 
  times = c(0, 24, 48), 
  type = 'infusion',
  cmt = 2,
  t_inf = 2
)
covariates <- NULL
tdm_data <- data.frame(
  t = c(2.5, 11.5), 
  dv = c(40, 14),
  type = "pk"
)
pd_data <- data.frame(
  t = c(3, 6, 9, 12) * 24, 
  dv = c(5, 1.5, .8, 2),
  type = "pd"
)
comb_data <- bind_rows(
  tdm_data, 
  pd_data
)

## Create combined dataset for Torsten/Stan to read:
data <- prepare_data(
  regimen,
  covariates = covariates, 
  data = comb_data,
  parameters = parameters,
  iiv = iiv,
  ltbs = list(pk = TRUE, pd = TRUE),
  ruv = ruv,
  dose_cmt = 2
)

## Sample from posterior
post <- get_mcmc_posterior(
  mod = mod,
  data = data,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.95,
  verbose = TRUE,
  skip_processing = TRUE
)

extract_map_estimates(post)

## Plot parameter distributions
plot_params(post)
