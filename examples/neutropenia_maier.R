## Neutropenia PK-PD example, taken to Maier et al. CPT 2020 and CPT-PSP 2022.

library(PKPDsim)
library(PKPDposterior)
library(pkpdneutropeniatemplate1)
library(dplyr)
library(ggplot2)
cmdstanr::set_cmdstan_path(
  path = file.path(Sys.getenv("STAN_PATH"), "cmdstan")
)
parameters <- list(
  V1 = 10.8, 
  V3 = 301,
  KMel = 0.667,  
  VMAXel = 35.9, 
  KMtr = 1.44, 
  VMAXtr = 175,
  k21 = 1.12,
  Q = 16.8, 
  ALPHA = 13.1, 
  MTT = 145, 
  CIRC0 = 6.48, 
  GAMMA = 0.257
)
iiv <- list(
  SLOPE = 1.0, MTT = 0.2, 
  CIRC0 = 0.5, GAMMA = 0.2
)
ruv <- list(
  pk = list(add = 0.2),
  pd = list(add = 0.3)
)
parameter_definitions <- list(
  "VMAXel" = "VMAXel",
  "KMel" = "KMel",
  "VMAXtr" = "VMAXtr",
  "KMtr" = "KMtr",
  "Q" = "Q",
  "V1" = "V1",
  "V3" = "V3",
  "k21" = "k21",
  "MTT" = "MTT",
  "CIRC0" = "CIRC0",
  "GAMMA" = "GAMMA",
  "ALPHA" = "ALPHA"
)
ode <- "
  real ktr = 4 / MTT;
  real x_ratio = 0.787;
  real ktr_prol = x_ratio * ktr;
  real ktr_stem = (1 - x_ratio) * ktr;
  real C1 = A[1] / V1;
  real k13 = Q/V1;
  real k31 = Q/V3;
  real conc;
  real Eff;
  real transit1;
  real transit2;
  real transit3;
  real circ;
  real stem;
  real prol;

  // PK paclitaxel
  // real VMAXel_i = VMAXel * (mean(BSA)/1.8)^1.14 * 1.07^SEX * (mean(AGE)/56)^-0.447 * (mean(BILI)/7)^-0.0942;
  real VMAXel_i = VMAXel;
  dAdt[1] = -(VMAXel_i/V1)/(KMel+C1) -(VMAXtr/V1)/(KMtr+C1) + k21*A[2] - k13*A[1] + k31*A[3];
  dAdt[2] =                           (VMAXtr/V1)/(KMtr+C1) - k21*A[2];
  dAdt[3] =  k13 * A[1] - k31 * A[3];

  // PD neutrophils
  Eff = ALPHA * C1; // slope model, not Emax
  stem = A[4] + CIRC0;
  prol = A[5] + CIRC0;
  transit1 = A[6] + CIRC0;
  transit2 = A[7] + CIRC0;
  transit3 = A[8] + CIRC0;
  circ = fmax(machine_precision(), A[9] + CIRC0);
  
  dAdt[4] = ktr_stem * stem * ((1 - Eff) * ((CIRC0 / circ)^GAMMA) - 1);
  dAdt[5] = ktr_prol * prol * ((1 - Eff) * ((CIRC0 / circ)^GAMMA) - 1) + ktr * stem;
  dAdt[6] = ktr * (prol - transit1);
  dAdt[7] = ktr * (transit1 - transit2);
  dAdt[8] = ktr * (transit2 - transit3);
  dAdt[9] = ktr * (transit3 - circ);
"

model <- new_stan_model(
  parameters = parameters,
  parameter_definitions = parameter_definitions,
  ode = ode,
  covariate_definitions = list(
    "BSA" = "real", 
    "AGE" = "real", 
    "SEX" = "real", 
    "BILI" = "real"
  ),
  solver = 'pmx_solve_rk45',
  obs_types = c("pk", "pd"),
  custom_ipred = list(
    "pk" = "A[1, ] ./ V1;",
    "pd" = "A[9, ] + theta[7];"
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
  iter_warmup = 250,
  iter_sampling = 250,
  adapt_delta = 0.95,
  verbose = TRUE
)

## Show posterior info
post

## Plot parameter distributions
plot_params(post)
