library(PKPDsim)
library(PKPDposterior)
library(pkpdneutropeniatemplate1)
library(dplyr)
library(ggplot2)
cmdstanr::set_cmdstan_path(
  path = file.path(Sys.getenv("STAN_PATH"), "cmdstan")
)
mapping <- list( # mapping between parameter names Stan vs PKPDsim
  "V1" = "V"
)

parameters <- list(CL = 5, V1 = 50,
                   SLOPE = 0.1, MTT = 100, 
                   CIRC0 = 5, GAMMA = 0.2)
iiv <- list(
  CL = 0.2, V1 = 0.5, 
  SLOPE = 1.0, MTT = 0.2, 
  CIRC0 = 0.5, GAMMA = 0.2
)
ruv <- list(
  pk = list(add = 0.2),
  pd = list(add = 0.3)
)
mapping <- list( # mapping between parameter names Stan vs PKPDsim
  "V1" = "V", 
  "gamma" = "GAMMA",
  "mtt" = "MTT",
  "alpha" = "SLOPE",
  "circ0" = "CIRC0"
)
parameter_definitions <- list(
  "CL" = "CL",
  "Q" = 0,
  "V1" = "V1",
  "V2" = 1,
  "KA" = 0,
  "mtt" = "MTT",
  "circ0" = "CIRC0",
  "gamma" = "GAMMA",
  "alpha" = "SLOPE"
)
ode <- "
  real ka  = 0;
  real k10 = CL / V1;
  real k12 = 0; // Q / V1;
  real k21 = 0; // Q / V2;
  real ktr = 4 / mtt;
  
  real conc;
  real EDrug;
  real transit1;
  real transit2;
  real transit3;
  real circ;
  real prol;
  
  dAdt[1] = -ka * A[1];
  dAdt[2] =  ka * A[1] - (k10 + k12) * A[2] + k21 * A[3];
  dAdt[3] = k12 * A[2] - k21 * A[3];
  conc = A[2] / V1;
  
  EDrug = alpha * conc; // slope model, not Emax
  prol = A[4] + circ0;
  transit1 = A[5] + circ0;
  transit2 = A[6] + circ0;
  transit3 = A[7] + circ0;
  circ = fmax(machine_precision(), A[8] + circ0); // Device for implementing a modeled 
  
  // initial condition
  dAdt[4] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
  dAdt[5] = ktr * (prol - transit1);
  dAdt[6] = ktr * (transit1 - transit2);
  dAdt[7] = ktr * (transit2 - transit3);
  dAdt[8] = ktr * (transit3 - circ);
"

model_file <- new_stan_model(
  parameters = parameters,
  parameter_definitions = parameter_definitions,
  ode = ode,
  covariate_definitions = NULL,
  n_cmt = 8,
  obs_types = c("pk", "pd"),
  custom_ipred = list(
    "pk" = "A[2, ] ./ V1;",
    "pd" = "A[8, ] + theta[7];"
  ),
  verbose = T
)

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
  init = parameters,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.95,
  verbose = TRUE,
  skip_processing = TRUE
)

extract_map_estimates(post)

## Plot parameter distributions
plot_params(post)
