library(PKPDsim)
library(PKPDposterior)
library(pkpdneutropeniatemplate1)
library(dplyr)
library(ggplot2)
cmdstanr::set_cmdstan_path(
  path = file.path(Sys.getenv("STAN_PATH"), "cmdstan")
)

mapping <- list( # mapping between parameter names Stan vs PKPDsim
  "V1" = "V", 
  "gamma" = "GAMMA",
  "mtt" = "MTT",
  "alpha" = "SLOPE",
  "circ0" = "CIRC0"
)

mod <- load_model(
  "pkpd_neutropenia_v1",
  force = T,
  verbose = T
)

## define init values (use population values): 
prior <- get_init(
  "pkpd_neutropenia_template1", 
  map = mapping,
  drop = c("KA", "Q", "V2")
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
  parameters = prior,
  iiv = list(
    CL = 0.2,
    V1 = 0.5,
    mtt = 0.2,
    circ0 = 0.2,
    alpha = 1.0,
    gamma = 0.5
  ),
  ltbs = list(pk = TRUE, pd = TRUE),
  ruv = list(
    pk = list(add = 0.2),
    pd = list(add = 0.3)
  ),
  dose_cmt = 2
)

## Sample from posterior
post <- get_mcmc_posterior(
  mod = mod,
  data = data,
  init = prior,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.95,
  verbose = TRUE,
  skip_processing = TRUE
)
extract_map_estimates(post)

## Plot parameter distributions
plot_params(post)

## Simulate from posterior
pred <- sim_from_draws(
  post, 
  model = pkpdneutropeniatemplate1::model(),
  map = mapping,
  parameters = list(KA = 0, Q = 0, V2 = 1),
  regimen = regimen,
  n = 200,
  t_obs = seq(0, 15*24, 6),
  summarize = T
)
plot_predictions(pred, obs = pd_data) +
  xlim(c(0.1, 300)) +
  scale_y_log10(limits = c(0.15, 8)) +
  vpc::theme_empty()

## simulate from prior
pred_prior <- sim_from_draws(
  post, 
  model = pkpdneutropeniatemplate1::model(),
  map = mapping,
  regimen = regimen,
  parameters = list(KA = 0, Q = 0, V2 = 1),
  n = 200,
  prior = TRUE,
  t_obs = seq(0, 15*24, 6),
  summarize = T
)
plot_predictions(pred_prior, obs = pd_data) +
  xlim(c(0.1, 300)) +
  scale_y_log10(limits = c(0.01, 8)) +
  vpc::theme_empty()

