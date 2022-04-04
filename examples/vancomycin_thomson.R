## Vancomycin PK model (Thomson et al)
## Implemented using analytic equation

library(PKPDsim)
library(PKPDposterior)
library(ggplot2)
library(pkvancothomson)

mapping <- list("V1" = "V")

# Compile or reload model
mod <- load_model(
  "pk_vanco_thomson", 
  force = T,
  verbose = T
)
 
## define init values (use population values): 
prior <- prior_from_PKPDsim_model(
  "pkvancothomson", 
  map = mapping, 
  drop = c("TH_CRCL", "TDM_INIT")
)

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
  parameters = prior,
  iiv = list(CL = 0.27, Q = 0.49, V1 = 0.15, V2 = 1.3),
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
post_vi <- get_mcmc_posterior(
  mod = mod,
  data = data,
  method = "vi",
  skip_processing = TRUE
)

## Show info about the posterior draws
post

## Plot parameter distributions
plot_params(post)

## Simulate from posterior and prior:
covariates$CL_HEMO <- new_covariate(0)
pred_post <- sim_from_draws(
  post, 
  model = pkvancothomson::model(),
  map = mapping,
  parameters = list(TH_CRCL = 0.0154, TDM_INIT = 0),
  regimen = regimen,
  covariates = covariates,
  n = 200,
  summarize = TRUE
)
pred_post_vi <- sim_from_draws(
  post_vi, 
  model = pkvancothomson::model(),
  map = mapping,
  parameters = list(TH_CRCL = 0.0154, TDM_INIT = 0),
  regimen = regimen,
  covariates = covariates,
  n = 200,
  summarize = TRUE
)
pred_prior <- sim_from_draws(
  post, 
  model = pkvancothomson::model(), 
  map = mapping,
  parameters = list(TH_CRCL = 0.0154, TDM_INIT = 0),
  regimen = regimen,
  covariates = covariates,
  prior = TRUE,
  n = 200,
  summarize = TRUE
)

## Plot confidence interval posterior and prior predictions
plot_predictions(pred_post, obs = tdm_data)
plot_predictions(pred_post_vi, obs = tdm_data)
plot_predictions(pred_prior, obs = tdm_data)

## Plot posterior AUC distribution
pred_post_full <- sim_from_draws( # don't summarize
  post, 
  map = mapping,
  parameters = list(TH_CRCL = 0.0154, TDM_INIT = 0),
  model = pkvancothomson::model(), 
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

## Probability that 400 < AUC < 700
pred_post_full %>%
  filter(t %in% c(36, 48)) %>%
  filter(comp == 3) %>%
  group_by(id) %>%
  tidyr::pivot_wider(names_from = t, values_from = y) %>%
  mutate(auc24 = 2 * (`48` - `36`)) %>%
  ungroup() %>%
  summarise(pta = sum(auc24 > 400 & auc24 < 700)/nrow(.))
  
