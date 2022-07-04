## Neutropenia PK-PD example
## Implemented using ODE system for both PK and PD

library(PKPDsim)
library(PKPDposterior)
library(pkpdneutropeniatemplate1)
library(dplyr)
library(ggplot2)

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
data <- new_stan_data(
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
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.95,
  verbose = TRUE
)
# saveRDS(post, "post.rds")
# post <- readRDS("post.rds")

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
  variable = "ANC",
  summarize = T
)
plot_predictions(pred, obs = pd_data) +
  xlim(c(0.1, 300)) +
  scale_y_log10(limits = c(0.3, 9)) +
  vpc::theme_empty()

comb_data$obs_type <- c(1,1,2,2,2,2)
library(pkpdneutropeniatemplate1)
indiv_est <- PKPDmap::get_map_estimates(
  model = pkpdneutropeniatemplate1::model(),
  data = comb_data %>% mutate(y = dv),
  parameters = pkpdneutropeniatemplate1::parameters(),
  covariates = covariates,
  omega = omega,
  error = list(add = c(0.2, 0.3), prop = c(0, 0)),
  regimen = regimen,
  obs_type_label = "obs_type",
  fixed = c("KA", "Q", "V2", "GAMMA"),
  verbose = T
)
res <- sim(
  pkpdneutropeniatemplate1::model(),
  omega = indiv_est$vcov_full,
  regimen = regimen,
  n = 200,
  t_obs = seq(0, 15*24, 6),
  parameters = indiv_est$parameters,
  covariates = covariates,
  output_include = list(variables = TRUE)
) %>%
  mutate(y = ANC)
ci <- c(0.05, 0.95)
pred_map <- res %>% 
  dplyr::filter(.data$comp == "obs") %>%
  dplyr::group_by(.data$t) %>%
  dplyr::summarise(
    ymedian = stats::median(.data$y, na.rm = TRUE), 
    ymin = stats::quantile(.data$y, ci[1], na.rm = TRUE),
    ymax = stats::quantile(.data$y, ci[2], na.rm = TRUE)
  ) %>%
  mutate(type = "map")
pred_comb <- bind_rows(
  pred %>% mutate(type = "hmc"), 
  pred_map
)
obs <- bind_rows(pd_data %>% mutate(type = "hmc"), 
                 pd_data %>% mutate(type = "map"))
plot_predictions(pred_comb, obs = obs) +
  xlim(c(0.1, 300)) +
  scale_y_log10(limits = c(0.3, 9)) +
  vpc::theme_empty() +
  facet_wrap(~type)

