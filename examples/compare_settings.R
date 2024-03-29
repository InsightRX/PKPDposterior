## Check what optimal default settings are

library(PKPDposterior)
library(PKPDsim)
library(pkvancothomson)
calc_stats_chains <- function(post) {
  post$draws_df %>%
    group_by(.chain) %>%
    select(!!names(parameters)) %>%
    tidyr::pivot_longer(cols = !!names(parameters)) %>%
    group_by(.chain, name) %>%
    summarise(value = mean(value)) %>%
    group_by(name) %>%
    summarise(rse = sd(value) / mean(value), max_dev = max(value - mean(value)))
}
calc_mean <- function(post) {
  as.data.frame(post$draws_df) %>%
    select(!!names(parameters)) %>%
    tidyr::pivot_longer(cols = !!names(parameters)) %>%
    group_by(name) %>%
    summarise(value = mean(value))
}
calc_stats <- function(post, ref_post) {
  mean_ref <- calc_mean(ref_post)
  data.frame(
    name = mean_ref$name,
    rse = (calc_mean(post)$value - mean_ref$value) / mean_ref$value
  )
}

## Vancomycin model (Thomson et al)
parameters <- list(CL = 2.99, Q = 2.28, V2 = 0.732, V1 = 0.675)
params_pkpdsim <- list(CL = 2.99, Q = 2.28, V2 = 0.732, V1 = 0.675, TH_CRCL = 0, TDM_INIT = 0)
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
regimen <- new_regimen(
  amt = 1500, 
  n = 4, 
  times = c(0, 12, 24, 36), 
  type = 'infusion',
  t_inf = 2
)
covariates <- list(
  WT = new_covariate(value = 70, unit = "kg"),
  CRCL = new_covariate(value = 5, unit = "l/hr"),
  CL_HEMO = new_covariate(0)
)
tdm_data <- sim(
  pkvancothomson::model(),
  regimen = regimen,
  parameters = remap(parameters, mapping),
  covariates = covariates,
  t_obs = c(2.5, 11.5),
  only_obs = TRUE
) %>% select(t, dv = y)

## Create combined dataset for Torsten/Stan to read:
data <- new_stan_data(
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

## 1. number of threads
dat <- data.frame()
warmup <- c(100, 200, 300, 500)
sampling <- c(300, 500, 800, 1000)
ref_post <- get_mcmc_posterior(
  mod = mod,
  data = data,
  iter_warmup = 2000,
  iter_sampling = 5000,
  adapt_delta = 0.95,
  chains = 1
)
reps <- 5
for(i in seq(warmup)) {
  for(j in seq(sampling)) {
    for(k in 1:reps) {
      post <- get_mcmc_posterior(
        mod = mod,
        data = data,
        iter_warmup = warmup[i],
        iter_sampling = sampling[j],
        adapt_delta = 0.95,
        chains = 1
      )
      dat <- bind_rows(
        dat,
        data.frame(
          calc_stats(post, ref_post),
          warmup = warmup[i], 
          sampling = sampling[j],
          iter = k
        )
      )
    }
  }
}

ggplot(dat, aes(x = sampling, y = rse, colour = as.factor(warmup))) +
  geom_line() +
  facet_wrap(~ name)


## 2. number of warmup

## 3. number of samples

## 4. thinning

## 5. adapt_delta
