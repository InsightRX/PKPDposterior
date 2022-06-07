## Check what optimal default settings are

library(PKPDposterior)
library(PKPDsim)
library(pkvancothomson)
library(dplyr)
library(ggplot2)

## define init values (use population values): 
prior <- prior_from_PKPDsim_model(
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

## Validate model vs PKPDsim implementation
validate_stan_model(
  stan_model = mod,
  pkpdsim_model = pkpdneutropeniatemplate1::model(),
  parameters = pkpdneutropeniatemplate1::parameters(),
  data = data,
  mapping = mapping
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
## 1. number of threads
warmup <- c(100, 200, 300, 500)
sampling <- c(300, 500, 800, 1000)
## Reference sampling:
ref_post <- get_mcmc_posterior(
  mod = mod,
  data = data,
  iter_warmup = 2000,
  iter_sampling = 2000,
  threads = 4,
  adapt_delta = 0.95,
  verbose = TRUE
)

reps <- 5
dat <- data.frame()
for(i in seq(warmup)) {
  for(j in seq(sampling)) {
    for(k in 1:reps) {
      post <- get_mcmc_posterior(
        mod = mod,
        data = data,
        iter_warmup = warmup[i],
        iter_sampling = sampling[j],
        adapt_delta = 0.95,
        chains = 1,
        seed = round(runif(1, 1, 982312))
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

dat %>%
  group_by(name, warmup, sampling) %>%
  summarise(mean_rse = mean(rse), sd_rse = sd(rse)) %>%
  ggplot(aes(x = sampling, y = mean_rse, colour = as.factor(warmup))) +
    geom_line() +
    geom_point() +
    facet_wrap(~ name)

dat %>%
  group_by(name, warmup, sampling) %>%
  summarise(mean_rse = mean(rse), sd_rse = sd(rse)) %>%
  ggplot(aes(x = sampling, y = sd_rse, colour = as.factor(warmup))) +
    geom_line() +
    geom_point() +
    facet_wrap(~ name)

## Repeat for number of chains = 4

reps <- 5
dat2 <- data.frame()
for(i in seq(warmup)) {
  for(j in seq(sampling)) {
    for(k in 1:reps) {
      post <- get_mcmc_posterior(
        mod = mod,
        data = data,
        iter_warmup = warmup[i],
        iter_sampling = sampling[j],
        adapt_delta = 0.95,
        chains = 4, # !!!
        seed = round(runif(1, 1, 982312))
      )
      dat2 <- bind_rows(
        dat2,
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

dat2 %>%
  group_by(name, warmup, sampling) %>%
  summarise(mean_rse = mean(rse), sd_rse = sd(rse)) %>%
  ggplot(aes(x = sampling, y = sd_rse, colour = as.factor(warmup))) +
  geom_line() +
  geom_point() +
  facet_wrap(~ name)

## 5. adapt_delta
reps <- 5
dat3 <- data.frame()
for(i in seq(warmup)) {
  for(j in seq(sampling)) {
    for(k in 1:reps) {
      post <- get_mcmc_posterior(
        mod = mod,
        data = data,
        iter_warmup = warmup[i],
        iter_sampling = sampling[j],
        adapt_delta = 0.8, # !!!
        chains = 1,
        seed = round(runif(1, 1, 982312))
      )
      dat3 <- bind_rows(
        dat3,
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

dat3 %>%
  group_by(name, warmup, sampling) %>%
  summarise(mean_rse = mean(rse), sd_rse = sd(rse)) %>%
  ggplot(aes(x = sampling, y = sd_rse, colour = as.factor(warmup))) +
  geom_line() +
  geom_point() +
  facet_wrap(~ name)
