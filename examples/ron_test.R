library(PKPDsim)
library(pkpdmcmctbd)
library(ggplot2)
library(dplyr)
library(tidyr)

# set Stan path
Sys.setenv(STAN_PATH = "/home/ron@insight-rx.com/git/Torsten")

# Compile or reload model
# this creates a binary in the installed package folder
# at ~/R/x86_64-pc-linux-gnu-library/4.1/pkpdmcmctbd/models/
mod <- load_model("pk2cmt", verbose = FALSE)

## define prior. 
## for our models this should be read from pkvancothomson::parameters() and 
## pkvancothomson::ruv(). Should create a get_prior() fucntion.
prior <- list( 
  CL = 7.43,
  ka = 1.08,
  Q = 28.0,
  sigma = 0.58,
  V1 = 78.4,
  V2 = 68.1
)
 
## Patient data using PKPDsim-style workflow:
# reg <- new_regimen()
# covs <- list(
#   WT = ..,
#   CRCL = ...
# )

## For now, using pre-created data from Torsten example
source("examples/example_data.R")

## Sample from posterior
chains <- 2
post <- get_mcmc_posterior(
  mod,
  nm_data = nm_data,
  prior = prior,
  chains = chains,
  parallel_chains = chains,
  iter_warmup = 500,
  iter_sampling = 500
  # regimen = reg,
  # covariates = covs,
  # data = tdm_data
)

## Plot parameters
par_table <- post$draws_df %>%
  select(CL, V1, Q, V2) %>%
  pivot_longer(cols = c(CL, V1, Q, V2))
prior_df <- data.frame(prior[c("CL", "V1", "Q", "V2")]) %>%
  pivot_longer(cols = c(CL, V1, Q, V2))
ggplot(par_table) +
  geom_histogram(aes(x = value)) + 
  facet_wrap(~name, scale = "free") +
  geom_vline(data=prior_df, aes(xintercept = value), colour = 'red') +
  irxreports::theme_irx_minimal()

## Simulate using PKPDsim to get posterior for DV

