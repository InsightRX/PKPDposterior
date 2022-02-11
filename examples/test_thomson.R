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
mod <- load_model(
  "pk_vanco_thomson", 
  force = TRUE,
  verbose = TRUE
)

## define prior. 
## for our models this should be read from pkvancothomson::parameters() and 
## pkvancothomson::ruv(). Should create a get_prior() fucntion.
prior <- get_init("pkvancothomson")
prior$TH_CRCL <- NULL
prior$TDM_INIT <- NULL

# ## Patient data using PKPDsim-style data:
# data <- prepare_data(
#   regimen = ,
#   covariates = ,
#   data = 
# )

## Sample from posterior
chains <- 2
post <- get_mcmc_posterior(
  mod,
  data = nm_data,
  init = prior,
  chains = chains,
  parallel_chains = chains,
  iter_warmup = 500,
  iter_sampling = 500
  # regimen = reg,
  # covariates = covs,
  # data = tdm_data
)
