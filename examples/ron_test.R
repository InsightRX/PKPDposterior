library(PKPDsim)
library(pkpdmcmctbd)

# set Stan path
Sys.setenv(STAN_PATH = "/home/ron@insight-rx.com/git/Torsten")

# Compile or reload model
# this creates a binary in the installed package folder
# at ~/R/x86_64-pc-linux-gnu-library/4.1/pkpdmcmctbd/models/
mod <- load_model("pk2cmt", verbose = FALSE)

## Patient data
reg <- new_regimen()
covs <- list(
  WT = ..,
  CRCL = ...
)
tdm_data <- data.frame()

## Sample from posterior
post <- get_mcmc_posteriori(
  mod,
  regimen = reg,
  covariates = covs,
  data = tdm_data
)

## Plotting

