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
prior$V1 <- prior$V
prior$V <- NULL
prior$TH_CRCL <- NULL
prior$TDM_INIT <- NULL

data <- list(
  addl = rep(0, 3), 
  ii = rep(0, 3),
  amt = c(1000, 0, 0),
  cmt = c(2, 2, 2),
  cObs = c(40, 10),
  evid = c(1, 0, 0),
  iObs = c(2, 3),
  nObs = 2,
  nt = 3,
  rate = c(500, 0, 0),
  ss = c(0, 0, 0),
  time = c(0, 2.5, 11.5),
  WT = c(70, 70, 70),
  CRCL = c(5, 5, 5)
)

## Sample from posterior
chains <- 1
post <- get_mcmc_posterior(
  mod,
  data = data,
  init = prior,
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
  select(CL, V1, Q, V2)
par_table_long <- par_table %>%
  pivot_longer(cols = c(CL, V1, Q, V2))
prior_df <- data.frame(prior[c("CL", "V1", "Q", "V2")]) %>%
  pivot_longer(cols = c(CL, V1, Q, V2))
ggplot(par_table_long) +
  geom_histogram(aes(x = value)) + 
  facet_wrap(~name, scale = "free") +
  geom_vline(data=prior_df, aes(xintercept = value), colour = 'red') +
  irxreports::theme_irx_minimal()

