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
  "pk2cmt", 
  force = TRUE,
  verbose = TRUE
)

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
 
# ## Patient data using PKPDsim-style data:
# data <- prepare_data(
#   regimen = ,
#   covariates = ,
#   data = 
# )

## For now, using pre-created data from Torsten example
source("examples/example_data.R")

## Sample from posterior
chains <- 2
post <- get_mcmc_posterior(
  mod,
  init = prior,
  data = nm_data,
  # prior = prior,
  chains = chains,
  parallel_chains = chains,
  iter_warmup = 500,
  iter_sampling = 500
)

## Plot parameters
par_table <- post$draws_df %>%
  mutate(KA = ka) %>%
  select(KA, CL, V1, Q, V2)
par_table_long <- par_table %>%
  pivot_longer(cols = c(KA,CL, V1, Q, V2))
prior_df <- data.frame(prior[c("ka", "CL", "V1", "Q", "V2")]) %>%
  mutate(KA = ka) %>%
  pivot_longer(cols = c(KA, CL, V1, Q, V2))
ggplot(par_table_long) +
  geom_histogram(aes(x = value)) + 
  facet_wrap(~name, scale = "free") +
  geom_vline(data=prior_df, aes(xintercept = value), colour = 'red') +
  irxreports::theme_irx_minimal()

## Simulate using PKPDsim to get posterior for DV
# recreate the model in PKPDsim
mod1 <- PKPDsim::new_ode_model(
  code = 
  "
  dAdt[0] = -KA*A[0]; \
  dAdt[1] =  KA*A[0] -(CL/V1)*A[1] - (Q/V1)*A[1] + (Q/V2)*A[2]; \
  dAdt[2] = +(Q/V1)*A[1] - (Q/V2)*A[2];\
", 
  obs = list(cmt = 2, scale = "V1"), 
  dose = list(cmt = 1, bioav = 1),
  parameters = c("KA", "CL", "V1", "Q", "V2"),
  covariates = NULL
)
# run simulation
res <- sim(
  ode = mod1,
  parameters_table = as.data.frame(par_table) %>% slice(1:100),
  regimen = PKPDsim::new_regimen(amt = 1000, n = 12, interval = 12),
  t_obs = seq(0, 72, .5),
  only_obs = TRUE
)
ggplot(res, aes(x = t, y = y, group = id)) +
  geom_line(alpha = 0.1) +
  geom_point()
  irxreports::theme_irx_minimal()
