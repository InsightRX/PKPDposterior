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
  force = F,
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
   dv = c(40, 12), 
   cmt = c(2, 2)
)
data <- prepare_data(
  regimen, 
  covariates, 
  tdm_data
)

# data <- list(
#   addl = rep(0, 3), 
#   ii = rep(0, 3),
#   amt = c(1000, 0, 0),
#   cmt = c(2, 2, 2),
#   cObs = c(40, 10),
#   evid = c(1, 0, 0),
#   iObs = c(2, 3),
#   nObs = 2,
#   nt = 3,
#   rate = c(500, 0, 0),
#   ss = c(0, 0, 0),
#   time = c(0, 2.5, 11.5),
#   WT = c(70, 70, 70),
#   CRCL = c(5, 5, 5)
# )

## Sample from posterior
chains <- 1
post <- get_mcmc_posterior(
  mod,
  data = data,
  init = prior,
  chains = chains,
  parallel_chains = chains,
  iter_warmup = 1500,
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

## Simulate using PKPDsim to get posterior for DV
# recreate the model in PKPDsim
library(pkvancothomson)
mod1 <- pkvancothomson::model()
# run simulation

par_table <- par_table %>%
  mutate(V = V1, TH_CRCL = 0.0154, TDM_INIT = 0)
covs <- list(
  WT = new_covariate(70),
  CRCL = new_covariate(5),
  CL_HEMO = new_covariate(0)
)
res <- sim(
  ode = mod1,
  parameters_table = as.data.frame(par_table) %>% slice(1:100),
  regimen = regimen,
  covariates = covs,
  t_obs = seq(0, 24, .5),
  only_obs = TRUE
)
ggplot(res, aes(x = t, y = y, group = id)) +
  geom_line(alpha = 0.25) +
  geom_point(data = tdm_data, mapping = aes(x = t, y = dv, group = NULL), colour = "red", size=2.5) +
  geom_point(data = tdm_data, mapping = aes(x = t, y = dv, group = NULL), colour = "white", size=1.5) +
  irxreports::theme_irx_minimal()
