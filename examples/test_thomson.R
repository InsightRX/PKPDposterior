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
library(pkvancothomson)
mod1 <- pkvancothomson::model()

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
  t_obs = seq(0, 48, .5),
  only_obs = FALSE
)

res %>%
  filter(comp == "obs") %>%
  ggplot(aes(x = t, y = y, group = id)) +
  geom_line(alpha = 0.25) +
  geom_point(data = tdm_data, mapping = aes(x = t, y = dv, group = NULL), colour = "red", size=2.5) +
  geom_point(data = tdm_data, mapping = aes(x = t, y = dv, group = NULL), colour = "white", size=1.5) +
  irxreports::theme_irx_minimal()

res %>%
  filter(t %in% c(36, 48)) %>%
  filter(comp == 3) %>%
  group_by(id) %>%
  tidyr::pivot_wider(names_from = t, values_from = y) %>%
  mutate(auc24 = 2 * (`48` - `36`)) %>%
  ggplot() + 
  aes(x = auc24) +
  geom_histogram()
  
