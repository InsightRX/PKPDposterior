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
  force = T,
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
   dv = c(40, 14), 
   cmt = c(2, 2)
)
data <- prepare_data(
  regimen, 
  covariates, 
  tdm_data
)

## Sample from posterior
post <- get_mcmc_posterior(
  mod,
  data = data,
  init = prior,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.95
)

#############################################################################
## Plots
#############################################################################

## Plot parameters
par_table <- post$draws_df              %>% select(CL, V1, Q, V2)
par_table_prior <- post$draws_df        %>% select(prior_CL, prior_V1, prior_Q, prior_V2)
names(par_table_prior) <- gsub("prior_", "", names(par_table_prior))
par_table_long <- par_table             %>% pivot_longer(cols = c(CL, V1, Q, V2))
par_table_prior_long <- par_table_prior %>% pivot_longer(cols = c(CL, V1, Q, V2))
ggplot(par_table_long) +
  geom_histogram(data = par_table_prior_long, aes(x = value), alpha = 0.2, fill="blue") + 
  geom_histogram(aes(x = value), alpha = 0.3, fill = "red") + 
  facet_wrap(~name, scale = "free") +
  irxreports::theme_irx_minimal()

## Simulate using PKPDsim to get posterior for concentration and AUC
library(pkvancothomson)
mod1 <- pkvancothomson::model()

par_table <- par_table %>% mutate(V = V1, TH_CRCL = 0.0154, TDM_INIT = 0)
par_table_prior <- par_table_prior %>% mutate(V = V1, TH_CRCL = 0.0154, TDM_INIT = 0)
covs <- list(
  WT = new_covariate(70),
  CRCL = new_covariate(5),
  CL_HEMO = new_covariate(0)
)
res <- sim(
  ode = mod1,
  parameters_table = as.data.frame(par_table), # %>% slice(1:100),
  regimen = regimen,
  covariates = covs,
  t_obs = seq(0, 48, .5),
  only_obs = FALSE
)
res_prior <- sim(
  ode = mod1,
  parameters_table = as.data.frame(par_table_prior), # %>% slice(1:100),
  regimen = regimen,
  covariates = covs,
  t_obs = seq(0, 48, .5),
  only_obs = FALSE
)

## Plot all posterior samples
res %>%
  filter(comp == "obs") %>%
  ggplot(aes(x = t, y = y, group = id)) +
    geom_line(alpha = 0.25) +
    geom_point(data = tdm_data, mapping = aes(x = t, y = dv, group = NULL), colour = "red", size=2.5) +
    geom_point(data = tdm_data, mapping = aes(x = t, y = dv, group = NULL), colour = "white", size=1.5) +
    irxreports::theme_irx_minimal()

## Plot confidence interval prior
res_prior %>%
  filter(comp == "obs") %>%
  group_by(t) %>%
  summarise(ymedian = median(y), 
            ymin = quantile(y, 0.05),
            ymax = quantile(y, 0.95)) %>%
  ggplot(aes(x = t, y = ymedian)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#cfcfcf") +
  geom_line(size = 1) +
  geom_point(data = tdm_data, mapping = aes(x = t, y = dv, group = NULL), colour = "red", size=2.5) +
  geom_point(data = tdm_data, mapping = aes(x = t, y = dv, group = NULL), colour = "white", size=1.5) +
  irxreports::theme_irx_minimal()

## Plot confidence interval posterior
res %>%
  filter(comp == "obs") %>%
  group_by(t) %>%
  summarise(ymedian = median(y), 
            ymin = quantile(y, 0.05),
            ymax = quantile(y, 0.95)) %>%
  ggplot(aes(x = t, y = ymedian)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#cfcfcf") +
    geom_line(size = 1) +
    geom_point(data = tdm_data, mapping = aes(x = t, y = dv, group = NULL), colour = "red", size=2.5) +
    geom_point(data = tdm_data, mapping = aes(x = t, y = dv, group = NULL), colour = "white", size=1.5) +
    irxreports::theme_irx_minimal()

## Plot AUC dist
res %>%
  filter(t %in% c(36, 48)) %>%
  filter(comp == 3) %>%
  group_by(id) %>%
  tidyr::pivot_wider(names_from = t, values_from = y) %>%
  mutate(auc24 = 2 * (`48` - `36`)) %>%
  ggplot() + 
  aes(x = auc24) +
  geom_histogram()

## Probability that 400 < AUC < 700
res %>%
  filter(t %in% c(36, 48)) %>%
  filter(comp == 3) %>%
  group_by(id) %>%
  tidyr::pivot_wider(names_from = t, values_from = y) %>%
  mutate(auc24 = 2 * (`48` - `36`)) %>%
  ungroup() %>%
  summarise(pta = sum(auc24 > 400 & auc24 < 700)/nrow(.))
  
