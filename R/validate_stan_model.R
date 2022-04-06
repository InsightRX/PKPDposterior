#' Validate the numerical simulation properties of a Stan model vs the model
#' written in PKPDsim.
#' 
#' @param stan_model Stan/Torsten model
#' @param pkpdsim_model PKPDsim implementation
validate_stan_model <- function(
  stan_model,
  pkpdsim_model,
  parameters,
  t_obs = c(3, 6, 24),
  regimen,
  covariates = NULL,
  mapping = NULL,
  dose_cmt = 2,
  n = 50,
  max_abs_delta = 1e-3,
  max_rel_delta = 1e-5,
  verbose = FALSE
) {
  
  # 1. simulate some dummy data in PKPDsim using a certain sampling schema. Needs to be reasonable data (so the sampler doesn’t crash), but otherwise can be anything.
  message("Simulating dummy data...")
  dummy_data <- PKPDsim::sim(
    pkpdsim_model,
    parameters = parameters,
    regimen = regimen,
    covariates = covariates,
    t_obs = t_obs,
    only_obs = TRUE,
    output_include = list()
  )
  
  # 2. pass this data into Stan and sample, e.g. 100 times. 
  #    We should extract both the posterior estimates for the PK parameters, 
  #    as well as the predictions for the samples.
  par_stan <- remap(parameters, mapping, reverse = TRUE)
  iiv <- par_stan
  for(key in names(iiv)) iiv[[key]] <- 0.2
  data <- prepare_data(
    regimen,
    covariates, 
    data = dummy_data %>% 
      dplyr::mutate(dv = y) %>%
      dplyr::select(t, dv),
    parameters = par_stan,
    iiv = iiv,
    ruv = list(
      prop = 0.1,
      add = 0.1
    ),
    dose_cmt = dose_cmt
  )
  
  # 3. Sample from the posterior using the prepared data
  message("Sampling from posterior...")
  suppressMessages({
    post <- get_mcmc_posterior(
      mod,
      data = data,
      iter_warmup = 500,
      iter_sampling = n,
      verbose = verbose
    )
  })
  draws <- as.data.frame(post$draws_df)
  
  # 4. The posterior parameter estimates can be fed into PKPDsim (use parameters_table), and used to simulate observations. 
  parameters_table <- remap(
    draws[, intersect(names(par_stan), names(draws))],
    mapping
  )
  for(key in names(parameters)) { # add fixed parameters
    if(is.null(parameters_table[[key]])) parameters_table[[key]] <- parameters[[key]]
  }
  message("Simulating posterior observations with PKPDsim...")
  simdata <- purrr::map_dfr(1:nrow(parameters_table), function(i) {
    sim(
      ode = pkpdsim_model,
      parameters = as.list(parameters_table[i,]),
      covariates = covariates,
      regimen = regimen,
      only_obs = TRUE,
      t_obs = t_obs
    )  
  })
  
  # 5. Observations from PKPDsim should be compared to those sampled from Stan. They should be equal.
  message("Comparing Stan and PKPDsim data...")
  comp <- draws[, grep("ipred_obs_pk", names(draws))] %>% 
    tidyr::pivot_longer(cols = names(.)) %>%
    dplyr::mutate(
      pkpdsim = simdata$y,
      delta = value - pkpdsim,
      delta_rel = ifelse(pkpdsim == 0, delta, delta / pkpdsim)
    )
  
  if(max(abs(comp$delta)) < max_abs_delta && max(abs(comp$delta_rel)) < max_rel_delta) {
    message("Succesful validation.")
  } else {
    message("Validation unsuccesful!")
  }
  message("Max absolute delta: ", max(abs(comp$delta)))
  message("Max relative delta: ", max(abs(comp$delta_rel)))
}
