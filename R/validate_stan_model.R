#' Validate the numerical simulation properties of a Stan model vs the model
#' written in PKPDsim.
#' 
#' @param stan_model Stan/Torsten model
#' @param pkpdsim_model PKPDsim implementation
#' 
#' @export
validate_stan_model <- function(
  stan_model,
  pkpdsim_model,
  data,
  parameters = NULL,
  # t_obs = c(3, 6, 24),
  # regimen,
  # covariates = NULL,
  mapping = NULL,
  # dose_cmt = 2,
  n = 50,
  obs_types = c("pk", "pd"),
  max_abs_delta = 1e-3,
  max_rel_delta = 1e-5,
  verbose = FALSE
) {
  
  message("Sampling from posterior...")
  suppressMessages({
    post <- get_mcmc_posterior(
      mod,
      data = data,
      iter_warmup = n,
      iter_sampling = n,
      verbose = verbose
    )
  })
  draws <- as.data.frame(post$draws_df)
  
  # 4. Convert sampled posterior parameters into parameters_table for PKPDsim
  par_stan <- gsub("theta_", "", names(data)[grep("^theta_", names(data))])
  parameters_table <- remap(
    draws[, intersect(par_stan, names(draws))], 
    mapping, 
    reverse = FALSE
  )
  for(key in names(parameters)) { # add fixed parameters
    if(is.null(parameters_table[[key]])) parameters_table[[key]] <- parameters[[key]]
  }
  message("Simulating posterior observations with PKPDsim...")
  
  ## Parse observation types
  obs_type <- rep(NA, length(data$time))
  for(key in obs_types) {
    obs_type[data[[paste0("i_obs_", key)]]] <- key
  }
  t_obs <- data$time[!is.na(obs_type)]
  obs_type <- obs_type[!is.na(obs_type)]
  unq_obs_type <- unique(obs_type)
  obs_type <- match(obs_type, unq_obs_type)
  
  ## Simulate
  simdata <- purrr::map_dfr(1:nrow(parameters_table), function(i) {
    sim(
      ode = pkpdsim_model,
      parameters = as.list(parameters_table[i,]),
      covariates = covariates,
      regimen = regimen,
      only_obs = TRUE,
      t_obs = t_obs,
      obs_type = obs_type,
      output_include = list(variable = TRUE)
    )
  })
  
  # 5. Observations from PKPDsim should be compared to those sampled from Stan. They should be equal.
  message("Comparing Stan and PKPDsim data...")
  comp <- draws[, grep("ipred_obs_", names(draws))] %>% 
    tidyr::pivot_longer(cols = names(.)) %>%
    dplyr::mutate(
      pkpdsim = simdata$y,
      delta = value - pkpdsim,
      delta_rel = ifelse(pkpdsim == 0, delta, delta / pkpdsim),
      idx = rep(seq(t_obs), n)
    )

  for(i in seq(unq_obs_type)) {
    tmp <- comp[comp$idx == i,]
    message("- max absolute delta (", unq_obs_type[i], "): ", max(abs(tmp$delta)))
    message("- max relative delta (", unq_obs_type[i], "): ", max(abs(tmp$delta_rel)))
  }
  if(max(abs(comp$delta)) < max_abs_delta && max(abs(comp$delta_rel)) < max_rel_delta) {
    message(testthat:::praise_emoji(), " ", crayon::green("Validation successful."))
  } else {
    message(crayon::red("Validation not succesful!"))
  }
}
