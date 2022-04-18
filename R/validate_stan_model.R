#' Validate the numerical simulation properties of a Stan model vs the model
#' written in PKPDsim.
#' 
#' @param stan_model Stan/Torsten model
#' @param pkpdsim_model PKPDsim implementation
#' @param data dataset (see [PKPDsim_to_stan_data()])
#' @param parameters list of parameter estimates that are used in simulation if 
#'   not already present in either the posterior or the prior draws. I.e. 
#'   specified parameters in this named list will not override parameters with 
#'   the same name in the samples.
#' @param covariates Optional. List of covariates used in PKPDsim model not 
#'   already specified in `data`.
#' @param mapping remap model parameters from Stan to PKPDsim syntax, specified 
#'   as a list. E.g. `map = list(V1 = V)`.
#' @param n number of patient datasets to use in validation. Default is 50.
#' @param max_abs_delta maximum allowed absolute delta between Stan and PKPDsim
#' @param max_rel_delta maximum allowed relative delta between Stan and PKPDsim,
#'   in which PKPDsim is treated as reference.
#' @param verbose verbose output?
#' 
#' @export
validate_stan_model <- function(
  stan_model,
  pkpdsim_model,
  data,
  parameters,
  covariates,
  mapping = NULL,
  n = 50,
  max_abs_delta = 1e-3,
  max_rel_delta = 1e-5,
  verbose = FALSE
) {
  
  if (verbose) message("Sampling from posterior...")
  post <- get_mcmc_posterior(
    stan_model,
    data = data,
    iter_warmup = 200,
    iter_sampling = n,
    verbose = FALSE
  )
  draws <- as.data.frame(post$draws_df)
  
  # Convert sampled posterior parameters into parameters_table for PKPDsim
  par_stan <- gsub("theta_", "", names(data)[grep("^theta_", names(data))])
  parameters_table <- remap(
    draws[, intersect(par_stan, names(draws))], 
    mapping, 
    reverse = FALSE
  )
  for(key in names(parameters)) { # add fixed parameters
    if(is.null(parameters_table[[key]])) parameters_table[[key]] <- parameters[[key]]
  }
  if (verbose) message("Simulating posterior observations with PKPDsim...")
  
  ## Parse observation types
  obs_types <- gsub("n_obs_", "", names(data)[grep("^n_obs_", names(data))])
  obs_type <- rep(NA, length(data$time))
  for(key in obs_types) {
    obs_type[data[[paste0("i_obs_", key)]]] <- key
  }
  t_obs <- data$time[!is.na(obs_type)]
  obs_type <- obs_type[!is.na(obs_type)]
  unq_obs_type <- unique(obs_type)
  obs_type <- match(obs_type, unq_obs_type)
  
  ## Simulate
  covariates_sim <- c(attr(data, "covariates"), covariates)
  simdata <- purrr::map_dfr(1:nrow(parameters_table), function(i) {
    PKPDsim::sim(
      ode = pkpdsim_model,
      parameters = as.list(parameters_table[i,]),
      covariates = covariates_sim,
      regimen = attr(data, "regimen"),
      only_obs = TRUE,
      t_obs = t_obs,
      obs_type = obs_type,
      output_include = list(variable = TRUE)
    )
  })
  
  # Observations from PKPDsim should be compared to those sampled from Stan. They should be equal.
  if (verbose) message("Comparing Stan and PKPDsim data...")
  comp <- draws[, grep("ipred_obs_", names(draws))] %>% 
    tidyr::pivot_longer(cols = dplyr::starts_with("ipred_obs_")) %>%
    dplyr::mutate(
      pkpdsim = simdata$y,
      delta = .data$value - .data$pkpdsim,
      delta_rel = ifelse(.data$pkpdsim == 0, .data$delta, .data$delta / .data$pkpdsim),
      idx = rep(seq(t_obs), n)
    )

  ## Report outcome
  for(i in seq(unq_obs_type)) {
    tmp <- comp[comp$idx == i,]
    message("Stan vs PKPDsim model predictions:")
    message("- max absolute \U{0394} (", unq_obs_type[i], "): ", max(abs(tmp$delta)))
    message("- max relative \U{0394} (", unq_obs_type[i], "): ", max(abs(tmp$delta_rel)))
  }
  if(max(abs(comp$delta)) < max_abs_delta && max(abs(comp$delta_rel)) < max_rel_delta) {
    message(
      praise_emoji(), 
      " ", 
      crayon::green("Validation successful.")
    )
  } else {
    warning("Validation not succesful!")
  }
}

#' Praise using emoji.
#' 
#' Taken from testthat. Not checking for utf8 compatibility.
praise_emoji <- function() {
  sample(c(
    "\U0001f600", # smile
    "\U0001f973", # party face
    "\U0001f638", # cat grin
    "\U0001f308", # rainbow
    "\U0001f947", # gold medal
    "\U0001f389", # party popper
    "\U0001f38a" # confetti ball
  ), 1)
}
