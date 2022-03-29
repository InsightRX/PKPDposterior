#' Parse model definitions into a list that can be used to generate a Stan 
#' model
#' 
#' @inheritParams new_stan_model
#' 
parse_model_definitions <- function(
  parameters,
  parameter_definitions,
  covariate_definitions,
  solver,
  obs_types,
  n_theta,
  n_cmt,
  scale
) {
  
  def <- list()
  cr <- "\n"
  
  ## Define observation variables and data
  obs_types_txt <- c()
  for(type in obs_types) {
    obs_types_txt <- c(
      obs_types_txt,
      "int<lower = 1> n_t;  // number of events",
      paste0("int<lower = 1> n_obs_", type, ";  // number of observation"),
      paste0("int<lower = 1> i_obs_", type, "[n_obs_", type, "];  // index of observation")
    )
  }
  def[["observation_variables"]] <- obs_types_txt
  def[["observation_data"]] <- paste0("vector<lower = 0>[n_obs_", obs_types, "] dv_", obs_types, ";  // observed concentration (Dependent Variable)")

  ## define model parameters in data section
  def[["population_parameters"]] <- paste0("real<lower=0> theta_", names(parameters), ";")
  def[["iiv_parameters"]] <- paste0("real<lower=0> omega_", names(parameters), ";")
  def[["error_parameters"]] <- paste0(
    "int<lower = 0, upper = 1> ltbs_", obs_types, ";", cr, "  ",
    "real<lower = 0> ruv_prop_", obs_types, ";", cr, "  ",
    "real<lower = 0> ruv_add_", obs_types, ";"
  )
  
  ## Input data, hardcoded
  def[["input_data"]] <- c(
    "input_data",
    "int<lower = 1> cmt[n_t];",
    "int  evid[n_t];",
    "int  addl[n_t];",
    "int  ss[n_t];",
    "real amt[n_t];",
    "real time[n_t];", 
    "real rate[n_t];",
    "real ii[n_t];"
  )
  
  ## Parse covariates
  if(!is.null(covariate_definitions)) {
    def[["covariate_data"]] <- paste0(
      as.character(covariate_definitions),
      " ",
      names(covariate_definitions),
      "[n_t];"
    )
  }

  ## Transformed data block:
  def[["log_transform_observations"]] <- paste0("vector[n_obs_", obs_types, "] log_dv_", obs_types, " = log(dv_", obs_types, ");")
  cmt_size <- list(
    "pmx_solve_onecpt" = 2,
    "pmx_solve_twocpt" = 3,
    "pmx_solve_threecpt" = 4
  )
  theta_size <- list(
    "pmx_solve_onecpt" = 3,
    "pmx_solve_twocpt" = 5,
    "pmx_solve_threecpt" = 5
  )
  if(!is.null(theta_size[[solver]])) {
    n_theta <- theta_size[[solver]]
  }
  if(!is.null(cmt_size[[solver]])) {
    n_cmt <- cmt_size[[solver]]
  }
  if(is.null(n_theta)) {
    stop("Please specify size of theta vector using `n_theta`.")
  }
  if(is.null(n_cmt)) {
    stop("Please specify size of compartment vector using `n_cmt`.")
  }
  def[["model_numbers"]] <- c(
    paste0("int n_theta = ", n_theta, ";"),
    paste0("int n_cmt = ", n_cmt,";")
  )
  
  ## parameters block
  def[["parameter_definitions"]] <- paste0("real<lower=0> ", names(parameters))

  ## Transformed parameters block:
  def[["solver_call"]] <- "A = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta);"
  def[["ipred_definition"]] = paste0("ipred_", obs_types, " = A[2, :] ./ ", scale[obs_types], ";")
  
  ## Parse parameters for analytic equations
  param_req <- list(
    "pmx_solve_twocpt" = c("CL", "Q", "V1", "V2", "KA")
  )
  if(solver %in% names(param_req)) { # check all parameters specified, if analytic model
    if(! all(param_req[[solver]] %in% names(parameter_definitions))) {
      stop(paste0(
        "Please specify parameter definitions for all parameters in this model: ", 
        paste0(param_req[[solver]], collapse=", "), "."
      ))
    }
  }
  
  def[["transformed_parameters_definitions"]] <- c(
    paste0("row_vector<lower = 0>[n_t] ipred_", obs_types, ";"),
    paste0("vector<lower = 0>[n_obs_", obs_types, "] ipred_obs_", obs_types, ";"),
    "matrix<lower = 0>[n_cmt, n_t] A;"
  )
  
  def[["pk_block"]] <- c(
    "array[n_t, n_theta] real theta;",
    paste0(
      "theta[j, ", 1:5, "] = ",
      parameter_definitions[unlist(param_req[solver])]
    )
  )
  
  def[["likelihood_parameters"]] <- paste(names(parameters), " ~ lognormal(log(theta_", names(parameters), ", omega_", names(parameters), ");")

  def[["likelihood_observed_data"]] <- paste0(
    "if(ltbs_", obs_types, ") {", cr, "  ",
    "  log_dv_", obs_types, " ~ normal(log(ipred_obs_", obs_types, "), ruv_add_", obs_types, "); ", cr, "  ",
    "} else {", cr, "  ",
    "  dv_", obs_types, " ~ normal(ipred_obs_", obs_types, ", (ruv_prop_", obs_types, " * ipred_obs_", obs_types, " + ruv_add_", obs_types,"));", cr, "  ",
    "}"
  )
  
  def[["sample_prior"]] <- paste("real prior_", names(parameters), " = lognormal_rng(log(theta_", names(parameters), ", omega_", names(parameters), ");")

  def[["simulate_posterior_ruv"]] <- c(
    "real ipred_ruv_pk[n_obs_pk];",
    paste0(
       "for(i in 1:n_obs_", obs_types, "){", cr, "  ",
       "  ipred_ruv_", obs_types, "[i] = normal_rng(ipred_obs_", obs_types, "[i], (ruv_prop_", obs_types, " * ipred_obs_", obs_types, "[i] + ruv_add_", obs_types, "));", cr, "  ",
       "}"
    ) 
  )

  def
  
}