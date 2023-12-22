#' Parse model definitions into a list that can be used to generate a Stan 
#' model
#' 
#' @inheritParams new_stan_model
#' @keywords internal
parse_model_definitions <- function(
  parameters,
  parameter_definitions,
  fixed,
  variable_definitions,
  ode,
  covariate_definitions,
  solver,
  obs_types,
  obs_cmt,
  scale,
  custom_ipred
) {
  
  def <- list()
  cr <- "\n"
  n_cmt <- get_n_cmt(solver, ode)
  n_theta <- length(parameter_definitions)
  
  if(!is.null(ode)) {
    def[["ode_function"]] <- new_ode_function(
      ode = ode,
      parameter_names = names(parameter_definitions),
      n_cmt = n_cmt
    )  
  }
  
  ## Define observation variables and data
  def[["observation_variables"]] <- c(
    "int<lower = 1> n_t;  // number of events",
    paste0("int<lower = 1> n_obs_", obs_types, ";  // number of observation"),
    paste0("int<lower = 1> i_obs_", obs_types, "[n_obs_", obs_types, "];  // index of observation")
  )
  def[["observation_data"]] <- paste0(
    "vector<lower = 0>[n_obs_", obs_types, "] dv_", obs_types, ";", 
    "  // observed concentration (Dependent Variable)"
  )

  ## define model parameters in data section
  # parameters_sampled <- setdiff(parameters, fixed)
  if(!is.null(fixed)) { # fixed parameters need to be defined differently, to avoid being sampled
    parameters_sampled <- parameters[! names(parameters) %in% fixed]
  } else {
    parameters_sampled <- parameters
  }
  def[["population_parameters"]] <- paste0("real theta_", names(parameters), ";")
  def[["iiv_parameters"]] <- paste0("real<lower=0> omega_", names(parameters_sampled), ";")
  def[["error_parameters"]] <- paste0(
    "int<lower = 0, upper = 1> ltbs_", obs_types, ";", cr, "  ",
    "real<lower = 0> ruv_prop_", obs_types, ";", cr, "  ",
    "real<lower = 0> ruv_add_", obs_types, ";"
  )
  
  ## Input data, hardcoded
  def[["input_data"]] <- c(
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
  def[["log_transform_observations"]] <- paste0(
    "vector[n_obs_", obs_types, "] log_dv_", obs_types, " = log(dv_", obs_types, ");"
  )
  theta_size <- switch(
    solver,
    "pmx_solve_onecpt" = 3,
    "pmx_solve_twocpt" = 5
  )
  if(!is.null(theta_size)) {
    if(n_theta != theta_size) {
      stop("Length of `parameter_definitions` vector different from expected length based on chosen analytic solver. Please revise `parameter_definitions` or change solver.")
    } 
  }
  if(is.null(n_theta)) {
    stop("Please specify size of theta vector using `n_theta`.")
  }
  def[["model_numbers"]] <- c(
    paste0("int n_theta = ", n_theta, ";"),
    paste0("int n_cmt = ", n_cmt,";")
  )
  
  ## parameters block
  def[["parameter_definitions"]] <- paste0("real<lower=0> ", names(parameters_sampled), ";")

  ## Transformed parameters block:
  solver_args <- switch(
    solver,
    "pmx_solve_onecpt" = "time, amt, rate, ii, evid, cmt, addl, ss, theta",
    "pmx_solve_twocpt" = "time, amt, rate, ii, evid, cmt, addl, ss, theta",
    "pmx_solve_rk45" = "ode, n_cmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5",
    "pmx_solve_adams" = "ode, n_cmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5",
    "pmx_solve_bdf" = "ode, n_cmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5"
  )
  def[["solver_call"]] <- paste0("A = ", solver, "(", solver_args, ");")
  if(!is.null(custom_ipred)) { # custom definition of individual predictions
    ipred <- paste0("ipred_", obs_types, " = ", custom_ipred[obs_types], ";")
  } else {
    ipred <- paste0("ipred_", obs_types, " = A[", obs_cmt[obs_types], ", :] ./ ", scale[obs_types], ";")
  }
  def[["ipred_definition"]] = c(
    ipred,
    paste0(
      "ipred_obs_", obs_types, " = ipred_", obs_types, "'[i_obs_", obs_types, "];", 
      " // predictions only for observed data records"
    )
  )

  ## Parse parameters for analytic equations
  param_req <- list(
    "pmx_solve_onecpt" = c("CL", "V", "KA"),
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
  
  if(!is.null(variable_definitions)) {
    def[["variable_definitions"]] <- paste0("real ", names(variable_definitions), ";")
  }
  
  def[["transformed_parameters_definitions"]] <- c(
    paste0("row_vector<lower = 0>[n_t] ipred_", obs_types, ";"),
    paste0("vector<lower = 0>[n_obs_", obs_types, "] ipred_obs_", obs_types, ";"),
    "matrix[n_cmt, n_t] A;"
  )
  
  if(!is.null(fixed)) {
    def[["fixed_parameters"]] <- paste0(
      "real ", fixed, " = theta_", fixed, ";"
    )
  }

  if(is.null(ode)) {
    # make sure to have the order correct for the analytic solver to understand
    pk_equations <- parameter_definitions[unlist(param_req[solver])] 
    variable_definitions_code <- NULL
    if(!is.null(variable_definitions)) {
      variable_definitions_code <- paste0(
        "  ", names(variable_definitions), " = ", as.character(variable_definitions), ";" 
      )
    }
    def[["pk_block"]] <- c(
      "array[n_t, n_theta] real theta;",
      "for(j in 1:n_t) {",
      variable_definitions_code,
      paste0(
        "  theta[j, ", seq(pk_equations), "] = ",
        pk_equations,
        ";"
      ),
      "}"
    )
  } else {
    pk_equations <- parameter_definitions # for ODEs the user defines the order
    def[["pk_block"]] <- c(
      "real theta[n_theta];",
      "for(j in 1:n_t) {",
      paste0(
        "  theta[", seq(pk_equations), "] = ",
        pk_equations,
        ";"
      ),
      "}"
    )
  }
  
  def[["likelihood_parameters"]] <- paste0(
    names(parameters_sampled), " ~ lognormal(log(theta_", 
    names(parameters_sampled), "), omega_", 
    names(parameters_sampled), ");"
  )

  def[["likelihood_observed_data"]] <- paste0(
    "if(ltbs_", obs_types, ") {", cr, "  ",
    "  log_dv_", obs_types, " ~ normal(log(ipred_obs_", obs_types, "), ruv_add_", obs_types, "); ", cr, "  ",
    "} else {", cr, "  ",
    "  dv_", obs_types, " ~ normal(ipred_obs_", obs_types, ", (sqrt(pow(ruv_prop_", obs_types, " * ipred_obs_", obs_types, ", 2) + pow(ruv_add_", obs_types,", 2)));", cr, "  ",
    "}"
  )
  
  def[["sample_prior"]] <- paste0(
    "real prior_", names(parameters_sampled),
    " = lognormal_rng(log(theta_", names(parameters_sampled), 
    "), omega_", names(parameters_sampled), ");"
  )

  def[["simulate_posterior_ruv"]] <- c(
    paste0("real ipred_ruv_", obs_types, "[n_obs_", obs_types, "];"),
    paste0(
       "for(i in 1:n_obs_", obs_types, "){", cr, "  ",
       "  ipred_ruv_", obs_types, "[i] = normal_rng(ipred_obs_", obs_types, "[i], (sqrt(pow(ruv_prop_", obs_types, " * ipred_obs_", obs_types, "[i], 2) + pow(ruv_add_", obs_types, ", 2)));", cr, "  ",
       "}"
    ) 
  )

  def
  
}

#' Get the number of compartments, either from analytical solver definition,
#' or from the dAdt definitions in the ODE block.
#' 
#' @inheritParams new_stan_model
#' @keywords internal
get_n_cmt <- function(solver, ode = NULL) {
  if(is.null(ode)) {
    if(solver == "pmx_solve_twocpt") {
      n_cmt <- 3
    } else {
      n_cmt <- 2
    }
  } else {
    ode_string <- paste0(ode, collapse="")
    dadt_matches <- regmatches(
      ode_string,
      gregexpr("dAdt\\[.\\]", ode_string)
    )[[1]]
    cmts <- as.numeric(gsub("\\D+", "", dadt_matches))
    n_cmt <- max(cmts, na.rm = TRUE)
  }
  n_cmt
}
