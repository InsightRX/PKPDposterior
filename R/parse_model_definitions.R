#' Parse model definitions into a list that can be used to generate a Stan 
#' model
#' 
#' @inheritParams new_stan_model
#' 
parse_model_definitions <- function(
  parameters,
  parameter_definitions,
  covariate_definitions,
  solver
) {
  
  def <- list()
  
  ## 
  

    # // save observations to variables:
    #   ipred_pk = A[2, :] ./ (V1 * mean(WT)); // predictions for all event records
    #   ipred_obs_pk = ipred_pk'[i_obs_pk];          // predictions only for observed data records

  def[["solver_call"]] <- "A = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta);"
  def[["ipred_definition"]] = "ipred_pk = A[2, :] ./ scale;"
  
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
  def[["pk_block"]] <- paste0(
    "theta[j, ", 1:5, "] = ",
    parameter_definitions[unlist(param_req[solver])]
  )
  
  ## Parse covariates
  covariate_def <- NULL
  if(!is.null(covariate_definitions)) {
    def[["covariate_data"]] <- paste0(
      as.character(covariate_definitions),
      " ",
      names(covariate_definitions),
      "[n_t];"
    )
  }
  
  def
  
}