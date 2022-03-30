#' Create a new stan model using a template
#' 
#' @param parameters list of parameter values, e.g. `list(CL = 5, V = 50)`.
#' @param parameter_definitions list of parameter definitions. See example 
#' below.
#' @param ode code for differential equation system, if not the analytic solver
#' is used.
#' @param covariate_definitions list of covariate definitions. See example
#' below.
#' @param solver Torsten solver to use. Starts with `pmx_solve_`.
#' @param template internal Stan template file to use. Defaults to 
#' `template.stan`
#' @param obs_types vector of observed data types. If `NULL`, will use `"pk"`.
#' @param obs_cmt observation compartment, or list of observation compartments 
#' for each observation type when multiple observation types. E.g. `2` or 
#' `list(pk = 2, pd = 8)`.
#' @param scale scale definition. Can be either parameter expression or number.
#' @param custom_ipred list of custom definitions for individual predictions for 
#' each observation type. Sometimes this is needed when more flexibility is 
#' needed than the `scale` parameter provides.
#' @param n_cmt number of compartments. If analytic equation, will use number 
#' of compartments for solver.
#' @param n_theta number of thetas. If analytic equation, will use number of 
#' thetas for solver.
#' @param file file to save model to. If NULL, will save to temporary file.
#' @param verbose verbosity
#' @param return_code return code? Default FALSE, i.e. return the model 
#' filename.
#' 
#' @returns filename
#' 
#' @examples
#' 
#' parameters <- list(CL = 5, V1 = 50, Q = 5, V2 = 100)
#' iiv <- list(CL = 0.27, Q = 0.49, V1 = 0.15, V2 = 1.3)
#' ruv <- list(add = 1.6, prop = 0.15)
#' 
#' model_file <- new_stan_model(
#'   parameters = parameters,
#'   parameter_definitions = list(
#'    "CL" = "CL * (1.0 + 0.0154 * ((CRCL[j] * 16.6667) - 66.0))",
#'    "Q"  = "Q",
#'    "V1" = "V1 * WT[j]",
#'    "V2" = "V2 * WT[j]",
#'    "KA" = "0" 
#'   ),
#'   covariate_definitions = list(
#'     "CRCL" = "real", # CrCl in L/hr
#'     "WT" = "real"    # WT in kg
#'   ),
#'   obs_cmt = 2,
#'   scale = "(V1 * mean(WT))",
#'   verbose = TRUE
#' )
#' 
#' @export
new_stan_model <- function(
  parameters,
  parameter_definitions,
  ode = NULL,
  covariate_definitions = NULL,
  solver = c("pmx_solve_onecpt", "pmx_solve_twocpt", "pmx_solve_rk45", 
             "pmx_solve_adams", "pmx_solve_bdf"),
  template = "template.stan",
  obs_types = NULL,
  obs_cmt = NULL,
  scale = NULL,
  custom_ipred = NULL,
  n_cmt = NULL,
  n_theta = NULL,
  file = NULL,
  verbose = FALSE,
  return_code = FALSE
) {
  
  solver <- match.arg(solver)
  if(!is.null(ode) && ! solver %in% c("pmx_solve_rk45", "pmx_solve_adams", "pmx_solve_bdf")) {
    message("When an ODE block is supplied, a Torsten ODE-solver is required. Switching to `pmx_solve_rk45`.")
    solver <- "pmx_solve_rk45"
  }
  if(is.null(obs_types)) obs_types <- "pk"
  if(is.null(scale) && is.null(custom_ipred)) {
    stop("Argument `scale` or `custom_ipred` is required.")
  }
  if(is.null(obs_cmt) && is.null(custom_ipred)) {
    stop("Argument `obs_cmt` or `custom_ipred` is required.")
  }
  if(!is.null(scale)) {
    if(class(scale) == "character") {
      scale_char <- scale
      scale <- list()
      for(key in obs_types) {
        if(verbose) message(paste0("Assuming `scale` definition (`", scale_char, "`) for observation type `", key, "`"))
        scale[[key]] <- scale_char
      }
    }
  }
  if(!is.null(obs_cmt)) {
    if(class(obs_cmt) %in% c("numeric", "integer")) {
      obs_cmt_num <- obs_cmt
      obs_cmt <- list()
      for(key in obs_types) {
        if(verbose) message(paste0("Assuming `obs_cmt` definition (`", obs_cmt_num, "`) for observation type `", key, "`"))
        obs_cmt[[key]] <- obs_cmt_num
      }
    }
  }

  ## read template file
  template_file <- system.file(paste0("models/", template), package = "PKPDposterior")
  if(!file.exists(template_file)) {
    templates <- dir(system.file("models", package = "PKPDposterior"), pattern = "templ")
    templates <- gsub("\\.stan", "", templates)
    stop("Requested template not found, please specify one of available templates: ", paste0())
  }
  template_code <- readLines(template_file)
  
  if(is.null(n_theta)) {
    n_theta <- length(parameter_definitions)
  }

  def <- parse_model_definitions(
    parameters = parameters,
    parameter_definitions = parameter_definitions,
    ode = ode,
    covariate_definitions = covariate_definitions,
    solver = solver,
    obs_types = obs_types,
    obs_cmt = obs_cmt,
    n_theta = n_theta,
    n_cmt = n_cmt,
    scale = scale,
    custom_ipred
  )

  def[["ode_function"]] <- new_ode_function(
    ode = ode,
    parameters = names(parameter_definitions),
    n_cmt = n_cmt
  )
  
  model_code <- generate_stan_code(
    template_code = template_code,
    definitions = def
  )

  if(verbose) {
    writeLines(model_code)
  }

  ## Save to file and return filename
  if(return_code) {
    return(model_code)
  } 
  if(is.null(file)) {
    file <- tempfile(fileext = ".stan")
  }
  
  writeLines(model_code, con=file)
  file
  
}
