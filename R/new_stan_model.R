#' Create a new stan model using a template
#' 
#' @param parameters list of parameter values, e.g. `list(CL = 5, V = 50)`.
#' @param parameter_definitions list of parameter definitions. See example 
#' below.
#' @param variable_definitions for models using the analytical solvers, extra 
#' variable definitions evaluated at each step through the event table. 
#' See busulfan example in vignettes for example.
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
#' @param verbose verbosity
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
#'   solver = "pmx_solve_twocpt",
#'   obs_cmt = 2,
#'   scale = "(V1 * mean(WT))",
#'   verbose = TRUE
#' )
#' 
#' @export
new_stan_model <- function(
  parameters,
  parameter_definitions,
  fixed = NULL,
  variable_definitions = NULL,
  ode = NULL,
  covariate_definitions = NULL,
  solver = c("pmx_solve_onecpt", "pmx_solve_twocpt", "pmx_solve_rk45", "pmx_solve_adams", "pmx_solve_bdf"),
  template = "template.stan",
  obs_types = NULL,
  obs_cmt = NULL,
  scale = NULL,
  custom_ipred = NULL,
  verbose = FALSE
) {
  solver <- match.arg(solver)
  if(!is.null(ode) && ! solver %in% c("pmx_solve_rk45", "pmx_solve_adams", "pmx_solve_bdf")) {
    warning(
      "When an ODE block is supplied, a Torsten ODE-solver is required. ",
      "Switching to `pmx_solve_rk45`."
    )
    solver <- "pmx_solve_rk45"
  }
  if(is.null(obs_types)) obs_types <- "pk"
  if(is.null(scale) && is.null(custom_ipred)) {
    stop("Argument `scale` or `custom_ipred` is required.")
  }
  if(is.null(obs_cmt) && is.null(custom_ipred)) {
    if(solver %in% c("pmx_solve_onecpt", "pmx_solve_twocpt")) {
      obs_cmt <- 2 # central compartment is 2 for both these solvers
    } else {
      stop("Argument `obs_cmt` or `custom_ipred` is required.")
    }
  }
  if(!is.null(scale) && inherits(scale, "character")) {
    if (verbose) message("assuming scale definition ", scale, " for all observation types")
    scale <- as.list(setNames(rep(scale, length(obs_types)), obs_types))
  }
  if (isTRUE(mode(obs_cmt) == "numeric")) {
    if(verbose) message("assuming observation type ", obs_cmt, " for all observation types")
    obs_cmt <- as.list(setNames(rep(obs_cmt, length(obs_types)), obs_types))
  }

  ## read template file
  template_file <- system.file(paste0("models/", template), package = "PKPDposterior")
  if(!file.exists(template_file)) {
    templates <- dir(system.file("models", package = "PKPDposterior"), pattern = "templ")
    templates <- gsub("\\.stan", "", templates)
    stop(
      "Requested template not found, please specify one of available templates: ", 
      paste0(templates, collapse = ", ")
      )
  }
  template_code <- readLines(template_file)

  def <- parse_model_definitions(
    parameters = parameters,
    parameter_definitions = parameter_definitions,
    fixed = fixed,
    variable_definitions = variable_definitions,
    ode = ode,
    covariate_definitions = covariate_definitions,
    solver = solver,
    obs_types = obs_types,
    obs_cmt = obs_cmt,
    scale = scale,
    custom_ipred
  )

  model_code <- generate_stan_code(
    template_code = template_code,
    definitions = def
  )

  model_code

}
