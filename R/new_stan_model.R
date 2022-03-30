#' Create a new stan model using a template
#' 
#' @param parameters list of parameter values, e.g. `list(CL = 5, V = 50)`.
#' @param parameter_definitions list of parameter definitions. See example 
#' below.
#' @param covariate_definitions list of covariate definitions. See example
#' below.
#' @param solver Torsten solver to use. Starts with `pmx_solve_`.
#' @param template internal Stan template file to use. Defaults to 
#' `template.stan`
#' @param obs_types vector of observed data types. If `NULL`, will use `"pk"`.
#' @param scale scale definition. Can be either parameter expression or number.
#' @param n_cmt number of compartments. If analytic equation, will use number 
#' of compartments for solver.
#' @param n_theta number of thetas. If analytic equation, will use number of 
#' thetas for solver.
#' @param file file to save model to. If NULL, will save to temporary file.
#' @param verbose verbosity
#' 
#' @returns filename
#' 
#' @example
#' 
#' parameters <- list(CL = 5, V1 = 50, Q = 5, V2 = 100)
#' iiv <- list(CL = 0.27, Q = 0.49, V1 = 0.15, V2 = 1.3)
#' ruv <- list(add = 1.6, prop = 0.15)
#' 
#' model_file <- new_stan_model(
#' parameters = parameters,
#' parameter_definitions = list(
#'  "CL" = "CL * (1.0 + 0.0154 * ((CRCL[j] * 16.6667) - 66.0))",
#'  "Q"  = "Q",
#'  "V1" = "V1 * WT[j]",
#'  "V2" = "V2 * WT[j]",
#'  "KA" = "0" 
#' ),
#' covariate_definitions = list(
#'   "CRCL" = "real", # CrCl in L/hr
#'   "WT" = "real"    # WT in kg
#' ),
#' scale = "(V1 * mean(WT))"
#' )
#' 
#' @export
new_stan_model <- function(
  parameters,
  parameter_definitions,
  covariate_definitions = NULL,
  solver = c("pmx_solve_onecpt", "pmx_solve_twocpt", "pmx_solve_rk45"),
  template = "template.stan",
  obs_types = NULL,
  scale,
  n_cmt = 5,
  n_theta = NULL,
  file = NULL,
  verbose = FALSE
) {
  
  solver <- match.arg(solver)
  if(is.null(obs_types)) obs_types <- "pk"
  if(class(scale) == "character") {
    scale_char <- scale
    scale <- list()
    for(key in obs_types) {
      message(paste0("Assuming `scale` definition (`", scale_char, "`) for observation type `", key, "`"))
      scale[[key]] <- scale_char
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
  
  ## Determine n_theta

  def <- parse_model_definitions(
    parameters = parameters,
    parameter_definitions = parameter_definitions,
    covariate_definitions = covariate_definitions,
    solver = solver,
    obs_types = obs_types,
    n_theta = n_theta,
    n_cmt = n_cmt,
    scale = scale
  )
  
  model_code <- generate_stan_model(
    template_code = template_code,
    definitions = def
  )
  
  if(verbose) {
    writeLines(model_code)
  }

  ## Save to file and return filename
  if(is.null(file)) {
    file <- tempfile(fileext = ".stan")
  }
  writeLines(model_code, con=file)

  file
  
}
