#' Create a new stan model using a template
#' 
#' @export
new_stan_model <- function(
  parameters,
  parameter_definitions,
  covariate_definitions = NULL,
  solver = c("pmx_solve_twocpt", "pmx_solve_rk45"),
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
    file <- tempfile()
  }
  writeLines(model_code, con=file)

  file
  
}
