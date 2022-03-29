#' Create a new stan model using a template
#' 
#' @export
new_stan_model <- function(
  parameters,
  parameter_definitions,
  covariate_definitions = NULL,
  solver = c("pmx_solve_twocpt", "pmx_solve_rk45"),
  template = "template.stan"
) {
  
  solver <- match.arg(solver)
  
  ## read template file
  template_file <- system.file(paste0("models/", template), package = "PKPDposterior")
  if(!file.exists(template_file)) {
    templates <- dir(system.file("models", package = "PKPDposterior"), pattern = "templ")
    templates <- gsub("\\.stan", "", templates)
    stop("Requested template not found, please specify one of available templates: ", paste0())
  }
  template_code <- readLines(template_file)

  def <- parse_model_definitions(
    parameters = parameters,
    parameter_definitions = parameter_definitions,
    covariate_definitions = covariate_definitions,
    solver = solver
  )
  
  model_code <- generate_stan_model(
    template_code = template_code,
    definitions = def
  )
  
  writeLines(model_code)

}
