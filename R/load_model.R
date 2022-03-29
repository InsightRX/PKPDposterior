#' Initialize a model
#' If model is not compiled yet, than compile it. Otherwise used the binary.
#' 
#' @param model model to be compiled, from the models available in the package in `inst/models`
#' @param force force recompile
#' @param verbose show output from Stan / cmdstanr? Defaults to `FALSE`.
#' @param ... passed onto 
#' 
#' @export
load_model <- function(
  model,
  force = FALSE,
  verbose = FALSE,
  ...
) {
  
  messages <- ifelse(verbose, function(x) { x }, suppressMessages)
  
  ## Initialize
  messages(
    cmdstanr::set_cmdstan_path(
      path = file.path(Sys.getenv("STAN_PATH"), "cmdstan")
    )
  )

  ## Compile Stan/Torsten model, or re-use old model
  model_file <- system.file("models", paste0(model, ".stan"), package = "PKPDposterior")
  if(model_file == "") {
    stop("The requested model was not found.")
  }
  if(force) {
    unlink(system.file("models", model, package = "PKPDposterior"))
  }
  messages(
    stan_out <- cmdstanr::cmdstan_model(
      stan_file = model_file
    )
  )
  stan_out
}
  
