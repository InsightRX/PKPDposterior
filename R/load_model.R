#' Initialize a model
#' If model is not compiled yet, than compile it. Otherwise used the binary.
#' 
#' @param model_file model file to be compiled, either full path to a model,
#' or else will look in `models` folder inside the installed package folder.
#' @param force force recompile
#' @param verbose show output from Stan / cmdstanr? Defaults to `FALSE`.
#' @param ... passed onto 
#' 
#' @export
load_model <- function(
  model_file,
  force = FALSE,
  verbose = FALSE,
  ...
) {
  
  messages <- ifelse(verbose, function(x) { x }, suppressMessages)

  ## Compile Stan/Torsten model, or re-use old model
  if(!file.exists(model_file)) {
    model_file <- system.file(
      "models",
      paste0(model_file, ".stan"),
      package = "PKPDposterior"
    )
    if(model_file == "") {
      stop("The requested model was not found.")
    } else if(verbose) {
      message("Model found in internal library.")
    }
  }

  messages(
    stan_out <- cmdstanr::cmdstan_model(
      stan_file = model_file
    )
  )
  stan_out
}
  
