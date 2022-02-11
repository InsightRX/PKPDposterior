#' Draw samples from the posterior given a Stan/Torsten model
#' and patient data
#' 
#' @param mod compiled Stan model
#' @param regimen drug regimen data created using `PKPDsim::new_regimen()`
#' @param covariates patient covariate data created using `PKPDsim::covariates()`
#' @param output_dir 
#' @param ... arguments passed to `mod$sample()`
#' 
#' Common arguments to cmdstanr include:
#' - chains
#' - parallel_chains
#' - iter_warmup
#' - iter_sampling
#' - thim
#' - adapt_delta
#' 
#' @export
get_mcmc_posterior <- function(
  mod,
  regimen,
  covariates,
  data,
  output_dir = getwd(),
  ...
) {
  
  ## Check model OK?
  if(! "CmdStanModel" %in% class(mod)) {
    stop("Supplied model is not a valid cmdstanr model. Please use `load_model()` to load/compile models.")
  }
  
  ## Check data OK?
  ## TODO
  
  ## convert regimen, covariates, and data to NONMEM-like data
  data <- prepare_data(
    regimen,
    covariates,
    data
  )
  
  res <- mod$sample(
    data = data,
    output_dir = output_dir,
    ...
  )
  
  ## Post-processing & return
  res
  
}