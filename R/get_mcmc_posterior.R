#' Draw samples from the posterior given a Stan/Torsten model
#' and patient data
#' 
#' @param mod compiled Stan model
#' @param regimen drug regimen data created using `PKPDsim::new_regimen()`
#' @param covariates patient covariate data created using `PKPDsim::covariates()`
#' @param data data.frame with observed data.
#' @param nm_data NONMEM-style dataset. If `nm_data` supplied, then `regimen`, `covariates`, and `data` are not required.
#' @param seed seed for sampling
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
  prior,
  regimen = NULL,
  covariates = NULL,
  data = NULL,
  nm_data = NULL, 
  seed = 12345,
  output_dir = getwd(),
  ...
) {
  
  ## Check model OK?
  if(! "CmdStanModel" %in% class(mod)) {
    stop("Supplied model is not a valid cmdstanr model. Please use `load_model()` to load/compile models.")
  }
  
  ## convert regimen, covariates, and data to NONMEM-like data
  if(!is.null(nm_data)) { # NONMEM-style data, can be passed directly to cmdstanr
    stan_data <- nm_data  
    ## TODO: check input data
  } else {
    if(is.null(regimen) || is.null(covariates) || is.null(data)) {
      stop("Arguments `regimen`, `covariates`, and `data` are required. Alternatively, a NONMEM-style dataset can be provided using the `nm_data` argument")
    }
    stan_data <- prepare_data(
      regimen,
      covariates,
      data
    )
  }

  res <- mod$sample(
    data = stan_data,
    init = function() { prior },
    seed = seed,
    output_dir = output_dir,
    ...
  )
  
  ## Post-processing
  out <- list(
    raw = res,
    draws_df =  posterior::as_draws_df(res$draws())
  )

  ## return
  out
}