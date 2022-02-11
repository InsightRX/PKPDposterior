#' Draw samples from the posterior given a Stan/Torsten model
#' and patient data
#' 
#' @param mod compiled Stan model
#' @param init initial parameter set for sampler
#' @param data dataset (see [prepare_data()])
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
#' @seealso [prepare_data()]
get_mcmc_posterior <- function(
  mod,
  init,
  data,
  seed = 12345,
  output_dir = getwd(),
  ...
) {
  
  ## Check model OK?
  if(! "CmdStanModel" %in% class(mod)) {
    stop("Supplied model is not a valid cmdstanr model. Please use `load_model()` to load/compile models.")
  }

  res <- mod$sample(
    data = data,
    init = function() { init },
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
