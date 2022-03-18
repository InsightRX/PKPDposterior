#' Draw samples from the posterior given a Stan/Torsten model
#' and patient data
#' 
#' @param mod compiled Stan model
#' @param init initial parameter set for sampler
#' @param data dataset (see [prepare_data()])
#' @param seed seed for sampling
#' @param refresh show output from sampler. Default is 0, meaning no output 
#'   from sampler is shown.
#' @param output_dir output directory
#' @param verbose verbosity
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
  output_dir = tempdir(),
  verbose = TRUE,
  refresh = 0,
  ...
) {
  
  ## Check model OK?
  if(! "CmdStanModel" %in% class(mod)) {
    stop("Supplied model is not a valid cmdstanr model. Please use `load_model()` to load/compile models.")
  }

  suppressMessages( # don't output divergencies etc, we can diagnose that at the end.
    res <- mod$sample(
      data = data,
      init = function() { init },
      seed = seed,
      refresh = refresh,
      output_dir = output_dir,
      chains = 1,
      ...
    )
  )

  ## Post-processing
  out <- list(
    raw = res,
    draws_df =  posterior::as_draws_df(res$draws()),
    settings = list(
      init = init,
      seed = seed
    )
  )
  out$map <- extract_map_estimates(out)
  
  ## Model diagnostics
  out$observed_post <- extract_from_draws(
    out, 
    data,
    filter = "cHatObs",
    verbose
  )
  out$observed_prior <- extract_from_draws(
    out, 
    data,
    filter = "cObsPred",
    verbose = FALSE
  )
  
  ## Sampler diagnostics
  out$sampler_diagnostics <- posterior::summarise_draws(
    res$sampler_diagnostics()
  )

  ## Finalize & return
  class(out) <- c("PKPDposterior", "list")
  out
}
