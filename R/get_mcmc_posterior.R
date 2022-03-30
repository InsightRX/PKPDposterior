#' Draw samples from the posterior given a Stan/Torsten model
#' and patient data
#' 
#' @param mod compiled Stan model
#' @param init initial parameter set for sampler
#' @param data dataset (see [prepare_data()])
#' @param seed seed for sampling
#' @param chains number of MCMC chains to simulate, passed on to Stan model
#' @param refresh show output from sampler. Default is 0, meaning no output 
#'   from sampler is shown.
#' @param output_dir output directory
#' @param verbose verbosity
#' @param skip_processing will return the raw output from the sampler, without
#'   further processing. This can sometimes be useful when errors occur in the
#'   processing stage.
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
  chains = 1,
  output_dir = tempdir(),
  verbose = TRUE,
  refresh = 0,
  skip_processing = FALSE,
  ...
) {
  
  ## Check model OK?
  if(! "CmdStanModel" %in% class(mod)) {
    stop(
      "Supplied model is not a valid cmdstanr model. ", 
      "Please use `load_model()` to load/compile models."
    )
  }
  if(verbose) (
    res <- mod$sample(
      data = data,
      init = function() { init },
      seed = seed,
      refresh = refresh,
      output_dir = output_dir,
      chains = chains,
      ...
    )
  ) else {
    suppressMessages(
      res <- mod$sample(
        data = data,
        init = function() { init },
        seed = seed,
        refresh = refresh,
        output_dir = output_dir,
        chains = chains,
        ...
      )
    )
  }

  ## Post-processing
  out <- list(
    raw = res,
    draws_df =  posterior::as_draws_df(res$draws()),
    settings = list(
      init = init,
      seed = seed
    )
  )
  
  if(skip_processing) {
    return(out)
  }
  
  out$map <- extract_map_estimates(out)
    
  ## Model diagnostics
  out$observed_post <- extract_from_draws(
    out, 
    data,
    filter = "ipred_obs",
    verbose
  )

  ## Sampler diagnostics
  out$sampler_diagnostics <- posterior::summarise_draws(
    res$sampler_diagnostics()
  )

  ## Finalize & return
  class(out) <- c("PKPDposterior", "list")
  out
}
