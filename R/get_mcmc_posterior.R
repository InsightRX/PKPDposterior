#' Draw samples from the posterior given a Stan/Torsten model
#' and patient data
#' 
#' @param mod compiled Stan model
#' @param data dataset (see [new_stan_data()])
#' @param seed seed for sampling
#' @param chains number of MCMC chains to simulate, passed on to Stan model
#' @param output_dir output directory
#' @param method either `hmc` (Hamilton Monte Carlo, using NUTS by default) or 
#' `vi` (variational inference).
#' @param verbose verbosity
#' @param skip_processing will return the raw output from the sampler, without
#'   further processing. This can sometimes be useful when errors occur in the
#'   processing stage.
#' @param iter_warmup number of warmup (burn-in) iterations
#' @param iter_sampling number of regular sampling iterations
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
#' @seealso [new_stan_data()]
get_mcmc_posterior <- function(
  mod,
  data,
  seed = 12345,
  chains = 1,
  iter_warmup = 500,
  iter_sampling = 500,
  output_dir = tempdir(),
  verbose = TRUE,
  method = c("hmc", "vi"),
  skip_processing = FALSE,
  ...
) {
  
  method <- match.arg(method)
  sample_func <- list(
    "hmc" = "sample",
    "vi" = "variational" 
  )
    
  ## Check model OK?
  if(! "CmdStanModel" %in% class(mod)) {
    stop(
      "Supplied model is not a valid cmdstanr model. ", 
      "Please use `load_model()` to load/compile models."
    )
  }
  
  ## Check that input data is OK for model
  check_stan_input_data(code = mod$code(), data = data)
 
  ## get parameter initial values from data object
  init <- data$stan_data[names(data$stan_data)[grep("theta_", names(data$stan_data))]]
  names(init) <- gsub("theta_", "", names(init))
  
  refresh <- ifelse(verbose, 50, 0)
  
  if(method == "vi") {
    run_cmdstanr <- function() {
      res <- mod$variational(
        data = data$stan_data,
        seed = seed,
        refresh = refresh,
        output_dir = output_dir,
        ...
      )
    }
  } else {
    run_cmdstanr <- function() {
      mod$sample(
        data = data$stan_data,
        init = function() { init },
        seed = seed,
        refresh = refresh,
        output_dir = output_dir,
        chains = chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        ...
      )
    }
  }
  
  ## Run
  res <- withCallingHandlers({
      run_cmdstanr()
    },
    message = function(c) if (!verbose) tryInvokeRestart("muffleMessage")
  )

  ## Post-processing
  out <- list(
    raw = res,
    draws_df =  posterior::as_draws_df(res$draws()),
    settings = list(
      init = init,
      seed = seed
    ),
    data = data
  )
  
  if(skip_processing) {
    return(out)
  }
  
  out$map <- extract_map_estimates(out)
    
  ## Model diagnostics
  out$observed_post <- extract_from_draws(out, verbose)

  ## Sampler diagnostics
  if(!is.null(res$sampler_diagnostics)) { # not available for VI method
    out$sampler_diagnostics <- posterior::summarise_draws(
      res$sampler_diagnostics()
    )
  }

  ## Finalize & return
  class(out) <- c("PKPDposterior", "list")
  out
}
