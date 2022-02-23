#' Extra MAP estimates from the posterior
#' 
#' Reports the mode for the posterior of each parameter
#' @param post posterior estimate object from `get_mcmc_posterior()`.
#' @param parameters parameters to fetch.
#' 
#' @returns list of modes of posterior (MAP estimate) for the parameters
#' 
#' @export
extract_map_estimates <- function(
  post,
  parameters = NULL
) {
  
  if(is.null(parameters)) {
    par_samples <- post$draws_df
  } else {
    par_fetch <- parameters[parameters %in% names(post$draws_df)]
    par_samples <- post$draws_df[, par_fetch]
  }

  lapply(par_samples, get_mode)
  
}
