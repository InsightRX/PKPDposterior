#' Extra MAP estimates from the posterior
#' 
#' Reports the mode for the posterior of each parameter
#' @param post posterior estimate object from `get_mcmc_posterior()`.
#' 
#' @export
extract_map_estimates <- function(
  post
) {
  
  par_fetch <- parameters[parameters %in% names(post$draws_df)]
  par_samples <- post$draws_df[, par_fetch]
  
  lapply(par_samples, get_mode)
  
}
