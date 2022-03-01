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

  if(is.null(parameters)) { # attempt to identify parameters
    parameters <- names(post$draws_df)
    idx <- c(
      stringr::str_detect(parameters, "^lp__") |
      stringr::str_detect(parameters, "^x\\[") |
      stringr::str_detect(parameters, "^cHat") |
      stringr::str_detect(parameters, "^cObs") |
      stringr::str_detect(parameters, "^theta\\[") |
      stringr::str_detect(parameters, "^prior_") |
      stringr::str_detect(parameters, "\\.") 
    )
    par_fetch <- !idx
  } else {
    par_fetch <- parameters[parameters %in% names(post$draws_df)]
  }
  suppressWarnings(
    par_samples <- post$draws_df[,!idx]
  )
  lapply(par_samples, get_mode)
  
}
