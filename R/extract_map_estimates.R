#' Extra MAP estimates from the posterior
#' 
#' Reports the mode for the posterior of each parameter
#' @param post posterior estimate object from `get_mcmc_posterior()`.
#' @param parameters parameters to fetch.
#' @param ... arguments passed on to `get_mode` for calculating the parameter
#'   mode from the sampled values.
#' 
#' @return list of modes of posterior (MAP estimate) for the parameters
#' @keywords internal
extract_map_estimates <- function(
  post,
  parameters = NULL,
  ...
) {

  if(is.null(parameters)) { # attempt to identify parameters
    parameters <- names(post$draws_df)
    idx <- c(
      grepl("^lp__", parameters) |
      grepl("^A\\[", parameters) |
      grepl("^ipred", parameters) |
      grepl("^dv", parameters) |
      grepl("^theta\\[", parameters) |
      grepl("^prior_", parameters) |
      grepl("\\.", parameters) 
    )
    par_fetch <- !idx
  } else {
    par_fetch <- parameters[parameters %in% names(post$draws_df)]
  }
  suppressWarnings(
    par_samples <- post$draws_df[, par_fetch]
  )

  lapply(par_samples, get_mode, ...)
  
}
