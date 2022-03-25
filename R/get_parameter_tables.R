#' Get parameter tables from posterior and prior
#' 
#' @param post posterior object from `get_mcmc_posterior`
#' @param long return tables in long format? Default is `FALSE`
#' 
#' @export
get_parameter_tables <- function(
  post,
  long = FALSE
) {
  if(is.null(post$settings) || is.null(post$draws_df)) {
    stop("Provided posterior object does not contain expected info.")
  }
  params <- names(post$settings$init)
  prior_params <- paste0("prior_", params)
  suppressWarnings({
    par_table <- post$draws_df %>% 
      dplyr::select(!!params)
    par_table_prior <- NULL
    if(all(prior_params %in% names(post$draws_df))) {
      par_table_prior <- post$draws_df %>% 
        dplyr::select(!!prior_params)
      names(par_table_prior) <- gsub("prior_", "", names(par_table_prior))
    }
  })
  
  if(!long) {
    tables <- list(
      posterior = par_table,
      prior = par_table_prior
    )
  } else {
    tables <- list(
      posterior = par_table %>% 
        tidyr::pivot_longer(cols = c(!!params))
    )
    if(!is.null(par_table_prior)) {
      tables$prior = par_table_prior %>% 
        tidyr::pivot_longer(cols = c(!!params))
    }
  }
  
  tables
  
}