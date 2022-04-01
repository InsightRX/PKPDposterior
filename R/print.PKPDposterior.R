#' Print function for PKPDposterior object
#' 
#' @param x posterior object (from get_mcmc_posterior())
#' @param ... arguments passed onto print function
#' 
#' @export
print.PKPDposterior <- function(x, ...) {
  cat("Parameters:\n")
  par_table <- data.frame(
    posterior_mode = round(unlist(x$map), 3)
  )
  post_info <- posterior::summarise_draws(
    posterior::subset_draws(
      x$draws_df, 
      names(x$map)
    )
  )
  par_table <- dplyr::bind_cols(
    par_table,
    post_info[, c("mean", "median", "sd", "q5", "q95", "rhat")] %>%
      round(3)
  )
  print(par_table, ...)
  
  print_cols <- c("time", "dv", "mean", "loc", "pct", "pct5", "pct95")
    cat("\nObserved data, posterior:\n")
  print(x$observed_post[, print_cols], ...)

  cat("\nObserved data, prior:\n")
  print(x$observed_prior[, print_cols], ...)
  
}
