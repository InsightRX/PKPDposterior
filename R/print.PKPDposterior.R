#' Print function for PKPDposterior object
#' 
#' @param x posterior object (from get_mcmc_posterior())
#' @param ... arguments passed onto print function
#' 
#' @export
#' @keywords internal
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
  
  cat("\nObserved data, posterior:\n")
  print_cols <- c("time", "dv", "mean", "loc", "pct", "pct5", "pct95")
  if("type" %in% names(x$observed_post)) print_cols <- c("type", print_cols)
  print(x$observed_post[, print_cols], ...)

}
