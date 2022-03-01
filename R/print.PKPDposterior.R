#' Print function for PKPDposterior object
#' 
#' @param post posterior object (from get_mcmc_posterior())
#' 
#' @export
print.PKPDposterior <- function(post) {
  cat("Parameters:\n")
  par_table <- data.frame(
    posterior_mode = round(unlist(post$map), 3)
  )
  post_info <- summarise_draws(subset_draws(post$draws_df, names(post$map)))
  par_table <- bind_cols(
    par_table,
    post_info[, c("mean", "median", "sd", "q5", "q95", "rhat")] %>%
      round(3)
  )
  print(par_table)
  
  cat("\nObserved data, posterior:\n")
  print(post$observed_post[, c("time", "dv", "mean", "loc", "pct", "pct5", "pct95")])

  cat("\nObserved data, prior:\n")
  print(post$observed_prior[, c("time", "dv", "mean", "loc", "pct", "pct5", "pct95")])
  
}
