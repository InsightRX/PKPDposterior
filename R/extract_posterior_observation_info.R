#' Extract info for observed data from posterior draws.
#' Includes e.g. the quantile of the observed data in the posterior draws
#' and the meanm, 95% CI, and standard deviations of the posterior.
#' 
#' @param post posterior object from CmdStanR
#' @param data data.frame with observed data
#' @param verbose verbosity
#' 
extract_posterior_observation_info <- function(
  post, 
  data,
  verbose = TRUE
) {
  post_info <- subset_draws(post$draws_df, "cHatObs") %>%
    summarise_draws(c("mean", "median", "sd"))
  obs_data <- data.frame(
    time = data$time,
    dv = data$DV,
    evid = data$evid
  ) %>%
    dplyr::filter(evid == 0) %>%
    dplyr::bind_cols(post_info[,-1])
  
  ## add percentile of observed in 
  for(i in 1:nrow(obs_data)) {
    obs_data$pct[i] <- get_quantile(
      obs_data$dv[i], 
      obs_post[,i]
    )
  }
  
  ## add quantiles
  obs_data$pct2.5 <- as.numeric(apply(obs_post, 2, quantile, 0.025))
  obs_data$pct5 <- as.numeric(apply(obs_post, 2, quantile, 0.05))
  obs_data$pct10 <- as.numeric(apply(obs_post, 2, quantile, 0.1))
  obs_data$pct25 <- as.numeric(apply(obs_post, 2, quantile, 0.25))
  obs_data$pct75 <- as.numeric(apply(obs_post, 2, quantile, 0.75))
  obs_data$pct90 <- as.numeric(apply(obs_post, 2, quantile, 0.90))
  obs_data$pct95 <- as.numeric(apply(obs_post, 2, quantile, 0.95))
  obs_data$pct97.5 <- as.numeric(apply(obs_post, 2, quantile, 0.975))
  
  ## warnings of poor fit
  if(any(obs_data$pct > 0.99) || any(obs_data$pct < 0.01)) {
    if(verbose) message("One or more observed data points were at the edges of the posterior (1st percentile)")
  }
  
  obs_data
  
}
