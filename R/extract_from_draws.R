#' Extract info for observed data from posterior draws.
#' Includes e.g. the quantile of the observed data in the posterior draws
#' and the mean, 95% CI, and standard deviations of the posterior.
#' 
#' @param post posterior object from CmdStanR
#' @param verbose verbosity
#' 
extract_from_draws <- function(
  post,
  verbose = TRUE
) {
  
  obs_types <- gsub("dv_", "", names(post$data$stan_data)[grep("^dv_", names(post$data$stan_data))])
  obs_data_all <- purrr::map_dfr(
    obs_types, 
    calc_stats_for_draws_predictions, 
    post = post
  )

  ## warnings of poor fit
  if(any(obs_data_all$pct > 0.99) || any(obs_data_all$pct < 0.01)) {
    if(verbose) message(
      "One or more observed data points were at the edges of the posterior ",
      "(1st percentile)"
    )
  }
  
  obs_data_all
  
}

#' Calculates statistics on draws from posterior for predictions (e.g. 
#' concentrations)
#' 
#' @param post posterior object
#' @param obs_type observation type label, default is `pk`
#' 
#' @returns data frame with statistics for each observation timepoint, such as
#' mean, median, sd, confidence intervals etc.
#' 
calc_stats_for_draws_predictions <- function(
  post,
  obs_type = "pk"
) {
  obs_post <- as.data.frame(post$draws_df) # needs to be a standard data.frame
  obs_post <- obs_post[, grepl(paste0("ipred_obs_", obs_type), names(obs_post))]
  post_info <- obs_post %>%
    posterior::summarise_draws(c("mean", "median", "sd"))
  obs_data <- data.frame(
    time = post$data$stan_data$time,
    dv = post$data$stan_data$DV,
    evid = post$data$stan_data$evid
  ) %>%
    dplyr::slice(post$data$stan_data[[paste0("i_obs_", obs_type)]]) %>%
    dplyr::mutate(type = !!obs_type) %>%
    dplyr::bind_cols(post_info)
  for(i in 1:nrow(obs_data)) {
    obs_data$pct[i] <- get_quantile(
      obs_data$dv[i], 
      as.numeric(obs_post[ ,i])
    )
  }
  
  ## add visualization of percentiles for the print function
  obs_data$loc <- plot_percentile(obs_data$pct)
  
  ## add quantiles
  obs_data$pct2.5 <- as.numeric(apply(obs_post, 2, "quantile", 0.025))
  obs_data$pct5 <- as.numeric(apply(obs_post, 2, "quantile", 0.05))
  obs_data$pct10 <- as.numeric(apply(obs_post, 2, "quantile", 0.1))
  obs_data$pct25 <- as.numeric(apply(obs_post, 2, "quantile", 0.25))
  obs_data$pct75 <- as.numeric(apply(obs_post, 2, "quantile", 0.75))
  obs_data$pct90 <- as.numeric(apply(obs_post, 2, "quantile", 0.90))
  obs_data$pct95 <- as.numeric(apply(obs_post, 2, "quantile", 0.95))
  obs_data$pct97.5 <- as.numeric(apply(obs_post, 2, "quantile", 0.975))
  
  obs_data
}