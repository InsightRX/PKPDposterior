#' Plot predictions from posterior (or prior) simulations
#' Requires as input a summarized simulation data.frame, e.g. obtained using
#' `sim_from_draws(post, summarize=TRUE)`
#' 
#' @param sim_data data.frame with simulated predictions
#' @param obs observed data. A data.frame with columns `t` and `dv`, e.g.
#' `tdm_data <- data.frame(t = c(2.5, 11.5),  dv = c(40, 14))`
#' 
#' @export
plot_predictions <- function(
  sim_data,
  obs = NULL
) {
  if(!(all(c("t", "ymedian", "ymin", "ymax") %in% names(sim_data)))) {
    stop("Requires a simulation dataset created using sim_from_draws()")
  }
  p <- sim_data %>%
    ggplot(aes(x = .data$t, y = .data$ymedian)) +
    geom_ribbon(aes(ymin = .data$ymin, ymax = .data$ymax), fill = "#cfcfcf") +
    geom_line(size = 1)
  if(!is.null(obs)) {
    p <- p +
      geom_point(
        data = obs, 
        mapping = aes(x = .data$t, y = .data$dv, group = NULL), 
        colour = "black", 
        size = 2.5
      ) +
      geom_point(
        data = obs, 
        mapping = aes(x = .data$t, y = .data$dv, group = NULL), 
        colour = "white", 
        size = 1.5
      )
  }
  p
}
