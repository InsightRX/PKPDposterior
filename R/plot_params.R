#' Plot posterior distributions of model parameters
#' If the model also supplies prior parameters, specified e.g. as `prior_CL` 
#' for the `CL` parameter, it will plot them as well.
#' 
#' @param post posterior object from `get_mcmc_posterior()`
#' 
#' @return a ggplot2 object with posterior distributions (histograms) faceted 
#' by parameter.
#'  
#' @export
plot_params <- function(post) {
  params <- names(post$settings$init)
  prior_params <- paste0("prior_", params)
  suppressWarnings({
    par_table <- post$draws_df %>% 
      dplyr::select(!!params)
    par_table_long <- par_table %>% 
      tidyr::pivot_longer(cols = c(!!params))
    par_table_prior <- NULL
    if(all(prior_params %in% names(post$draws_df))) {
      par_table_prior <- post$draws_df %>% 
        dplyr::select(!!prior_params)
      names(par_table_prior) <- gsub("prior_", "", names(par_table_prior))
      par_table_prior_long <- par_table_prior %>% 
        tidyr::pivot_longer(cols = c(!!params))
    }
  })
  suppressMessages({ # histogram function is too chatty
    p <- ggplot2::ggplot(par_table_long)
    if(!is.null(par_table_prior_long)) {
      p <- p + ggplot2::geom_histogram(
        data = par_table_prior_long, 
        ggplot2::aes(x = value), 
        alpha = 0.2, 
        fill = "black"
      )
    }
    p <- p +
      ggplot2::geom_histogram(ggplot2::aes(x = value), alpha = 0.4, fill = "#3366ee") + 
      ggplot2::facet_wrap(~name, scale = "free") +
      irxreports::theme_irx_minimal()
  })

  p
}
