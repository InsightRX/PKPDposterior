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
  
  par_table_long <- get_parameter_tables(post, long = TRUE)

  suppressMessages({ # histogram function is too chatty
    p <- ggplot2::ggplot()
    if(!is.null(par_table_long$prior)) {
      p <- p + ggplot2::geom_histogram(
        data = par_table_long$prior, 
        ggplot2::aes(x = value), 
        alpha = 0.2, 
        fill = "black"
      )
    }
    p <- p +
      ggplot2::geom_histogram(
        data = par_table_long$posterior,
        ggplot2::aes(x = value), alpha = 0.4, fill = "#3366ee"
      ) + 
      ggplot2::facet_wrap(~name, scale = "free") +
      irxreports::theme_irx_minimal()
  })

  p
}
