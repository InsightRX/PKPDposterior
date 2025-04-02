#' Plot posterior distributions of model parameters
#' If the model also supplies prior parameters, specified e.g. as `prior_CL` 
#' for the `CL` parameter, it will plot them as well.
#' 
#' @param post posterior object from `get_mcmc_posterior()`
#' @param params vector of parameter names to extract from table with 
#' posterior draws. If `NULL`, will take default parameters from `init` object.
#' 
#' @return a ggplot2 object with posterior distributions (histograms) faceted 
#' by parameter.
#'  
#' @export
#' 
plot_params <- function(
  post,
  params = NULL
) {
  
  par_table_long <- get_parameter_tables(
    post,
    params = params,
    long = TRUE
  )
  pars <- par_table_long$posterior %>%
    dplyr::mutate(sample = "posterior")
  if(!is.null(par_table_long$prior)) {
    pars <- dplyr::bind_rows(
      pars,
      par_table_long$prior %>%
        dplyr::mutate(sample = "prior")
    )
  }

  suppressMessages({ # histogram function is too chatty
    p <- ggplot2::ggplot() + 
      ggplot2::geom_histogram(
        data = pars,
        ggplot2::aes(
          x = .data$value, 
          fill = .data$sample,
          colour = .data$sample
        ), 
        alpha = 0.4
      ) + 
      ggplot2::facet_wrap(~name, scale = "free")
  })

  p
}
