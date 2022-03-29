#' Simulate from posterior or prior draws using PKPDsim::sim()
#' 
#' @param post posterior object from `get_mcmc_posterior()`
#' @param model PKPDsim model function
#' @param map list with parameter name remapping between Stan model and PKPDsim model.
#' @param parameters list of parameter estimates that are used in simulation if 
#' not already present in either the posterior or the prior draws. I.e. 
#' specifiied parameters in this named list will not override parameters with 
#' the same name in the samples.
#' @param prior simulate from prior? Default `FALSE`, i.e. simulate from 
#' posterior
#' @param n number of parameter draws / simulations to do. Defaults to `NULL`,
#' i.e. use all draws.
#' @param summarize should data be summarized to median and confidence interval? Default is `FALSE`.
#' Defaults to `FALSE`, i.e. to return all simulated observations from 
#' parameter draws.
#' @param ci confidence interval to use when `summary=TRUE`. Defaults to 
#' `c(0.05, 0.95)`
#' 
#' @export
sim_from_draws <- function(
  post,
  model,
  map = NULL,
  parameters = NULL,
  prior = FALSE,
  n = NULL,
  summarize = FALSE,
  ci = c(0.05, 0.95),
  ...
) {
  
  par_table <- get_parameter_tables(post, long = FALSE)
  if(!is.null(n)) {
    par_table$posterior <- par_table$posterior %>%
      dplyr::slice(1:n)
    if(prior) {
      par_table$prior <- par_table$prior %>%
        dplyr::slice(1:n)
    }
  }
  
  if(!is.null(map)) {
    par_table$posterior <- remap(par_table$posterior, map, reverse = FALSE)
    if(prior) {
      par_table$prior <- remap(par_table$prior, map, reverse = FALSE)
    }
  }

  if(! all(attr(model, "parameters") %in% c(names(par_table$posterior), names(parameters)))) {
    stop("Not all parameters are available in posterior draws or specified using the `parameters` argument. If parameter names are different in Stan vs PKPDsim model, please use `map` argument to translate parameter names.")
  }

  if(prior) {
    for(key in names(parameters)) {
      if(is.null(par_table$prior[[key]])) {
        message(paste0("Using fixed estimate for ", key))
        par_table$prior[[key]] <- parameters[[key]]
      }
    }
    res <- PKPDsim::sim(
      ode = model,
      parameters_table = as.data.frame(par_table$prior),
      ...
    )
  } else {
    if(!is.null(parameters)) {
      for(key in names(parameters)) {
        if(is.null(par_table$posterior[[key]])) {
          message(paste0("Using fixed estimate for ", key))
          par_table$posterior[[key]] <- parameters[[key]]
        }
      }
    }
    res <- PKPDsim::sim(
      ode = model,
      parameters_table = as.data.frame(par_table$posterior),
      output_include = list(variables = TRUE),
      ...
    )
  }
  
  if(!summarize) {
    res
  } else {
    res %>% 
      dplyr::filter(comp == "obs") %>%
      dplyr::group_by(t) %>%
      dplyr::summarise(
        ymedian = stats::median(y, na.rm = TRUE), 
        ymin = stats::quantile(y, ci[1], na.rm = TRUE),
        ymax = stats::quantile(y, ci[2], na.rm = TRUE)
      )
  }

}