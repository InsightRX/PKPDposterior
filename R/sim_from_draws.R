#' Simulate from posterior or prior draws using PKPDsim::sim()
#' 
#' @param post posterior object from `get_mcmc_posterior()`
#' @param model PKPDsim model function
#' @param map list with parameter name remapping between Stan model and 
#'   PKPDsim model.
#' @param prior simulate from prior (`TRUE`) or posterior (`FALSE`) 
#' @param n number of parameter draws / simulations to do. A value of `NULL` 
#'   indicates all draws should be used.
#' @param variable use a variable from the simulated data, defaults to `y`.
#' @param summarize should data be summarized to median and confidence interval? 
#'   Defaults to `FALSE`, i.e. to return all simulated observations from 
#'   parameter draws.
#' @param ci confidence interval to use when `summary=TRUE`. 
#' @param ... arguments passed on to `PKPDsim::sim`
#' 
#' @export
sim_from_draws <- function(
  post,
  model,
  map = NULL,
  prior = FALSE,
  n = NULL,
  variable = "y",
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
  
  fixed <- post$data$fixed
  pkpdsim_parameter_names <- attr(model, "parameters")
  if(length(grep("kappa_", pkpdsim_parameter_names)) > 0) { ## filter out IOV parameters
    pkpdsim_parameter_names <- pkpdsim_parameter_names[-grep("kappa_", pkpdsim_parameter_names)]
  }

  ## add IOV parameters, if not already specified
  iov <- attr(model, "iov")
  par_type <- ifelse(prior, "prior", "posterior")
  if(!is.null(iov)) {
    for(key in names(iov$cv)) {
      for(i in seq(1, length(iov$bins)-1)) {
        if(is.null(par_table[[par_type]][[paste0("kappa_", key, "_", i)]])) {
          par_table[[par_type]][[paste0("kappa_", key, "_", i)]] <- 0
        }
      }
    }
  }
  ## add fixed parameters
  if(!is.null(fixed)) {
    for(key in fixed) {
      message(paste0("Using fixed estimate for ", key))
      par_table[[par_type]][[key]] <- post$data$parameters[[key]]
    }
  }
  
  res <- PKPDsim::sim(
    ode = model,
    parameters_table = as.data.frame(par_table[[par_type]]),
    output_include = list(variables = TRUE),
    ...
  )
  
  res$y <- res[[variable]]
  
  if(!summarize) {
    res
  } else {
    res %>% 
      dplyr::filter(.data$comp == "obs") %>%
      dplyr::group_by(.data$t) %>%
      dplyr::summarise(
        ymedian = stats::median(.data$y, na.rm = TRUE), 
        ymin = stats::quantile(.data$y, ci[1], na.rm = TRUE),
        ymax = stats::quantile(.data$y, ci[2], na.rm = TRUE)
      )
  }

}
