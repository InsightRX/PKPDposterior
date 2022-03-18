#' Simulate from posterior or prior draws using PKPDsim::sim()
#' 
#' @param post posterior object from `get_mcmc_posterior()`
#' @param model PKPDsim model function
#' @param prior simulate from prior? Default `FALSE`, i.e. simulate from 
#' posterior
#' @param n number of parameter draws / simulations to do. Defaults to `NULL`,
#' i.e. use all draws.
#' 
#' @export
sim_from_draws <- function(
  post,
  model,
  prior = FALSE,
  n = NULL,
  ...
) {
  
  par_table <- get_parameter_tables(post, long = FALSE)
  if(!is.null(n)) {
    par_table$posterior <- par_table$posterior %>%
      dplyr::slice(1:n)
    par_table$prior <- par_table$prior %>%
      dplyr::slice(1:n)
  }

  if(prior) {
    par_table$prior <- par_table$prior %>% 
      dplyr::mutate(V = V1, TH_CRCL = 0.0154, TDM_INIT = 0)
    res <- PKPDsim::sim(
      ode = model,
      parameters_table = as.data.frame(par_table$prior),
      ...
    )
  } else {
    par_table$posterior <- par_table$posterior %>% 
      dplyr::mutate(V = V1, TH_CRCL = 0.0154, TDM_INIT = 0)
    res <- PKPDsim::sim(
      ode = model,
      parameters_table = as.data.frame(par_table$posterior),
      ...
    )
  }
  
  res
  
}