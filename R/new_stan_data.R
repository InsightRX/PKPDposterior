#' Create a data object for use in get_mcmc_posterior()
#' 
#' @param regimen Regimen object (created by [PKPDsim::new_regimen()])
#' @param covariates List of covariate objects created by 
#'   [PKPDsim::new_covariate()]
#' @param data Data frame with columns `t`, `dv`, and `cmt`
#' @param dose_cmt Specify what dose compartment. Observation compartment in 
#' dataset is irrelevant, handled in model.
#' @param parameters list of population parameters, e.g. `list(CL = 5, V = 50)`
#' @param fix optional, vector of parameters to fix. These parameters will be
#' still be sampled from the posterior but from a negligibly narrow distribution
#' around the central value specified in `parameters`.
#' @param iiv list of inter-individual variability for parameters. Should have 
#' exact same list elements as `parameters`, and magnitude supplied on SD scale. 
#' @param ruv magnitude of residual unexplained variability (RUV). Should be a 
#' list specifying proportional and/or additive error magnitude on standard 
#' deviation scale, e.g. `list("prop" = 0.1, "add" = 1)`. If `ltbs` is TRUE, 
#' should specify only an `add` part, which applies an additive error on the 
#' log-scale (which then becomes an approximate proportional error).
#' @param ltbs use log-transform-both-sides approach for observations? Default 
#' is `FALSE`.
#' @param verbose verbosity
#' 
#' @return Named list suitable for passing on to Torsten.
#' 
#' @export
#' 
#' @examples
#' regimen <- PKPDsim::new_regimen(
#'   amt = 1500, 
#'   n = 4, 
#'   times = c(0, 12, 24, 36), 
#'   type = 'infusion'
#' )
#' covariates <- list(
#'   WT = PKPDsim::new_covariate(
#'     value = c(150, 149.5),
#'     times = c(0, 30),
#'     unit = "kg"
#'   ),
#'   CRCL = PKPDsim::new_covariate(
#'     value = c(6.5, 6.7),
#'     times = c(0, 12),
#'     unit = "l/hr"
#'   )
#' )
#' tdm_data <- data.frame(
#'   t = c(1, 2), 
#'   dv = c(900, 800), 
#'   cmt = c(2, 2)
#' )
#' new_stan_data(
#'   regimen, 
#'   covariates, 
#'   tdm_data,
#'   parameters = list(CL = 5, V = 50),
#'   iiv = list(CL = 0.1, V = 0.2),
#'   ruv = list(prop = 0.1, add = 1)
#' )
#' 
new_stan_data <- function(
  regimen, 
  covariates, 
  data,
  parameters,
  fix = NULL,
  iiv,
  ruv,
  dose_cmt = 1,
  ltbs = FALSE,
  verbose = FALSE
) {
  
  stan_data <- PKPDsim_to_stan_data(
    regimen = regimen, 
    covariates = covariates, 
    data = data,
    parameters = parameters,
    fix = fix,
    iiv = iiv,
    ruv = ruv,
    dose_cmt = dose_cmt,
    ltbs = ltbs,
    verbose = verbose
  )
  
  list(
    parameters = parameters,
    fix = fix,
    regimen = regimen,
    covariates = covariates,
    data = data,
    iiv = iiv,
    ruv = ruv,
    stan_data = stan_data
  )
  
}