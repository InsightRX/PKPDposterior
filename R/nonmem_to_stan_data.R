#' Reshape NONMEM data into PKPDposterior-ready or Stan-ready data
#'
#' @inherit new_stan_data return
#' @inheritParams new_stan_data
#' @param nm NONMEM data frame for a single patient, with columns for TIME,
#'   EVID, DV, and AMT. Observations in DV are assumed to be taken from
#'   compartment 2 unless otherwise specified in an (optional) `"CMT"`
#'   column.
#' @param covariate_cols names of columns in `nm` that should be used as
#'   covariates
#' @param covariates_implementation named list indicating which implementation
#'   method should be used for each covariate. See [PKPDsim::new_covariate] for
#'   options.
#' @export

nonmem_to_stan_data <- function(
  nm,
  covariate_cols,
  parameters,
  iiv,
  ruv,
  covariates_implementation = list(),
  dose_cmt = 1,
  ltbs = FALSE,
  verbose = FALSE
){
  if (!"ID" %in% colnames(nm)) {
    nm$ID <- 1
  } else if (length(unique(nm[["ID"]])) > 1) {
    stop("Please supply NONMEM data for one patient at a time.")
  }
  stopifnot(all(covariate_cols %in% colnames(nm)))
  stopifnot(all(c("EVID", "TIME", "DV") %in% colnames(nm)))
  stopifnot(any(nm$EVID == 0))
  stopifnot(any(nm$EVID == 1))
  if (!"CMT" %in% colnames(nm)) {
    nm$CMT <- 2
  }
  
  reg <- PKPDsim::nm_to_regimen(nm)
  if (is.null(covariate_cols) || length(covariate_cols) == 0) {
    covs <- NULL
  } else {
    covs <- PKPDsim::covariates_table_to_list(
      nm[c("ID", "TIME", covariate_cols)],
      covariates_implementation
    )[[1]]
  }
  tdms <- nm[nm$EVID == 0, c("TIME", "DV", "CMT")]
  colnames(tdms) <- c("t", "dv", "cmt")
  
  new_stan_data(reg, covs, tdms, parameters, iiv, ruv, dose_cmt, ltbs, verbose) 
}