#' Prepare data object for use in get_mcmc_posterior()
#'
#' @param regimen Regimen object (created by [PKPDsim::new_regimen()])
#' @param covariates List of covariate objects created by 
#'   [PKPDsim::new_covariate()]
#' @param data Data frame with columns `t`, `dv`, and `cmt`
#' @param dose_cmt Specify what dose compartment. Observation compartment in 
#' dataset is irrelevant, handled in model.
#' @param parameters list of population parameters, e.g. `list(CL = 5, V = 50)`
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
#' @return Named list suitable for passing on to Torsten.
#' @export
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
#' PKPDsim_to_stan_data(
#'   regimen, 
#'   covariates, 
#'   tdm_data,
#'   parameters = list(CL = 5, V = 50),
#'   iiv = list(CL = 0.1, V = 0.2),
#'   ruv = list(prop = 0.1, add = 1)
#' )

PKPDsim_to_stan_data <- function(
  regimen, 
  covariates, 
  data,
  parameters,
  iiv,
  ruv,
  dose_cmt = 1,
  ltbs = FALSE,
  verbose = FALSE
) {
  ## Convert regimen, covariates, tdm data
  reg <- regimen_to_nm(
    regimen, 
    dose_cmt = dose_cmt
  )
  cov <- covariates_to_nm(covariates)
  
  ## Parse observed data
  # sets <- list()
  # if(!is.null(data$type)) {
  #   types <- unique(data$type) 
  #   for(key in types) {
  #     sets[[key]] <- tdm_to_nm(data[data$type == "key",])
  #     obs <- bind_rows()
  #   }
  # } else {
  #   
  # }
  obs <-  tdm_to_nm(data)

  ## Combine into one dataset
  nm_data <- dplyr::bind_rows(reg, obs, cov) %>%
    dplyr::arrange(.data$ID, .data$TIME, .data$EVID) %>%
    ## Fill in timevarying covariates
    tidyr::fill(names(covariates), .direction = "downup")

  ## Not using ADDL, SS, II
  nm_data$addl <- 0
  nm_data$ss <- 0
  nm_data$ii <- 0
  if(is.null(nm_data$TYPE)) nm_data$TYPE <- "pk" # force an observation type
  out <- as.list(nm_data)

  ## Lowercase some names
  lowercase <- intersect(
    names(out),
    c("TIME", "EVID", "AMT", "CMT", "SS", "II", "ADDL", "RATE")
  )
  names(out)[which(names(out) %in% lowercase)] <- tolower(lowercase)

  ## Population parameters and IIV
  for(key in names(parameters)) {
    out[[paste0("theta_", key)]] <- parameters[[key]]
    if(is.null(iiv[[key]])) {
      stop("`iiv` object requires same list elements as `parameters` object.")
    }
    out[[paste0("omega_", key)]] <- iiv[[key]]
  }
  
  ## Additional info
  types <- setdiff(out$TYPE, NA)
  if(verbose) {
    message("Parsing observation types: ", paste0(types, collapse = ", "))
  }
  if(class(ltbs) == "logical") { # make it a list of lists
    if(length(types) > 1) {
      message("Assuming `ltbs=", dput(ltbs), "` for all observation types.")
    }
    ltbs_list <- list()
    for(key in types) ltbs_list[[key]] <- ltbs
    ltbs <- ltbs_list
  }
  if(class(ruv[[1]]) != "list") { # make it a list of lists
    ruv_list <- list()
    for(key in types) ruv_list[[key]] <- ruv
    ruv <- ruv_list
  }
  if(!all(types %in% names(ruv))) {
    stop(
      "With multiple observation types, ",
      "`ruv` must be specified as a list of lists for each observation type."
    )
  }
  for (key in types) {
    out[[paste0("dv_", key)]] <- nm_data %>%
      dplyr::filter(.data$TYPE == key) %>%
      dplyr::filter(.data$EVID == 0) %>%
      dplyr::pull(.data$DV)
    out[[paste0("i_obs_", key)]] <- which(out$evid == 0 & out$TYPE == key)
    out[[paste0("n_obs_", key)]] <- length(out[[paste0("i_obs_", key)]])
    
    ## error model:
    if (ltbs[[key]] &&
        (is.null(ruv[[key]]$add) || !is.null(ruv[[key]]$prop))) {
      stop(
        "With LTBS, additive error magnitude needs to be specified ", 
        "and proportional error cannot be specified. "
      )
    }
    out[[paste0("ruv_prop_", key)]] <-
      ifelse(!is.null(ruv[[key]]$prop), ruv[[key]]$prop, 0)
    out[[paste0("ruv_add_", key)]] <-
      ifelse(!is.null(ruv[[key]]$add), ruv[[key]]$add, 0)
    out[[paste0("ltbs_", key)]] <- ltbs[[key]]
  }
  out$TYPE <- NULL
  out$n_t <- nrow(nm_data)
  
  out
}

#' Convert list of covariate objects to NONMEM formatted data
#'
#' @inheritParams PKPDsim_to_stan_data
covariates_to_nm <- function(covariates) {
  if(is.null(covariates)) {
    return(NULL)
  }
  dat <- dplyr::bind_rows(covariates, .id = "cov") %>%
    dplyr::select(.data$cov, .data$value, .data$times) %>%
    tidyr::pivot_wider(names_from = .data$cov, values_from = .data$value) %>%
    dplyr::rename(TIME = .data$times)
  dat$ID <- 1
  dat$EVID <- 2
  dat$MDV <- 0
  dat$AMT <- 0
  dat$RATE <- 0
  dat$DV <- 0
  dat$CMT <- 2 # TODO: don't hard code
  dat
}

#' Convert TDM data to NONMEM format
#'
#' @inheritParams PKPDsim_to_stan_data
tdm_to_nm <- function(data) {
  if(is.null(data$cmt)) {
    data$CMT <- 1 # irrelevant, handled in Stan model
  }
  data$EVID <- 0
  data$ID <- 1
  data$MDV <- 0
  data$AMT <- 0
  data$RATE <- 0
  names(data)[names(data) == "t"] <- "TIME"
  names(data) <- toupper(names(data))
  data
}