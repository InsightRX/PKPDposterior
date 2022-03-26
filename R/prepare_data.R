#' Prepare data
#'
#' @param regimen Regimen object (created by [PKPDsim::new_regimen()])
#' @param covariates List of covariate objects created by [PKPDsim::new_covariate()]
#' @param tdm_data Data frame with columns t, dv, and cmt
#' @param dose_cmt Specify what dose compartment. Observation compartment in dataset is irrelevant, handled in model.
#' @return Named list suitable for passing on to Torsten
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
#' prepare_data(regimen, covariates, tdm_data)
prepare_data <- function(
  regimen, 
  covariates, 
  data,
  dose_cmt = 1
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
  out <- as.list(nm_data)

  ## Lowercase some names
  lowercase <- names(out) %in% c("TIME", "EVID", "AMT", "CMT", "SS", "II", "ADDL", "RATE")
  names(out)[lowercase] <- tolower(names(out)[lowercase])

  ## Additional info
  types <- unique(nm_data$TYPE)
  types <- types[!is.na(types)]
  if(!is.null(types)) {
    message(paste0("Parsing multiple observation types: ", paste0(types, collapse = ", ")))
    for(key in types) {
      out[[paste0("cObs", key)]] <- nm_data %>%
        dplyr::filter(.data$TYPE == key) %>%
        dplyr::filter(.data$EVID == 0) %>%
        dplyr::pull(.data$DV)
      out[[paste0("iObs", key)]] <- which(out$evid == 0 & out$TYPE == key)
      out[[paste0("nObs", key)]] <- length(out[[paste0("iObs", key)]])
    }
    out$TYPE <- NULL
  } else {
    out$cObs <- nm_data %>%
      dplyr::filter(.data$EVID == 0) %>%
      dplyr::pull(.data$DV)
    out$iObs <- which(out$evid == 0)
    out$nObs <- length(out$iObs)
  }
  out$nt <- nrow(nm_data)
  
  out
}

#' Convert list of covariate objects to NONMEM formatted data
#'
#' @inheritParams prepare_data
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
#' @inheritParams prepare_data
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
