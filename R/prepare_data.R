#' Prepare data
#'
#' @param regimen Regimen object (created by [PKPDsim::new_regimen()])
#' @param covariates List of covariate objects created by [PKPDsim::new_covariate()]
#' @param tdm_data Data frame with columns t, dv, and cmt
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
#'   WT = PKPDsim::new_covariate(value = 150, unit = "kg"),
#'   CRCL = PKPDsim::new_covariate(value = 6.5, unit = "l/hr")
#' )
#' tdm_data <- data.frame(
#'   t = c(1, 2), 
#'   dv = c(900, 800), 
#'   cmt = c(1, 1)
#' )
#' prepare_data(regimen, covariates, tdm_data)
prepare_data <- function(regimen, covariates, tdm_data) {
  ## Convert regimen, covariates, tdm data
  reg <- regimen_to_nm(regimen)
  cov <- covariates_to_nm(covariates)
  tdm <- tdm_to_nm(tdm_data)

  ## Combine into one dataset
  nm_data <- rbind(reg, tdm)
  nm_data <- merge(nm_data, cov, by = "ID")
  nm_data <- nm_data[order(nm_data$ID, nm_data$TIME), ]
  out <- as.list(nm_data)

  ## Lowercase some names
  lowercase <- names(out) %in% c("TIME", "EVID", "AMT", "CMT", "SS", "II", "ADDL", "RATE")
  names(out)[lowercase] <- tolower(names(out)[lowercase])
  names(out)[names(out) == "DV"] <- "cObs"

  ## Not using ADDL, SS, II
  out$addl <- 0
  out$ss <- 0
  out$ii <- 0

  ## Additional info
  out$nt <- nrow(nm_data)
  out$iObs <- which(out$evid == 0)
  out$nObs <- length(out$iObs)

  out
}

#' Convert list of covariate objects to NONMEM formatted data
#'
#' @inheritParams prepare_data
covariates_to_nm <- function(covariates) {
  dat <- data.frame(lapply(covariates, `[[`, "value"))
  dat$ID <- 1
  dat
}

#' Convert TDM data to NONMEM format
#'
#' @inheritParams prepare_data
tdm_to_nm <- function(tdm_data) {
  tdm_data$EVID <- 0
  tdm_data$ID <- 1
  tdm_data$MDV <- 0
  tdm_data$AMT <- 0
  names(tdm_data)[names(tdm_data) == "t"] <- "TIME"
  names(tdm_data) <- toupper(names(tdm_data))
  tdm_data
}
