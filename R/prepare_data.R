#' Prepare data
#'
#' @param dat Output from [irxnlme::irxanalytics_to_nm_data()]
#' @return Named list suitable for passing on
#' @export
#' @examples
#' dat <- data.frame(
#'   ID = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
#'   OCC = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
#'   TIME = c(0, 19.1234, 42.5678, 43.1958, 66.9230, 78.5511, 78.5514, 91.1947),
#'   EVID = c(1, 1, 0, 1, 1, 2, 2, 0),
#'   MDV = c(1, 1, 0, 1, 1, 1, 1, 0),
#'   DV = c(0, 0, 8.4, 0, 0, 0, 0, 12.6),
#'   AMT = c(1500, 1000, 0, 1250, 1250, 0, 0, 0),
#'   LENGTH = c(1.5, 1, 0, 1, 1, 0, 0, 0),
#'   RATE = c(1000, 1000, 0, 1250, 1250, 0, 0, 0),
#'   SEX = c(1, 1, 1, 1, 1, 1, 1, 1),
#'   AGE = 82,
#'   WEIGHT = 77,
#'   CREAT = 1,
#'   HEIGHT = 150
#' )
#' prepare_data(dat)
prepare_data <- function(dat) {
  out <- as.list(dat)

  ## Lowercase some names
  lowercase <- names(out) %in% c("TIME", "EVID", "AMT", "CMT", "SS", "II", "ADDL", "RATE")
  names(out)[lowercase] <- tolower(names(out)[lowercase])
  names(out)[names(out) == "DV"] <- "cObs"

  ## Not using ADDL, SS, II
  out$addl <- 0
  out$ss <- 0
  out$ii <- 0

  ## Additional info
  out$nt <- nrow(dat)
  out$iObs <- which(out$evid == 0)
  out$nObs <- length(out$iObs)

  out
}
