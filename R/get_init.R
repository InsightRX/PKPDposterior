#' Get prior for Torsten from a PKPDsim model package
#' 
#' @param model package name
#' @export
get_init <- function(
  model
) {
  model <- gsub("_", "", model)
  par <- get("parameters", asNamespace(model))()
  ruv <- get("ruv", asNamespace(model))()
  par$sigma <- ruv$prop
  par$TDM_INIT <- NULL
  attr(par, "units") <- NULL
  par
}
