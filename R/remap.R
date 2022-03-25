#' Remap named list or data.frame using a translation map
#' @param x names list or data.frame
#' @param map map used to translate elements of x, e.g. `list("CLi" = "CL", "V" = "V1")`
#' 
#' @export
remap <- function(x, map) {
  for(key in names(map)) {
    x[[key]] <- x[[map[[key]]]]
    x[[map[[key]]]] <- NULL
  }
  x
}
