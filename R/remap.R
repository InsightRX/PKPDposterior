#' Remap named list or data.frame using a translation map
#' @param x names list or data.frame
#' @param map map used to translate elements of x, e.g. `list("CLi" = "CL", "V"
#'   = "V1")`
#' @param reverse interpret map as reverse
#' 
#' @export
remap <- function(x, map, reverse = FALSE) {
  for(key in names(map)) {
    if(reverse) {
      x[[key]] <- x[[map[[key]]]]
      x[[map[[key]]]] <- NULL
    } else {
      x[[map[[key]]]] <- x[[key]]
      x[[key]] <- NULL
    }
  }
  x
}
