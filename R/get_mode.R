#' Simple mode estimation function from density
#' 
#' This function implements a simple estimation of model using the peak of the
#' univariate kernel density. More advanced methods are available in `modeest` 
#' package.
#' 
#' @param x vector of continuous values to estimate mode from
#' @param ... function arguments passed on to [stats::density]
#' 
#' @export
get_mode <- function(x, ...) { # simple mode estimation
  xdens <- stats::density(x, ...)
  xdens$x[which.max(xdens$y)]
}
