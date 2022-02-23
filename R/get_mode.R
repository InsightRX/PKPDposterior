#' Simple mode estimation function from density
#' 
#' More advanced methods are available in `modeest` package,
#' but likely not necessary for our purpose.
#' 
#' @param x vector of continuous values to estimate mode from
#' 
#' @export
get_mode <- function(x) { # simple mode estimation
  xdens <- density(x)
  xdens$x[which.max(xdens$y)]
}
