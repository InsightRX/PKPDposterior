#' Get quantile for a given new value, given a vector of values.
#' 
#' @param x value to calculate quantile for
#' @param y vector of values, posterior draws in most cases
#' 
#' @export
get_quantile <- function(x, y) {
  ecdf(y)(x)
}
