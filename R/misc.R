#' Get quantile for a given new value, given a vector of values.
#' 
#' @param x value to calculate quantile for
#' @param y vector of values, posterior draws in most cases
#' @keywords internal
get_quantile <- function(x, y) {
  stats::ecdf(y)(x)
}
