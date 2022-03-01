#' Visual representation of percentile estimate
#' 
#' Adapted from PKPDmap::plot_eta()
#' 
#' @param x percentile of parameter or observation in posterior
#' @export
plot_percentile<- function(x) {
  if(length(x) > 1) {
    return(unlist(lapply(x, "plot_percentile")))
  }
  dummy <- "-----"
  x <- x - 0.5
  if(x == 0) {
    return(paste0(dummy, "o", dummy))
  }
  position <- max(min(round(abs(x * 10)), 4), -4)
  scale <- paste0(
    paste0(rep("-", position), collapse=""), 
    "o", 
    paste0(rep("-", max(4-position, 0)), collapse=""), 
    collapse = "")
  out <- paste0(
    dummy, 
    "|",
    scale,
    collapse = ""
  )
  if(x < 0) { # reverse
    splits <- strsplit(out, "")[[1]]
    reversed <- rev(splits)
    out <- paste(reversed, collapse = "")
  }
  out
}
