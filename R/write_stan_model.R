#' Write a Stan model code to file
#' @param code Stan code (vector of character strings).
#' @param file filename to write to. If no name specified, will create temporary
#' file.
#' @return filename of saved model
#' @export
write_stan_model <- function(
  code,
  file = NULL
) {
  if(is.null(file)) {
    file <- tempfile(fileext = ".stan")  
  }
  writeLines(code, con = file)
  
  file
}
