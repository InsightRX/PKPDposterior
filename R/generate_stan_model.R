#' Generate a Stan model from model definitions and template code.
#' This is an internal function. For creating a Stan model, please use 
#' `new_stan_model()`.
#' 
#' @param template_code template model code
#' @param definitions list of model definitions, parsed using 
#' `parse_model_definitions()`.
#' 
generate_stan_model <- function(
  template_code,
  definitions
) {

  cr <- "\n"
  for(key in names(definitions)) {
    idx <- grep(paste0("// ", key), template_code)
    if(length(idx) > 0) {
      template_code[[idx]] <- paste0(
        template_code[[idx]], cr,
        paste0(
          "  ", 
          definitions[[key]], 
          cr, 
          collapse = ""
        ),
        collapse = ""
      )
    } else {
      stop(paste0("Definition for `", key, "` not specified."))
    }
  }
  
  template_code
  
}
