#' Generate a Stan model code from model definitions and template code.
#' This is an internal function. For creating a Stan model, please use 
#' `new_stan_model()`.
#' 
#' @param template_code template model code
#' @param definitions list of model definitions, parsed using 
#' `parse_model_definitions()`.
#' 
generate_stan_code <- function(
  template_code,
  definitions
) {

  cr <- "\n"
  no_indent <- c("ode_function")
  for(key in names(definitions)) {
    idx <- grep(paste0("// ", key), template_code)
    indent <- ifelse(key %in% no_indent, "", "  ")
    if(length(idx) == 1) {
      template_code[[idx]] <- paste0(
        template_code[[idx]], cr,
        paste0(
          indent,
          definitions[[key]], 
          cr, 
          collapse = ""
        ),
        collapse = ""
      )
    } else {
      stop(paste0("Tag for `", key, "` not found in template file."))
    }
  }
  
  template_code
  
}
