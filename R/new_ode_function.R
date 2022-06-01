#' Generates the ODE function that is used in a Stan/Torsten model
#' 
#' @param ode Stan code for ODE block
#' @param parameter_names vector of parameter names. Needs to be in order defined by the 
#' chosen Torsten analytic solver.
#' @param n_cmt number of compartments in the ODE
#' 
#' @return Character vector containing valid Stan code with ODE function
#' 
#' @export
new_ode_function <- function(
  ode,
  parameter_names,
  n_cmt
) {
  
  code <- c(
    "functions{",
    "  vector ode(real t, vector A, real[] theta, real[] x_r, int[] x_i) {",
    "    // define_ode_state",
    "    // define_ode_parameters  ",
    "    // ode_block",
    "    return dAdt;",
    "  }",
    "}")
  
  def <- list()
  def[["define_ode_state"]] <- paste0("  vector[", n_cmt, "] dAdt;")
  if(class(parameter_names) != "character") {
    stop("`parameters` needs to be either a list or a character vector.")
  }
  def[["define_ode_parameters"]] <- paste0(
    "  real ", parameter_names, " = theta[", 1:length(parameter_names), "];"
  )
  def[["ode_block"]] <- paste(" ", ode)
  
  generate_stan_code(
    template_code = code,
    definitions = def
  )
  
}
