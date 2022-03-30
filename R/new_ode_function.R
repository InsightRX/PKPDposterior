#' Generates the ODE function that is used in a Stan/Torsten model
#' 
#' @param ode Stan code for ODE block
#' @param parameters list or vector of parameters. Needs to be in order!
#' @param n_cmt number of compartments in the ODE
#' 
#' @returns Stan code with ODE function
#' 
#' @export
new_ode_function <- function(
  ode,
  parameters,
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
  if(class(parameters) == "list") {
    par_names <- names(parameters)
  } else if(class(parameters) == "character") {
    par_names <- parameters
  } else {
    stop("`parameters` needs to be either a list or a character vector.")
  }
  def[["define_ode_parameters"]] <- paste0(
    "  real ", par_names, " = theta[", 1:length(par_names), "];"
  )
  def[["ode_block"]] <- paste(" ", ode)
  
  stan_code <- generate_stan_code(
    template_code = code,
    definitions = def
  )

  stan_code
  
}