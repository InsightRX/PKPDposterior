#' Check that the supplied input data for Stan is OK for the model
#' 
#' @param code Stan model code
#' @param data input data, generated using `new_stan_data()`, which includes
#' an element `stan_data` that is checked in this function.
#' 
check_stan_input_data <- function(
  code,
  data
) {
  
  ## Check that all observation types defined in model are supplied in input data:
  obs_types_mod <- lapply(
    code[grep("\\[n_obs_(.*)\\]", code)], 
    function(x) { gsub(".*\\[n_obs_(.*)\\].*", "\\1", x) }
  ) %>% 
    unlist() %>% 
    unique()
  obs_types_data <- lapply(
    names(data$stan_data)[grep("n_obs_(.*)", names(data$stan_data))], 
    function(x) { gsub("n_obs_(.*)", "\\1", x) }
  ) %>% 
    unlist() %>% 
    unique()

  if(!all(obs_types_mod %in% obs_types_data)) {
    stop(
      "Not all observation types defined in the model were supplied in the input data. Missing: ",
      paste0(setdiff(obs_types_mod, obs_types_data), collapse = " ")
    )
  }
  if(!all(obs_types_data %in% obs_types_mod)) {
    message(
      "The input data has observation type(s) that are not defined in model. Undefined types: ", 
      paste0(setdiff(obs_types_data, obs_types_mod), collapse = " "), 
      ". Sampling will continue but undefined observation types are not included in the likelihood."
    )
  }
  
}