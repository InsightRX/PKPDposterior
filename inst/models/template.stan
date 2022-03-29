// Two compartment model using built-in analytical solution

data{
  
  // observation_data  

  // observation_variables

  // population_parameters

  // iiv_parameters

  // error_parameters

  // input_data
  
  // covariate_data

}

transformed data{

  // log_transform_observations

  // model_numbers

}

parameters{

  // parameter_definitions

}

transformed parameters{
  
  // transformed_parameters

  // pk_block

  // solver_call

  real scale = theta[3];
  
  // ipred_definition

}

model{
  
  // likelihood_parameters

  // likelihood_observed_data

}

generated quantities{

  // sample_prior

  // simulate_posterior_ruv

}
