// Two compartment model using built-in analytical solution

data{
  
  // observation_data  
  vector<lower = 0>[n_obs_pk] dv_pk;  // observed concentration (Dependent Variable)
  
  // observation_variables
  int<lower = 1> n_t;  // number of events
  int<lower = 1> n_obs_pk;  // number of observation
  int<lower = 1> i_obs_pk[n_obs_pk];  // index of observation
  
  // population_parameters
  real<lower = 0> theta_CL;
  real<lower = 0> theta_Q;
  real<lower = 0> theta_V1;
  real<lower = 0> theta_V2;
  
  // iiv_parameters
  real<lower = 0> omega_CL;
  real<lower = 0> omega_Q;
  real<lower = 0> omega_V1;
  real<lower = 0> omega_V2;
  
  // error_parameters
  int<lower = 0, upper = 1> ltbs_pk; // should log-transform-both-sides be used for observations? (boolean)
  real<lower = 0> ruv_prop_pk;
  real<lower = 0> ruv_add_pk;
  
  // input_data
  int<lower = 1> cmt[n_t];
  int  evid[n_t];
  int  addl[n_t];
  int  ss[n_t];
  real amt[n_t];
  real time[n_t];
  real rate[n_t];
  real ii[n_t];
  
  // covariate_input
  real WT[n_t];
  real CRCL[n_t];
  
}

transformed data{
  // log_transform_observations
  vector[n_obs_pk] log_dv_pk = log(dv_pk);
  
  // model_numbers
  int n_theta = 5;  // number of ODE parameters in Two Compartment Model
  int n_cmt = 3;  // number of compartments in model
}

parameters{
  // parameter_definitions
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
}

transformed parameters{
  
  // transformed_parameters
  array[n_t, n_theta] real theta;
  row_vector<lower = 0>[n_t] ipred_pk;
  vector<lower = 0>[n_obs_pk] ipred_obs_pk;
  matrix<lower = 0>[n_cmt, n_t] A;
  
  // pk_block
  for(j in 1:n_t) {
    theta[j, 1] = CL * (1.0 + 0.0154 * ((CRCL[j] * 16.6667) - 66.0));
    theta[j, 2] = Q;
    theta[j, 3] = V1 * WT[j];
    theta[j, 4] = V2 * WT[j];
    theta[j, 5] = 0;
  }
  
  // solver_call
  A = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta);
  
  // ipred_definition
  ipred_pk = A[2, :] ./ (V1 * mean(WT)); // predictions for all event records
  ipred_obs_pk = ipred_pk'[i_obs_pk]; // predictions only for observed data records

}

model{
  
  // likelihood_parameters
  CL     ~ lognormal(log(theta_CL), omega_CL);
  Q      ~ lognormal(log(theta_Q), omega_Q);
  V1     ~ lognormal(log(theta_V1), omega_V1);
  V2     ~ lognormal(log(theta_V2), omega_V2);
  
  // likelihood_observed_data
  if(ltbs_pk) {
    log_dv_pk ~ normal(log(ipred_obs_pk), ruv_add_pk);
  } else {
    dv_pk ~ normal(ipred_obs_pk, (ruv_prop_pk * ipred_obs_pk + ruv_add_pk));
  }
  
}

generated quantities{

  // sample_prior
  real prior_CL = lognormal_rng(log(theta_CL), omega_CL);
  real prior_Q  = lognormal_rng(log(theta_Q),  omega_Q);
  real prior_V1 = lognormal_rng(log(theta_V1), omega_V1);
  real prior_V2 = lognormal_rng(log(theta_V2), omega_V2);
  
  // simulate_posterior_ruv
  real ipred_ruv_pk[n_obs_pk];
  for(i in 1:n_obs_pk){
    ipred_ruv_pk[i] = normal_rng(ipred_obs_pk[i], (ruv_prop_pk * ipred_obs_pk[i] + ruv_add_pk));
  }
  
}
