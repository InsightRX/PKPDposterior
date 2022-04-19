// Two compartment model using built-in analytical solution

// ode_function (optional)
functions{
  vector ode(real t, vector A, real[] theta, real[] x_r, int[] x_i) {
    // define_ode_state
    vector[8] dAdt;

    // define_ode_parameters  
    real CL = theta[1];
    real Q = theta[2];
    real V = theta[3];
    real V2 = theta[4];
    real KA = theta[5];
    real MTT = theta[6];
    real CIRC0 = theta[7];
    real GAMMA = theta[8];
    real ALPHA = theta[9];

    // ode_block
    
  real ka  = 0;
  real k10 = CL / V;
  real k12 = 0;
  real k21 = 0;
  real ktr = 4 / MTT;
  
  real conc;
  real EDrug;
  real transit1;
  real transit2;
  real transit3;
  real circ;
  real prol;
  
  dAdt[1] = -KA * A[1];
  dAdt[2] =  KA * A[1] - (k10 + k12) * A[2] + k21 * A[3];
  dAdt[3] = k12 * A[2] - k21 * A[3];
  conc = A[2] / V;
  
  EDrug = ALPHA * conc; // slope model, not Emax
  prol = A[4] + CIRC0;
  transit1 = A[5] + CIRC0;
  transit2 = A[6] + CIRC0;
  transit3 = A[7] + CIRC0;
  circ = fmax(machine_precision(), A[8] + CIRC0); // Device for implementing a modeled 
  
  // initial condition
  dAdt[4] = ktr * prol * ((1 - EDrug) * ((CIRC0 / circ)^GAMMA) - 1);
  dAdt[5] = ktr * (prol - transit1);
  dAdt[6] = ktr * (transit1 - transit2);
  dAdt[7] = ktr * (transit2 - transit3);
  dAdt[8] = ktr * (transit3 - circ);


    return dAdt;
  }
}


data{
  
  // observation_variables
  int<lower = 1> n_t;  // number of events
  int<lower = 1> n_obs_pk;  // number of observation
  int<lower = 1> n_obs_pd;  // number of observation
  int<lower = 1> i_obs_pk[n_obs_pk];  // index of observation
  int<lower = 1> i_obs_pd[n_obs_pd];  // index of observation


  // observation_data  
  vector<lower = 0>[n_obs_pk] dv_pk;  // observed concentration (Dependent Variable)
  vector<lower = 0>[n_obs_pd] dv_pd;  // observed concentration (Dependent Variable)


  // population_parameters
  real<lower=0> theta_CL;
  real<lower=0> theta_V;
  real<lower=0> theta_SLOPE;
  real<lower=0> theta_MTT;
  real<lower=0> theta_CIRC0;
  real<lower=0> theta_GAMMA;


  // iiv_parameters
  real<lower=0> omega_CL;
  real<lower=0> omega_V;
  real<lower=0> omega_SLOPE;
  real<lower=0> omega_MTT;
  real<lower=0> omega_CIRC0;
  real<lower=0> omega_GAMMA;


  // error_parameters
  int<lower = 0, upper = 1> ltbs_pk;
  real<lower = 0> ruv_prop_pk;
  real<lower = 0> ruv_add_pk;
  int<lower = 0, upper = 1> ltbs_pd;
  real<lower = 0> ruv_prop_pd;
  real<lower = 0> ruv_add_pd;


  // input_data
  int<lower = 1> cmt[n_t];
  int  evid[n_t];
  int  addl[n_t];
  int  ss[n_t];
  real amt[n_t];
  real time[n_t];
  real rate[n_t];
  real ii[n_t];

  
  // covariate_data

}

transformed data{

  // log_transform_observations
  vector[n_obs_pk] log_dv_pk = log(dv_pk);
  vector[n_obs_pd] log_dv_pd = log(dv_pd);


  // model_numbers
  int n_theta = 9;
  int n_cmt = 8;


}

parameters{

  // parameter_definitions
  real<lower=0> CL;
  real<lower=0> V;
  real<lower=0> SLOPE;
  real<lower=0> MTT;
  real<lower=0> CIRC0;
  real<lower=0> GAMMA;


}

transformed parameters{
  
  // transformed_parameters_definitions
  row_vector<lower = 0>[n_t] ipred_pk;
  row_vector<lower = 0>[n_t] ipred_pd;
  vector<lower = 0>[n_obs_pk] ipred_obs_pk;
  vector<lower = 0>[n_obs_pd] ipred_obs_pd;
  matrix[n_cmt, n_t] A;


  // pk_block
  real theta[n_theta];
  for(j in 1:n_t) {
    theta[1] = CL;
    theta[2] = 0;
    theta[3] = V;
    theta[4] = 1;
    theta[5] = 0;
    theta[6] = MTT;
    theta[7] = CIRC0;
    theta[8] = GAMMA;
    theta[9] = SLOPE;
  }


  // solver_call
  A = pmx_solve_rk45(ode, n_cmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);


  // ipred_definition
  ipred_pk = A[2, ] ./ V;;
  ipred_pd = A[8, ] + theta[7];;
  ipred_obs_pk = ipred_pk'[i_obs_pk]; // predictions only for observed data records
  ipred_obs_pd = ipred_pd'[i_obs_pd]; // predictions only for observed data records


}

model{
  
  // likelihood_parameters
  CL ~ lognormal(log(theta_CL), omega_CL);
  V ~ lognormal(log(theta_V), omega_V);
  SLOPE ~ lognormal(log(theta_SLOPE), omega_SLOPE);
  MTT ~ lognormal(log(theta_MTT), omega_MTT);
  CIRC0 ~ lognormal(log(theta_CIRC0), omega_CIRC0);
  GAMMA ~ lognormal(log(theta_GAMMA), omega_GAMMA);


  // likelihood_observed_data
  if(ltbs_pk) {
    log_dv_pk ~ normal(log(ipred_obs_pk), ruv_add_pk); 
  } else {
    dv_pk ~ normal(ipred_obs_pk, (ruv_prop_pk * ipred_obs_pk + ruv_add_pk));
  }
  if(ltbs_pd) {
    log_dv_pd ~ normal(log(ipred_obs_pd), ruv_add_pd); 
  } else {
    dv_pd ~ normal(ipred_obs_pd, (ruv_prop_pd * ipred_obs_pd + ruv_add_pd));
  }


}

generated quantities{

  // sample_prior
  real prior_CL = lognormal_rng(log(theta_CL), omega_CL);
  real prior_V = lognormal_rng(log(theta_V), omega_V);
  real prior_SLOPE = lognormal_rng(log(theta_SLOPE), omega_SLOPE);
  real prior_MTT = lognormal_rng(log(theta_MTT), omega_MTT);
  real prior_CIRC0 = lognormal_rng(log(theta_CIRC0), omega_CIRC0);
  real prior_GAMMA = lognormal_rng(log(theta_GAMMA), omega_GAMMA);


  // simulate_posterior_ruv
  real ipred_ruv_pk[n_obs_pk];
  real ipred_ruv_pd[n_obs_pd];
  for(i in 1:n_obs_pk){
    ipred_ruv_pk[i] = normal_rng(ipred_obs_pk[i], (ruv_prop_pk * ipred_obs_pk[i] + ruv_add_pk));
  }
  for(i in 1:n_obs_pd){
    ipred_ruv_pd[i] = normal_rng(ipred_obs_pd[i], (ruv_prop_pd * ipred_obs_pd[i] + ruv_add_pd));
  }


}
