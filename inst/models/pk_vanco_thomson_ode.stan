functions{
  
  vector ode (real t, vector A, real[] theta, real[] x_r, int[] x_i){

    real CL = theta[1];
    real Q  = theta[2];
    real V1 = theta[3];
    real V2 = theta[4];
    real ka = 0;
    
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    
    vector[3] dAdt;

    dAdt[1] = -ka*A[1];
    dAdt[2] =  ka*A[1] - (k10 + k12)*A[2] + k21*A[3];
    dAdt[3] =  k12*A[2] - k21*A[3];

    return dAdt;
  }
}

data {
  int<lower = 1> n_t;            // number of events
  int<lower = 1> n_obs_pk;          // number of observations
  int<lower = 1> i_obs_pk[n_obs_pk];    // index of observation
  
  // population parameters
  real<lower = 0> theta_CL;
  real<lower = 0> theta_Q;
  real<lower = 0> theta_V1;
  real<lower = 0> theta_V2;
  
  // inter-individual variability (SD scale)
  real<lower = 0> omega_CL;
  real<lower = 0> omega_Q;
  real<lower = 0> omega_V1;
  real<lower = 0> omega_V2;
  
  // error model
  int<lower = 0, upper = 1> ltbs_pk; // should log-transform-both-sides be used for observations? (boolean)
  real<lower = 0> ruv_prop_pk;
  real<lower = 0> ruv_add_pk;

  // NONMEM data
  int<lower = 1> cmt[n_t];
  int  evid[n_t];
  int  addl[n_t];
  int    ss[n_t];
  real  amt[n_t];
  real time[n_t];
  real rate[n_t];
  real   ii[n_t];
  real   WT[n_t];
  real CRCL[n_t];
  
  row_vector<lower = 0>[n_obs_pk] dv_pk;  // observed concentration (dependent variable)
}

transformed data {
  row_vector[n_obs_pk] log_dv_pk = log(dv_pk);
  int n_theta = 5;   // number of parameters
  int n_cmt = 3;     // number of compartments
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
}

transformed parameters {
  real theta[n_theta];
  row_vector<lower = 0>[n_t] ipred_pk;
  row_vector<lower = 0>[n_obs_pk] ipred_obs_pk;
  matrix<lower = 0>[3, n_t] A; 

  theta[1] = CL * (1.0 + 0.0154 * ((mean(CRCL) * 16.6667) - 66.0));;
  theta[2] = Q;
  theta[3] = V1 * mean(WT);
  theta[4] = V2 * mean(WT);
  theta[5] = 0;

  A = pmx_solve_rk45(ode, n_cmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  ipred_pk = A[2, ] ./ theta[3];

  for(i in 1:n_obs_pk){
    ipred_obs_pk[i] = ipred_pk[i_obs_pk[i]];  // predictions for observed data records
  }
}

model{
  // likelihood for parameters:
  CL     ~ lognormal(log(theta_CL), omega_CL);
  Q      ~ lognormal(log(theta_Q), omega_Q);
  V1     ~ lognormal(log(theta_V1), omega_V1);
  V2     ~ lognormal(log(theta_V2), omega_V2);
  
  // likelihood for observed data:
  if(ltbs_pk) {
    log_dv_pk ~ normal(log(ipred_obs_pk), ruv_add_pk);
  } else {
    dv_pk ~ normal(ipred_obs_pk, (ruv_prop_pk * ipred_obs_pk + ruv_add_pk));
  }
}

generated quantities{
  real ipred_ruv_pk[n_obs_pk];
  
  // sample prior:
  real prior_CL = lognormal_rng(log(theta_CL), omega_CL);
  real prior_Q  = lognormal_rng(log(theta_Q), omega_Q);
  real prior_V1 = lognormal_rng(log(theta_V1), omega_V1);
  real prior_V2 = lognormal_rng(log(theta_V2), omega_V2);
  
  // posterior:
  for(i in 1:n_obs_pk){
    ipred_ruv_pk[i] = normal_rng(ipred_obs_pk[i], (ruv_prop_pk * ipred_obs_pk[i] + ruv_add_pk));
  }
}
