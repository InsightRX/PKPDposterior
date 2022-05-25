// Two compartment model using built-in analytical solution

data{
  int<lower = 1> n_t;  // number of events
  int<lower = 1> n_obs_pk;  // number of observation
  int<lower = 1> i_obs_pk[n_obs_pk];  // index of observation
  
  // population parameters
  real<lower = 0> theta_CL;
  real<lower = 0> theta_Q;
  real<lower = 0> theta_V1;
  real<lower = 0> theta_V2;
  real<lower = 0> theta_TH_CRCL;
  real<lower = 0> theta_TDM_INIT; // not used, for compatibility
  
  // inter-individual variability (SD scale)
  real<lower = 0> omega_CL;
  real<lower = 0> omega_Q;
  real<lower = 0> omega_V1;
  real<lower = 0> omega_V2;
  real<lower = 0> omega_TH_CRCL; // not used, for compatibility
  real<lower = 0> omega_TDM_INIT; // not used, for compatibility
 
  // error model
  int<lower = 0, upper = 1> ltbs_pk; // should log-transform-both-sides be used for observations? (boolean)
  real<lower = 0> ruv_prop_pk;
  real<lower = 0> ruv_add_pk;

  // NONMEM data
  int<lower = 1> cmt[n_t];
  int  evid[n_t];
  int  addl[n_t];
  int  ss[n_t];
  real amt[n_t];
  real time[n_t];
  real rate[n_t];
  real ii[n_t];
  real WT[n_t];
  real CRCL[n_t];
  
  vector<lower = 0>[n_obs_pk] dv_pk;  // observed concentration (Dependent Variable)
}

transformed data{
  vector[n_obs_pk] log_dv_pk = log(dv_pk);
  int n_theta = 5;  // number of ODE parameters in Two Compartment Model
  int n_cmt = 3;  // number of compartments in model
}

parameters{
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
}

transformed parameters{
  array[n_t, n_theta] real theta;
  row_vector<lower = 0>[n_t] ipred_pk;
  vector<lower = 0>[n_obs_pk] ipred_obs_pk;
  matrix<lower = 0>[n_cmt, n_t] A;
  real<lower = 0> TH_CRCL = theta_TH_CRCL;
  real<lower = 0> TDM_INIT = theta_TDM_INIT;
  
  for(j in 1:n_t) {
    theta[j, 1] = CL * (1.0 + TH_CRCL * ((CRCL[j] * 16.6667) - 66.0));
    theta[j, 2] = Q;
    theta[j, 3] = V1 * WT[j];
    theta[j, 4] = V2 * WT[j];
    theta[j, 5] = 0;
  }

  // call to analytic solver:
  A = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta);

  // save observations to variables:
  ipred_pk = A[2, :] ./ (V1 * mean(WT)); // predictions for all event records
  ipred_obs_pk = ipred_pk'[i_obs_pk];          // predictions only for observed data records

}

model{
  // likelihood for parameters:
  CL      ~ lognormal(log(theta_CL), omega_CL);
  Q       ~ lognormal(log(theta_Q), omega_Q);
  V1      ~ lognormal(log(theta_V1), omega_V1);
  V2      ~ lognormal(log(theta_V2), omega_V2);
  
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
