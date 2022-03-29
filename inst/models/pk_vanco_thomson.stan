// Two compartment model using built-in analytical solution

data{
  int<lower = 1> n_t;  // number of events
  int<lower = 1> n_obs;  // number of observation
  int<lower = 1> i_obs[n_obs];  // index of observation
  
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
  int<lower = 0, upper = 1> ltbs; // should log-transform-both-sides be used for observations? (boolean)
  real<lower = 0> ruv_prop;
  real<lower = 0> ruv_add;

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
  
  vector<lower = 0>[n_obs] dv;  // observed concentration (Dependent Variable)
}

transformed data{
  vector[n_obs] log_dv = log(dv);
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
  row_vector<lower = 0>[n_t] ipred;
  vector<lower = 0>[n_obs] ipred_obs;
  matrix<lower = 0>[n_cmt, n_t] A;
  
  for(j in 1:n_t) {
    theta[j, 1] = CL * (1.0 + 0.0154 * ((CRCL[j] * 16.6667) - 66.0));
    theta[j, 2] = Q;
    theta[j, 3] = V1 * WT[j];
    theta[j, 4] = V2 * WT[j];
    theta[j, 5] = 0;
  }

  // call to analytic solver:
  A = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta);

  // save observations to variables:
  ipred = A[2, :] ./ (V1 * mean(WT)); // predictions for all event records
  ipred_obs = ipred'[i_obs];          // predictions only for observed data records

}

model{
  // likelihood for parameters:
  CL     ~ lognormal(log(theta_CL), omega_CL);
  Q      ~ lognormal(log(theta_Q), omega_Q);
  V1     ~ lognormal(log(theta_V1), omega_V1);
  V2     ~ lognormal(log(theta_V2), omega_V2);
  
  // likelihood for observed data:
  if(ltbs) {
    log_dv ~ normal(log(ipred_obs), ruv_add);
  } else {
    dv ~ normal(ipred_obs, (ruv_prop * ipred_obs + ruv_add));
  }
}

generated quantities{
  real ipred_ruv[n_obs];
  
  // sample prior:
  real prior_CL = lognormal_rng(log(theta_CL), omega_CL);
  real prior_Q  = lognormal_rng(log(theta_Q), omega_Q);
  real prior_V1 = lognormal_rng(log(theta_V1), omega_V1);
  real prior_V2 = lognormal_rng(log(theta_V2), omega_V2);
  
  // posterior:
  for(i in 1:n_obs){
    ipred_ruv[i] = normal_rng(ipred_obs[i], (ruv_prop * ipred_obs[i] + ruv_add));
  }
}
