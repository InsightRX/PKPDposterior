// The bioavailability & absorption lag model described by the Friberg et al 
// (AAC, 2012) publication has not been implemented. Please use this only for 
// infusion medication administration simulations.
functions{
  vector ode(real t, vector A, real[] theta, real[] x_r, int[] x_i){
    
    // define_ode_state
    vector[3] dAdt;
    
    // define_ode_parameters  
    real CL = theta[1];
    real Q = theta[2];
    real V = theta[3];
    real V2 = theta[4];
    real KA = theta[5];
    real KM = theta[6];
    real VMAX1 = theta[7];
    real T50 = theta[8];
    
    real k10 = CL / V;
    real k12 = Q / V;
    real k21 = Q / V2;
    real VMAXINH = exp(1.5) / (1 + exp(1.5));
    real Vmax = VMAX1 * (1 - VMAXINH * (t-1) / ((t-1) + (T50 - 1)));
    dAdt[1] = -KA*A[1];
    dAdt[2] =  KA*A[1] - (k10 + k12)*A[2] + k21*A[3] - (Vmax * A[2]/V)/(KM + A[2]/V);
    dAdt[3] = k12*A[2] - k21*A[3];

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
  real<lower = 0> theta_V;
  real<lower = 0> theta_V2;
  real<lower = 0> theta_KA;
  real<lower = 0> theta_KM;
  real<lower = 0> theta_VMAX1;
  real<lower = 0> theta_T50;
  real<lower = 0> theta_TLAG; // not used, for compatibility
  real<lower = 0> theta_F1; // not used, for compatibility
  real<lower = 0> theta_BCF; // not used, for compatibility
  
  // inter-individual variability (SD scale)
  real<lower = 0> omega_CL;
  real<lower = 0> omega_Q;
  real<lower = 0> omega_V;
  real<lower = 0> omega_V2;
  real<lower = 0> omega_KA;
  real<lower = 0> omega_KM;
  real<lower = 0> omega_VMAX1;
  real<lower = 0> omega_T50; // not used, for compatibility
  real<lower = 0> omega_TLAG; // not used, for compatibility
  real<lower = 0> omega_F1; // not used, for compatibility
  real<lower = 0> omega_BCF; // not used, for compatibility

  // error model
  int<lower = 0, upper = 1> ltbs_pk; // should log-transform-both-sides be used for observations? (boolean)
  real<lower = 0> ruv_prop_pk;
  real<lower = 0> ruv_add_pk;

  // NONMEM data
  int<lower = 1> cmt[n_t];
  int evid[n_t];
  int addl[n_t];
  int ss[n_t];
  real amt[n_t];
  real time[n_t];
  real rate[n_t];
  real ii[n_t];
  
  // covariates
  real WT[n_t];
  
  row_vector<lower = 0>[n_obs_pk] dv_pk;  // observed concentration (dependent variable)
}

transformed data {
  row_vector[n_obs_pk] log_dv_pk = log(dv_pk);
  int n_theta = 8;   // number of parameters
  int n_cmt = 3;   // number of compartments
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V;
  real<lower = 0> V2;
  real<lower = 0> KA;
  real<lower = 0> KM;
  real<lower = 0> VMAX1;
}

transformed parameters {
  real theta[n_theta];
  row_vector<lower = 0>[n_t] ipred_pk;
  row_vector<lower = 0>[n_obs_pk] ipred_obs_pk;
  matrix<lower = 0>[3, n_t] A; 
  real<lower = 0> T50 = theta_T50;
  real<lower = 0> TLAG = theta_TLAG;
  real<lower = 0> F1 = theta_F1;
  real<lower = 0> BCF = theta_BCF;

  theta[1] = CL * pow(mean(WT)/70, 0.75);
  theta[2] = Q * 1.637 * pow(mean(WT)/70, 0.75);
  theta[3] = V * mean(WT)/70;
  theta[4] = V2 * mean(WT)/70;
  theta[5] = KA;
  theta[6] = KM;
  theta[7] = VMAX1 * pow(mean(WT)/70, 0.75);
  theta[8] = T50;
  // theta[8] = F1;
  
  A = pmx_solve_rk45(ode, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  ipred_pk = A[2, ] ./ (V * mean(WT)/70);

  for(i in 1:n_obs_pk){
    ipred_obs_pk[i] = ipred_pk[i_obs_pk[i]];  // predictions for observed data records
  }
}

model{
  // informative prior
  CL ~ lognormal(log(theta_CL), omega_CL);
  Q  ~ lognormal(log(theta_Q), omega_Q);
  V ~ lognormal(log(theta_V), omega_V);
  V2 ~ lognormal(log(theta_V2), omega_V2);
  KA ~ lognormal(log(theta_KA), omega_KA);
  KM ~ lognormal(log(theta_KM), omega_KM);
  VMAX1 ~ lognormal(log(theta_VMAX1), omega_VMAX1);
  // F1 ~ lognormal(log(theta_F1), omega_F1);
  
  // likelihood for observed data:
  if(ltbs_pk) {
    log_dv_pk ~ normal(log(ipred_obs_pk), ruv_add_pk);
  } else {
    dv_pk ~ normal(ipred_obs_pk, (ruv_prop_pk * ipred_obs_pk + ruv_add_pk));
  }
}

generated quantities{
  real ipred_ruv_pk[n_obs_pk];

  // sample prior
  real prior_CL = lognormal_rng(log(theta_CL), omega_CL); //"parameters", "iiv" values
  real prior_Q  = lognormal_rng(log(theta_Q), omega_Q);
  real prior_V = lognormal_rng(log(theta_V), omega_V);
  real prior_V2 = lognormal_rng(log(theta_V2), omega_V2);
  real prior_KA = lognormal_rng(log(theta_KA), omega_KA);
  real prior_KM = lognormal_rng(log(theta_KM), omega_KM);
  real prior_VMAX1 = lognormal_rng(log(theta_VMAX1), omega_VMAX1);

  for(i in 1:n_obs_pk){
    ipred_ruv_pk[i] = normal_rng(ipred_obs_pk[i], (ruv_prop_pk * ipred_obs_pk[i] + ruv_add_pk));
  }
}
