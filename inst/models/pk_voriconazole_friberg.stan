functions{
  vector ode(real t, vector A, real[] theta, real[] x_r, int[] x_i){
    
    real CL = theta[1];
    real Q = theta[2];
    real V1 = theta[3];
    real V2 = theta[4];
    real KA = theta[5];
    real KM = theta[6];
    real VMAX1 = theta[7];
    
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    
    real VMAXINH = exp(1.5) / (1 + exp(1.5));
    real Vmax = VMAX1 * (1 - VMAXINH * (t-1) / ((t-1) + (2.41 - 1))); // T50 = 2.41
    
    vector[3] dAdt;

    dAdt[1] = -KA*A[1];
    dAdt[2] =  KA*A[1] - (k10 + k12)*A[2] + k21*A[3] - (Vmax * A[2]/V1)/(KM + A[2]/V1);
    dAdt[3] = k12*A[2] - k21*A[3];

    return dAdt;
  }
}

data {
  int<lower = 1> n_t;            // number of events
  int<lower = 1> n_obs;          // number of observations
  int<lower = 1> i_obs[n_obs];    // index of observation
  
  // population parameters
  real<lower = 0> theta_CL;
  real<lower = 0> theta_Q;
  real<lower = 0> theta_V1;
  real<lower = 0> theta_V2;
  real<lower = 0> theta_KA;
  real<lower = 0> theta_KM;
  real<lower = 0> theta_VMAX1;
  
  // inter-individual variability (SD scale)
  real<lower = 0> omega_CL;
  real<lower = 0> omega_Q;
  real<lower = 0> omega_V1;
  real<lower = 0> omega_V2;
  real<lower = 0> omega_KA;
  real<lower = 0> omega_KM;
  real<lower = 0> omega_VMAX1;

  // error model
  int<lower = 0, upper = 1> ltbs; // should log-transform-both-sides be used for observations? (boolean)
  real<lower = 0> ruv_prop;
  real<lower = 0> ruv_add;

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
  
  row_vector<lower = 0>[n_obs] dv;  // observed concentration (dependent variable)
}

transformed data {
  row_vector[n_obs] log_dv = log(dv);
  int n_theta = 7;   // number of parameters
  int n_cmt = 3;   // number of compartments
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
  real<lower = 0> KA;
  real<lower = 0> KM;
  real<lower = 0> VMAX1;
}

transformed parameters {
  real theta[n_theta];
  row_vector<lower = 0>[n_t] ipred;
  row_vector<lower = 0>[n_obs] ipred_obs;
  matrix<lower = 0>[3, n_t] A; 

  theta[1] = CL * pow(mean(WT)/70, 0.75);
  theta[2] = Q * 1.637 * pow(mean(WT)/70, 0.75);
  theta[3] = V1 * mean(WT)/70;
  theta[4] = V2 * mean(WT)/70;
  theta[5] = KA;
  theta[6] = KM;
  theta[7] = VMAX1 * pow(mean(WT)/70, 0.75);
  // theta[8] = F1;
  
  A = pmx_solve_rk45(ode, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  ipred = A[2, ] ./ V1;

  for(i in 1:n_obs){
    ipred_obs[i] = ipred[i_obs[i]];  // predictions for observed data records
  }
}

model{
  // informative prior
  CL ~ lognormal(log(theta_CL), omega_CL);
  Q  ~ lognormal(log(theta_Q), omega_Q);
  V1 ~ lognormal(log(theta_V1), omega_V1);
  V2 ~ lognormal(log(theta_V2), omega_V2);
  KA ~ lognormal(log(theta_KA), omega_KA);
  KM ~ lognormal(log(theta_KM), omega_KM);
  VMAX1 ~ lognormal(log(theta_VMAX1), omega_VMAX1);
  // F1 ~ lognormal(log(theta_F1), omega_F1);
  
  // likelihood for observed data:
  if(ltbs) {
    log_dv ~ normal(log(ipred_obs), ruv_add);
  } else {
    dv ~ normal(ipred_obs, (ruv_prop * ipred_obs + ruv_add));
  }
}

generated quantities{
  real ipred_ruv[n_obs];

  // sample prior
  real prior_CL = lognormal_rng(log(theta_CL), omega_CL); //"parameters", "iiv" values
  real prior_Q  = lognormal_rng(log(theta_Q), omega_Q);
  real prior_V1 = lognormal_rng(log(theta_V1), omega_V1);
  real prior_V2 = lognormal_rng(log(theta_V2), omega_V2);
  real prior_KA = lognormal_rng(log(theta_KA), omega_KA);
  real prior_KM = lognormal_rng(log(theta_KM), omega_KM);
  real prior_VMAX1 = lognormal_rng(log(theta_VMAX1), omega_VMAX1);

  for(i in 1:n_obs){
    ipred_ruv[i] = normal_rng(ipred_obs[i], (ruv_prop * ipred_obs[i] + ruv_add));
  }
}
