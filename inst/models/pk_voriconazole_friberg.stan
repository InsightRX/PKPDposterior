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

    return A;
  }
}

data {
  int<lower = 1> n_t;            // number of events
  int<lower = 1> n_obs;          // number of observations
  int<lower = 1> i_obs[n_obs];    // index of observation
  
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
  CL ~ lognormal(log(6.16), 0.5);
  Q  ~ lognormal(log(15.5), 0.42);
  V1 ~ lognormal(log(79), 0.136);
  V2 ~ lognormal(log(103), 0.77);
  KA ~ lognormal(log(1.19), 0.9);
  KM ~ lognormal(log(1.15), 1.0);
  VMAX1 ~ lognormal(log(114), 0.50);
  // F1 ~ lognormal(log(0.585), 1.0);
  
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
  real prior_CL = lognormal_rng(log(6.16), 0.50); //"parameters", "iiv" values
  real prior_Q  = lognormal_rng(log(15.5), 0.42);
  real prior_V1 = lognormal_rng(log(79), 0.136);
  real prior_V2 = lognormal_rng(log(103), 0.77);
  real prior_KA = lognormal_rng(log(1.19), 0.90);
  real prior_KM = lognormal_rng(log(1.15), 1.0);
  real prior_VMAX1 = lognormal_rng(log(114), 0.50);
  // real prior_F1 = lognormal_rng(log(0.585), 1.0);

  for(i in 1:n_obs){
    ipred_ruv[i] = normal_rng(ipred_obs[i], (ruv_prop * ipred_obs[i] + ruv_add));
  }
}
