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
  int<lower = 1> n_obs;          // number of observations
  int<lower = 1> i_obs[n_obs];    // index of observation
  
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
  
  row_vector<lower = 0>[n_obs] dv;  // observed concentration (dependent variable)
}

transformed data {
  row_vector[n_obs] log_dv = log(dv);
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
  row_vector<lower = 0>[n_t] ipred;
  row_vector<lower = 0>[n_obs] ipred_obs;
  matrix<lower = 0>[3, n_t] A; 

  theta[1] = CL * (1.0 + 0.0154 * ((mean(CRCL) * 16.6667) - 66.0));;
  theta[2] = Q;
  theta[3] = V1 * mean(WT);
  theta[4] = V2 * mean(WT);
  theta[5] = 0;

  A = pmx_solve_rk45(ode, n_cmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  ipred = A[2, ] ./ (V1 * mean(WT));

  for(i in 1:n_obs){
    ipred_obs[i] = ipred[i_obs[i]];  // predictions for observed data records
  }
}

model{
  // informative prior
  CL ~ lognormal(log(2.99), 0.27);
  Q ~ lognormal(log(2.28), 0.49);
  V1 ~ lognormal(log(0.675), 0.15);
  V2 ~ lognormal(log(0.732), 1.3);
  
  dv ~ normal(ipred_obs, (0.15 * ipred_obs + 1.6));
}

generated quantities{
  real ipred_ruv[n_obs];
  
  // sample prior:
  real prior_CL = lognormal_rng(log(2.99),  0.27);
  real prior_Q  = lognormal_rng(log(2.28),  0.49);
  real prior_V1 = lognormal_rng(log(0.675), 0.15);
  real prior_V2 = lognormal_rng(log(0.732), 1.3);
  
  // posterior:
  for(i in 1:n_obs){
    ipred_ruv[i] = normal_rng(ipred_obs[i], (0.15 * ipred_obs[i] + 1.6));
  }
}
