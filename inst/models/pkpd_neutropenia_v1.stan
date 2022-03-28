// PK-PD model for neutropenia based on Friberg et al.
// 
// Drug is fictitious, just a generic linear 2-cmt oral model, which can 
// also be used as 1-cmt and/or iv model.

functions{
  
  vector ode(real t, vector A, real[] theta, real[] x_r, int[] x_i) {
    real CL    = theta[1];
    real V1    = theta[3];
    real mtt   = theta[6];
    real circ0 = theta[7];
    real gamma = theta[8];
    real alpha = theta[9];
    // real Q    = theta[2];
    // real V2   = theta[4];
    // real ka   = theta[5];
    
    real ka  = 0;
    real k10 = CL / V1;
    real k12 = 0; // Q / V1;
    real k21 = 0; // Q / V2;
    real ktr = 4 / mtt;
    
    vector[8] dAdt;
    
    real conc;
    real EDrug;
    real transit1;
    real transit2;
    real transit3;
    real circ;
    real prol;
    
    dAdt[1] = -ka * A[1];
    dAdt[2] =  ka * A[1] - (k10 + k12) * A[2] + k21 * A[3];
    dAdt[3] = k12 * A[2] - k21 * A[3];
    conc = A[2] / V1;
    
    EDrug = alpha * conc; // slope model, not Emax
    prol = A[4] + circ0;
    transit1 = A[5] + circ0;
    transit2 = A[6] + circ0;
    transit3 = A[7] + circ0;
    circ = fmax(machine_precision(), A[8] + circ0); // Device for implementing a modeled 
    
    // initial condition
    dAdt[4] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
    dAdt[5] = ktr * (prol - transit1);
    dAdt[6] = ktr * (transit1 - transit2);
    dAdt[7] = ktr * (transit2 - transit3);
    dAdt[8] = ktr * (transit3 - circ);
    
    return dAdt;
  }
}

data {
  int<lower = 1> n_t;                  // number of event rows
  int<lower = 1> n_obs_pk;            // number of observations for PK
  int<lower = 1> n_obs_pd;            // number of observations for PD
  int<lower = 1> i_obs_pk[n_obs_pk];  // index of observations PK
  int<lower = 1> i_obs_pd[n_obs_pd];  // index of observations PD

  // NONMEM data
  int<lower = 1> cmt[n_t];
  int  evid[n_t];
  int  addl[n_t];
  int    ss[n_t];
  real  amt[n_t];
  real time[n_t];
  real rate[n_t];
  real   ii[n_t];
  
  vector<lower = 0>[n_obs_pk] dv_pk;
  vector<lower = 0>[n_obs_pd] dv_pd;
  
}

transformed data{
  vector[n_obs_pk] log_dv_pk = log(dv_pk);
  vector[n_obs_pd] log_dv_pd = log(dv_pd);
  
  int n_theta = 9;  // number of ODE parameters
  int n_cmt = 8;  // number of compartments

  real sigma_pk = 0.1;
  real sigma_pd = 0.3;
  
}

parameters{
  real<lower = 0> CL;
  // real<lower = 0> Q;
  real<lower = 0> V1;
  // real<lower = 0> V2;
  // real<lower = 0> ka;
  real<lower = 0> mtt;
  real<lower = 0> circ0;
  real<lower = 0> alpha;
  real<lower = 0> gamma;
}

transformed parameters{
  
  real theta[n_theta];
  row_vector<lower = 0>[n_t]      ipred_pk;
  row_vector<lower = 0>[n_t]      ipred_pd;
  row_vector<lower = 0>[n_obs_pk] ipred_obs_pk;
  row_vector<lower = 0>[n_obs_pd] ipred_obs_pd;
  matrix[n_cmt, n_t] A;

  theta[1] = CL;
  theta[2] = 0; // Q;
  theta[3] = V1;
  theta[4] = 1; // V2;
  theta[5] = 0; // ka;
  theta[6] = mtt;   //  mean transit time
  theta[7] = circ0; // baseline level of neutrophils
  theta[8] = gamma; // feedback
  theta[9] = alpha; // slope

  A = pmx_solve_rk45(ode, n_cmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  ipred_pk = A[2, ] ./ V1; // vector will all observation times, for PK
  ipred_pd = A[8, ] + theta[7]; // vector will all observation times, for PK

  ipred_obs_pk = ipred_pk[i_obs_pk]; // vector with only PK observations
  ipred_obs_pd = ipred_pd[i_obs_pd]; // vector with only PD observations
}

model{

  // Priors
  CL    ~ lognormal(log(5),   0.2);
  V1    ~ lognormal(log(50),  0.2);
  // Q     ~ lognormal(log(5),   0.2);
  // V2    ~ lognormal(log(100), 0.2);
  // ka    ~ lognormal(log(1),   0.2);

  mtt     ~ lognormal(log(100), 0.2);
  circ0   ~ lognormal(log(5),   0.1);
  alpha   ~ lognormal(log(0.2), 1.0);
  gamma   ~ lognormal(log(0.2), 0.2);

  // observed data likelihood
  log_dv_pk ~ normal(log(ipred_obs_pk), sigma_pk);
  log_dv_pd ~ normal(log(ipred_obs_pd), sigma_pd);

}

generated quantities{
  
  vector<lower = 0>[n_obs_pk] ipred_ruv_pk;
  vector<lower = 0>[n_obs_pd] ipred_ruv_pd;

  // Sample from prior:
  // PK
  real prior_CL = lognormal_rng(log(5), 0.2);
  real prior_V1 = lognormal_rng(log(50), 0.2);
  // real prior_ka = lognormal_rng(log(1), 0.2);
  // real prior_Q = lognormal_rng(log(5), 0.2);
  // real prior_V2 = lognormal_rng(log(100), 0.2);
  
  // PD
  real prior_mtt = lognormal_rng(log(100), 0.2);
  real prior_circ0 = lognormal_rng(log(5), 0.2);
  real prior_alpha = lognormal_rng(log(0.1), 0.2);
  real prior_gamma = lognormal_rng(log(0.2), 0.2);
  
  for(i in 1:n_obs_pk) {
    ipred_ruv_pk[i] = exp(normal_rng(log(fmax(machine_precision(), ipred_obs_pk[i])), sigma_pk)); // ipred with added residual variability for PK
  }
  
  for(i in 1:n_obs_pd) {
    ipred_ruv_pd[i] = exp(normal_rng(log(fmax(machine_precision(), ipred_obs_pd[i])), sigma_pd)); // ipred with added residual variability for PD
  }
}