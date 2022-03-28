// PK-PD model for neutropenia based on Friberg et al.
// Drug is fictitious, just a simple 2-cmt oral model.
functions{
  
  vector two_cmt_oral_neutropenia_ode(real t, vector x, real[] parms, real[] x_r, int[] x_i) {
    real CL = parms[1];
    // real Q = parms[2];
    real V1 = parms[3];
    // real V2 = parms[4];
    // real ka = parms[5];
    real mtt = parms[6];
    real circ0 = parms[7];
    real gamma = parms[8];
    real alpha = parms[9];
    
    real ka  = 0;
    real k10 = CL / V1;
    real k12 = 0; // Q / V1;
    real k21 = 0; // Q / V2;
    real ktr = 4 / mtt;
    
    vector[8] dxdt;
    
    real conc;
    real EDrug;
    real transit1;
    real transit2;
    real transit3;
    real circ;
    real prol;
    
    dxdt[1] = -ka * x[1];
    dxdt[2] =  ka * x[1] - (k10 + k12) * x[2] + k21 * x[3];
    dxdt[3] = k12 * x[2] - k21 * x[3];
    conc = x[2] / V1;
    EDrug = alpha * conc; // slope model, not Emax
    // EDrug = 0 * conc; // slope model, not Emax
    
    prol = x[4] + circ0;
    transit1 = x[5] + circ0;
    transit2 = x[6] + circ0;
    transit3 = x[7] + circ0;
    circ = fmax(machine_precision(), x[8] + circ0); // Device for implementing a modeled 
    
    // initial condition
    dxdt[4] = ktr * prol * ((1 - EDrug) * ((circ0 / circ)^gamma) - 1);
    dxdt[5] = ktr * (prol - transit1);
    dxdt[6] = ktr * (transit1 - transit2);
    dxdt[7] = ktr * (transit2 - transit3);
    dxdt[8] = ktr * (transit3 - circ);
    
    return dxdt;
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
  row_vector<lower = 0>[n_t] ipred_all_pk;
  row_vector<lower = 0>[n_obs_pk] ipred_pk;
  row_vector<lower = 0>[n_t] ipred_all_pd;
  row_vector<lower = 0>[n_obs_pd] ipred_pd;
  matrix[n_cmt, n_t] x;

  theta[1] = CL;
  theta[2] = 0; // Q;
  theta[3] = V1;
  theta[4] = 1; // V2;
  theta[5] = 0; // ka;
  theta[6] = mtt;   //  mean transit time
  theta[7] = circ0; // baseline level of neutrophils
  theta[8] = gamma; // feedback
  theta[9] = alpha; // slope

  x = pmx_solve_rk45(
    two_cmt_oral_neutropenia_ode, 
    8,
    time,
    amt, 
    rate, 
    ii, 
    evid, 
    cmt, 
    addl, 
    ss, 
    theta, 
    1e-5, 1e-8, 1e5
  );

  ipred_all_pk = x[2, ] ./ V1; // vector will all observation times, for PK
  ipred_all_pd = x[8, ] + theta[7]; // vector will all observation times, for PK

  ipred_pk = ipred_all_pk[i_obs_pk]; // vector with only PK observations
  ipred_pd = ipred_all_pd[i_obs_pd]; // vector with only PD observations
}

model{

  // Priors
  CL    ~ lognormal(log(5),   0.2);
  // Q     ~ lognormal(log(5),   0.2);
  V1    ~ lognormal(log(50),  0.2);
  // V2    ~ lognormal(log(100), 0.2);
  // ka    ~ lognormal(log(1),   0.2);

  mtt     ~ lognormal(log(100), 0.2);
  circ0   ~ lognormal(log(5),   0.1);
  alpha   ~ lognormal(log(0.2), 1.0);
  gamma   ~ lognormal(log(0.2), 0.2);

  // observed data likelihood
  log_dv_pk ~ normal(log(ipred_pk), sigma_pk);
  log_dv_pd ~ normal(log(ipred_pd), sigma_pd);

}

generated quantities{
  
  vector<lower = 0>[n_obs_pk] ipred_pk_ruv;
  vector<lower = 0>[n_obs_pd] ipred_pd_ruv;

  // sample prior
  real prior_CL = lognormal_rng(log(5), 0.2);
  real prior_V1 = lognormal_rng(log(50), 0.2);
  // real prior_ka = lognormal_rng(log(1), 0.2);
  // real prior_Q = lognormal_rng(log(5), 0.2);
  // real prior_V2 = lognormal_rng(log(100), 0.2);
  real prior_mtt = lognormal_rng(log(100), 0.2);
  real prior_circ0 = lognormal_rng(log(5), 0.2);
  real prior_alpha = lognormal_rng(log(0.1), 0.2);
  real prior_gamma = lognormal_rng(log(0.2), 0.2);
  
  for(i in 1:n_obs_pk) {
    ipred_pk_ruv[i] = exp(normal_rng(log(fmax(machine_precision(), ipred_pk[i])), sigma_pk)); // ipred with added residual variability for PK
  }
  
  for(i in 1:n_obs_pd) {
    ipred_pd_ruv[i] = exp(normal_rng(log(fmax(machine_precision(), ipred_pd[i])), sigma_pd)); // ipred with added residual variability for PD
  }
}