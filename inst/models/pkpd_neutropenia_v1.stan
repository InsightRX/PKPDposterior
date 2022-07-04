// PK-PD model for neutropenia based on Friberg et al.
// 
// Drug is fictitious, just a generic linear 2-cmt oral model, which can 
// also be used as 1-cmt and/or iv model.

functions{
  
  vector ode(real t, vector A, real[] theta, real[] x_r, int[] x_i) {
    real CL    = theta[1];
    real Q     = theta[2];
    real V    = theta[3];
    real V2    = theta[4];
    real KA    = theta[5];
    real MTT   = theta[6];
    real CIRC0 = theta[7];
    real GAMMA = theta[8];
    real SLOPE = theta[9];
    
    real k10 = CL / V;
    real k12 = Q / V;
    real k21 = Q / V2;
    real ktr = 4 / MTT;
    
    vector[8] dAdt;
    
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
    
    EDrug = SLOPE * conc; // slope model, not Emax
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

data {
  int<lower = 1> n_t;                  // number of event rows
  int<lower = 1> n_obs_pk;            // number of observations for PK
  int<lower = 1> n_obs_pd;            // number of observations for PD
  int<lower = 1> i_obs_pk[n_obs_pk];  // index of observations PK
  int<lower = 1> i_obs_pd[n_obs_pd];  // index of observations PD
  
  // population parameters
  real<lower = 0> theta_CL;
  real<lower = 0> theta_V;
  real<lower = 0> theta_MTT;
  real<lower = 0> theta_CIRC0;
  real<lower = 0> theta_SLOPE;
  real<lower = 0> theta_GAMMA;
  real<lower = 0> theta_Q;
  real<lower = 0> theta_V2;
  real<lower = 0> theta_KA;
  
  // inter-individual variability (SD scale)
  real<lower = 0> omega_CL;
  real<lower = 0> omega_V;
  real<lower = 0> omega_MTT;
  real<lower = 0> omega_CIRC0;
  real<lower = 0> omega_SLOPE;
  real<lower = 0> omega_GAMMA;
  real<lower = 0> omega_Q;
  real<lower = 0> omega_V2;
  real<lower = 0> omega_KA;

  // error model
  int<lower = 0, upper = 1> ltbs_pk; // should log-transform-both-sides be used for observations? (boolean)
  int<lower = 0, upper = 1> ltbs_pd; // should log-transform-both-sides be used for observations? (boolean)
  real<lower = 0> ruv_prop_pk;
  real<lower = 0> ruv_prop_pd;
  real<lower = 0> ruv_add_pk;
  real<lower = 0> ruv_add_pd;

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

}

parameters{
  real<lower = 0> CL;
  // real<lower = 0> Q;
  real<lower = 0> V;
  // real<lower = 0> V2;
  // real<lower = 0> KA;
  real<lower = 0> MTT;
  real<lower = 0> CIRC0;
  real<lower = 0> SLOPE;
  // real<lower = 0> GAMMA;
}

transformed parameters{
  
  real theta[n_theta];
  row_vector<lower = 0>[n_t]      ipred_pk;
  row_vector<lower = 0>[n_t]      ipred_pd;
  row_vector<lower = 0>[n_obs_pk] ipred_obs_pk;
  row_vector<lower = 0>[n_obs_pd] ipred_obs_pd;
  matrix[n_cmt, n_t] A;
  
  real<lower = 0> GAMMA = theta_GAMMA;
  real<lower = 0> Q = theta_Q;
  real<lower = 0> V2 = theta_V2;
  real<lower = 0> KA = theta_KA;

  theta[1] = CL;
  theta[2] = Q; // Q;
  theta[3] = V;
  theta[4] = V2; // V2;
  theta[5] = KA; // KA;
  theta[6] = MTT;   //  mean transit time
  theta[7] = CIRC0; // baseline level of neutrophils
  theta[8] = GAMMA; // feedback, GAMMA
  theta[9] = SLOPE; // slope

  A = pmx_solve_rk45(ode, n_cmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  ipred_pk = A[2, ] ./ V; // vector will all observation times, for PK
  ipred_pd = A[8, ] + theta[7]; // vector will all observation times, for PK

  ipred_obs_pk = ipred_pk[i_obs_pk]; // vector with only PK observations
  ipred_obs_pd = ipred_pd[i_obs_pd]; // vector with only PD observations
}

model{

  // Priors
  CL    ~ lognormal(log(theta_CL), omega_CL);
  V    ~ lognormal(log(theta_V), omega_V);

  MTT     ~ lognormal(log(theta_MTT), omega_MTT);
  CIRC0   ~ lognormal(log(theta_CIRC0), omega_CIRC0);
  SLOPE   ~ lognormal(log(theta_SLOPE), omega_SLOPE);

  // observed data likelihood
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
  
  vector<lower = 0>[n_obs_pk] ipred_ruv_pk;
  vector<lower = 0>[n_obs_pd] ipred_ruv_pd;

  // Sample from prior:
  // PK
  real prior_CL = lognormal_rng(log(theta_CL), omega_CL);
  real prior_V = lognormal_rng(log(theta_V), omega_V);

  // PD
  real prior_MTT = lognormal_rng(log(theta_MTT), omega_MTT);
  real prior_CIRC0 = lognormal_rng(log(theta_CIRC0), omega_CIRC0);
  real prior_SLOPE = lognormal_rng(log(theta_SLOPE), omega_SLOPE);
  // real prior_GAMMA = lognormal_rng(log(theta_GAMMA), omega_GAMMA);
  
  if(ltbs_pk) {
    for(i in 1:n_obs_pk) {
      ipred_ruv_pk[i] = exp(normal_rng(log(fmax(machine_precision(), ipred_obs_pk[i])), ruv_add_pk)); // ipred with added residual variability for PK
    }
  } else {
    for(i in 1:n_obs_pk) {
      ipred_ruv_pk[i] = normal_rng(ipred_obs_pk[i], (ruv_prop_pk * ipred_obs_pk[i] + ruv_add_pk));
    }
  }
  if(ltbs_pd) {
    for(i in 1:n_obs_pd) {
      ipred_ruv_pd[i] = exp(normal_rng(log(fmax(machine_precision(), ipred_obs_pd[i])), ruv_add_pd)); // ipred with added residual variability for PD
    }
  } else {
    for(i in 1:n_obs_pd) {
      ipred_ruv_pd[i] = normal_rng(ipred_obs_pd[i], (ruv_prop_pd * ipred_obs_pd[i] + ruv_add_pd));
    }
  }
}