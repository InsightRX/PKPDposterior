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
  int<lower = 1> nt;                // number of event rows
  int<lower = 1> nObsPK;            // number of observations
  int<lower = 1> nObsPD;            // number of observations
  int<lower = 1> iObsPK[nObsPK];    // index of observation
  int<lower = 1> iObsPD[nObsPD];    // index of observation

  // NONMEM data
  int<lower = 1> cmt[nt];
  int  evid[nt];
  int  addl[nt];
  int  ss[nt];
  real amt[nt];
  real time[nt];
  real rate[nt];
  real ii[nt];
  
  vector<lower = 0>[nObsPK] cObsPK;
  vector<lower = 0>[nObsPD] cObsPD;
  
}

transformed data{
  vector[nObsPK] logCObsPK = log(cObsPK);
  vector[nObsPD] logCObsPD = log(cObsPD);
  
  int nTheta = 9;  // number of ODE parameters
  int nCmt = 8;  // number of compartments

  real sigmaPK = 0.1;
  real sigmaPD = 0.3;
  
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
  
  real theta[nTheta];
  row_vector<lower = 0>[nt] cHatPK;
  row_vector<lower = 0>[nObsPK] cHatPKObs;
  row_vector<lower = 0>[nt] cHatPD;
  row_vector<lower = 0>[nObsPD] cHatPDObs;
  matrix[nCmt, nt] x;

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

  cHatPK = x[2, ] ./ V1;
  cHatPD = x[8, ] + theta[7];

  cHatPKObs = cHatPK[iObsPK];
  cHatPDObs = cHatPD[iObsPD];
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
  logCObsPK ~ normal(log(cHatPKObs), sigmaPK);
  logCObsPD ~ normal(log(cHatPDObs), sigmaPD);
  // cObsPD ~ normal(cHatPDObs, sigmaPD);

}

generated quantities{
  
  matrix[nCmt, nt] xPred;
  row_vector[nt] cHatPKPred;
  row_vector[nt] cHatPDPred;
  vector<lower = 0>[nObsPK] cHatPKObsCond;
  row_vector<lower = 0>[nObsPK] cHatPKObsPred;
  vector<lower = 0>[nObsPD] cHatPDObsCond;
  row_vector<lower = 0>[nObsPD] cHatPDObsPred;
  
  // sample prior
  real prior_CL = lognormal_rng(log(5), 0.2);
  // real prior_Q = lognormal_rng(log(5), 0.2);
  real prior_V1 = lognormal_rng(log(50), 0.2);
  // real prior_V2 = lognormal_rng(log(100), 0.2);
  real prior_mtt = lognormal_rng(log(100), 0.2);
  real prior_circ0 = lognormal_rng(log(5), 0.2);
  real prior_alpha = lognormal_rng(log(0.1), 0.2);
  real prior_gamma = lognormal_rng(log(0.2), 0.2);
  
  xPred = pmx_solve_rk45(
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

  cHatPKPred = x[2, ] ./ V1;
  cHatPDPred = x[8, ] + theta[7];
  
  for(i in 1:nObsPK) {
    cHatPKObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), cHatPKObs[i])), sigmaPK));
    cHatPKObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), cHatPKObsPred[i])), sigmaPK));
  }
  
  for(i in 1:nObsPD) {
    cHatPDObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), cHatPDObs[i])), sigmaPD));
    cHatPDObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), cHatPDObsPred[i])), sigmaPD));
  }
}