functions{
  vector ode_rhs(real t, vector x, real[] parms, real[] x_r, int[] x_i){
    real CL = parms[1];
    real Q = parms[2];
    real V1 = parms[3];
    real V2 = parms[4];
    real KA = parms[5];
    real KM = parms[6];
    real VMAX1 = parms[7];
    
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    
    real VMAXINH = exp(1.5) / (1 + exp(1.5));
    
    real Vmax = VMAX1 * (1 - VMAXINH * (t-1) / ((t-1) + (2.41 - 1))); // T50 = 2.41
    
    vector[3] y;

    y[1] = -KA*x[1];
    y[2] = KA*x[1] - (k10 + k12)*x[2] + k21*x[3] - (Vmax * x[2]/V1)/(KM + x[2]/V1);
    y[3] = k12*x[2] - k21*x[3];

    return y;
  }
}

data {
  int<lower = 1> nt;            // number of events
  int<lower = 1> nObs;          // number of observations
  int<lower = 1> iObs[nObs];    // index of observation
  
  // NONMEM data
  int<lower = 1> cmt[nt];
  int evid[nt];
  int addl[nt];
  int ss[nt];
  real amt[nt];
  real time[nt];
  real rate[nt];
  real ii[nt];
  
  // covariates
  real WT[nt];
  
  row_vector<lower = 0>[nObs] cObs;  // observed concentration (dependent variable)
}

transformed data {
  row_vector[nObs] logCObs = log(cObs);
  int nTheta = 7;   // number of parameters
  int nCmt = 3;   // number of compartments
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
  real theta[nTheta];
  row_vector<lower = 0>[nt] cHat;
  row_vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[3, nt] x; 

  theta[1] = CL * pow(mean(WT)/70, 0.75);
  theta[2] = Q * 1.637 * pow(mean(WT)/70, 0.75);
  theta[3] = V1 * mean(WT)/70;
  theta[4] = V2 * mean(WT)/70;
  theta[5] = KA;
  theta[6] = KM;
  theta[7] = VMAX1 * pow(mean(WT)/70, 0.75);
  // theta[8] = F1;
  
  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  cHat = x[2, ] ./ V1;

  for(i in 1:nObs){
    cHatObs[i] = cHat[iObs[i]];  // predictions for observed data records
  }
}

model{
  // informative prior
  CL ~ lognormal(log(6.16), 0.5);
  Q ~ lognormal(log(15.5), 0.42);
  V1 ~ lognormal(log(79), 0.136);
  V2 ~ lognormal(log(103), 0.77);
  KA ~ lognormal(log(1.19), 0.9);
  KM ~ lognormal(log(1.15), 1.0);
  VMAX1 ~ lognormal(log(114), 0.50);
  // F1 ~ lognormal(log(0.585), 1.0);
  
  cObs ~ normal(cHatObs, (0.3 * cHatObs + 0.01));
}

generated quantities{
  real cObsPred[nObs];

  // sample prior
  real prior_CL = lognormal_rng(log(6.16), 0.50); //"parameters", "iiv" values
  real prior_Q = lognormal_rng(log(15.5), 0.42);
  real prior_V1 = lognormal_rng(log(79), 0.136);
  real prior_V2 = lognormal_rng(log(103), 0.77);
  real prior_KA = lognormal_rng(log(1.19), 0.90);
  real prior_KM = lognormal_rng(log(1.15), 1.0);
  real prior_VMAX1 = lognormal_rng(log(114), 0.50);
  // real prior_F1 = lognormal_rng(log(0.585), 1.0);

  for(i in 1:nObs){
    cObsPred[i] = normal_rng(cHatObs[i], (0.3 * cHatObs[i] + 0.01));
  }
}