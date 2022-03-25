functions{
  vector ode_rhs(real t, vector x, real[] parms, real[] x_r, int[] x_i){
    real CL = parms[1];
    real Q = parms[2];
    real V1 = parms[3];
    real V2 = parms[4];
    real ka = 0;
    
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    
    vector[3] y;

    y[1] = -ka*x[1];
    y[2] =  ka*x[1] - (k10 + k12)*x[2] + k21*x[3];
    y[3] =  k12*x[2] - k21*x[3];

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
  real WT[nt];
  real CRCL[nt];
  
  row_vector<lower = 0>[nObs] cObs;  // observed concentration (dependent variable)
}

transformed data {
  row_vector[nObs] logCObs = log(cObs);
  int nTheta = 5;   // number of parameters
  int nCmt = 3;   // number of compartments
}

parameters {
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
}

transformed parameters {
  real theta[nTheta];
  row_vector<lower = 0>[nt] cHat;
  row_vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[3, nt] x; 

  theta[1] = CL * (1.0 + 0.0154 * ((mean(CRCL) * 16.6667) - 66.0));;
  theta[2] = Q;
  theta[3] = V1 * mean(WT);
  theta[4] = V2 * mean(WT);
  theta[5] = 0;

  x = pmx_solve_rk45(ode_rhs, 3, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  cHat = x[2, ] ./ (V1 * mean(WT));

  for(i in 1:nObs){
    cHatObs[i] = cHat[iObs[i]];  // predictions for observed data records
  }
}

model{
  // informative prior
  CL ~ lognormal(log(2.99), 0.27);
  Q ~ lognormal(log(2.28), 0.49);
  V1 ~ lognormal(log(0.675), 0.15);
  V2 ~ lognormal(log(0.732), 1.3);
  cObs ~ normal(cHatObs, (0.15 * cHatObs + 1.6));
}

generated quantities{
  real cObsPred[nObs];
  
  // sample prior
  real prior_CL = lognormal_rng(log(2.99), 0.27);
  real prior_Q = lognormal_rng(log(2.28), 0.49);
  real prior_V1 = lognormal_rng(log(0.675), 0.15);
  real prior_V2 = lognormal_rng(log(0.732), 1.3);

  for(i in 1:nObs){
    cObsPred[i] = normal_rng(cHatObs[i], (0.15 * cHatObs[i] + 1.6));
  }
}