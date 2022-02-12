// TwoCptModel.stan
// Run two compartment model using built-in analytical solution 
// Heavily anotated to help new users

// using proportional + additive error model!!

data{
  int<lower = 1> nt;  // number of events
  int<lower = 1> nObs;  // number of observation
  int<lower = 1> iObs[nObs];  // index of observation
  
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
  
  vector<lower = 0>[nObs] cObs;  // observed concentration (Dependent Variable)
}

transformed data{
  vector[nObs] logCObs = cObs;
  int nTheta = 5;  // number of ODE parameters in Two Compartment Model
  int nCmt = 3;  // number of compartments in model
}

parameters{
  real<lower = 0> CL;
  real<lower = 0> Q;
  real<lower = 0> V1;
  real<lower = 0> V2;
}

transformed parameters{
  array[nt, nTheta] real theta;  // ODE parameters
  row_vector<lower = 0>[nt] cHat;
  vector<lower = 0>[nObs] cHatObs;
  matrix<lower = 0>[nCmt, nt] x;
  
  for(j in 1:nt) {
    theta[j, 1] = CL * (1.0 + 0.0154 * ((CRCL[j] * 16.6667) - 66.0));
    theta[j, 2] = Q;
    theta[j, 3] = V1 * WT[j];
    theta[j, 4] = V2 * WT[j];
    theta[j, 5] = 0; //ka = 0, IV model
  }

  x = pmx_solve_twocpt(time, amt, rate, ii, evid, cmt, addl, ss, theta);
  
  cHat = x[2, :] ./ (V1 * mean(WT)); // we're interested in the amount in the second compartment

  cHatObs = cHat'[iObs]; // predictions for observed data recors
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
  
  for(i in 1:nObs){
    cObsPred[i] = cHatObs[i];
    cObsPred[i] += normal_rng(0, 0.15);
    cObsPred[i] += normal_rng(0, 0.1);
  }
}
