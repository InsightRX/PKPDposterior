//
// This Stan program defines the 2cmt IV vanco model
//originally published by Thomson et al.
//

functions {
  vector TwoCptIV_Model(real t, vector cmt, real[] parms, real[] x_r, int[] x_i){
    
    real CL  = parms[1];
    real V1  = parms[2];
    real Q   = parms[3];
    real V2  = parms[4];
    
    vector[2] dxdt_cmt;

    dxdt_cmt[1] = (Q/V2) * cmt[2] - (CL + Q) * cmt[1] / V1;
    dxdt_cmt[2] = (Q/V1) * cmt[1] - (Q/V2) * cmt[2];

    return dxdt_cmt;
  }

  real cv_to_sd(real cv) {
    return sqrt(log(cv^2+1));
  }

  real sigma(real conc){
    return (conc * 0.15 + 1.6);  // proportional and residual error
  }
}


// The input data is a vector 'y' of length 'N'.
data{
  int<lower = 1> nt;  // number of events
  int<lower = 1> nObs;  // number of observation
  int<lower = 1> iObs[nObs];  // index of observation
  
  // NONMEM data, general columns
  int<lower = 1> cmt[nt];
  int evid[nt];
  int addl[nt];
  int ss[nt];
  real amt[nt];
  real time[nt];
  real rate[nt];
  real ii[nt];
  
  // Covariates specific to model
  real<lower=0> CRCL[nt];
  real<lower=0> WT[nt];
  
  // Model parameters
  real Prior_CL;
  real Prior_V1;
  real Prior_Q;
  real Prior_V2;

  real Prior_CL_omega;
  real Prior_V1_omega;
  real Prior_Q_omega;
  real Prior_V2_omega;
  
  vector<lower = 0>[nObs] cObs;  // observed concentration (Dependent Variable)
}

// The parameters accepted by the model. Our model
// accepts four parameters: CL, V1, V2, Q
parameters {
  real<lower=0.0> CL;
  real<lower=0.0> V1;
  real<lower=0.0> Q;
  real<lower=0.0> V2;
}

transformed parameters {
  real theta[4];

  row_vector[nt] cHat;
  row_vector[nObs] cHatObs;

  real<lower=0> sigmaEPS[nObs];
  real CLi;
  real Vi;
  real V2i;
  
  CLi = CL .* (1 .+ 0.0154 .* (CRCL .* 16.66667 .- 66));
  V1i = V1 .* WT;
  V2i = V2 .* WT;

  theta[1] = CLi;
  theta[2] = V1i;
  theta[3] = Q;
  theta[4] = V2i;

  matrix[nCmt, nt] x = pmx_solve_rk45(TwoCptIV_Model, nCmt, time, amt, rate, ii, evid, cmt, addl, ss, theta, 1e-5, 1e-8, 1e5);

  cHat = x[1, ] ./ V1;

  for (i in 1:nObs) {
    cHatObs[i] = cHat[iObs[i]];
    sigmaEPS[i] = sigma(cHatObs[i]);
  }
}

// The model to be estimated
model {
  CL ~ lognormal(log(2.99), 0.27);
  V1 ~ lognormal(log(0.675), 0.15);
  Q  ~ lognormal(log(2.28),  0.49);
  V2 ~ lognormal(log(0.732), 1.30);

  for ( i in 1:nObs){
    cObs[i] ~ normal(cHatObs[i], sigmaEPS[i]);
  }
}
