// PK-PD model for neutropenia based on Friberg et al.
// Drug is fictitious, just a simple 2-cmt oral model.
functions{
  
  real[] twoCptNeutModelODE(real t,
                            real[] x,
                            real[] parms,
                            real[] rdummy,
                            int[] idummy){
    real CL = parms[1];
    real Q = parms[2];
    real V1 = parms[3];
    real V2 = parms[4];
    real ka = params[5];
    real mtt = parms[6];
    real circ0 = parms[7];
    real gamma = parms[8];
    real alpha = parms[9];
    
    real k10 = CL / V1;
    real k12 = Q / V1;
    real k21 = Q / V2;
    real ktr = 4 / mtt;
    
    real dxdt[8];
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

data{
  int<lower = 1> nt;
  int<lower = 1> nObsPK;
  int<lower = 1> nObsPD;
  int<lower = 1> iObsPK[nObsPK];
  int<lower = 1> iObsPD[nObsPD];
  real<lower = 0> amt[nt];
  int cmt[nt];
  int<lower = 0> evid[nt];
  real<lower = 0> time[nt];
  real<lower = 0> ii[nt];
  int<lower = 0> addl[nt];
  int<lower = 0> ss[nt];
  real rate[nt];
  vector<lower = 0>[nObsPK] cObs;
  vector<lower = 0>[nObsPD] neutObs;

  // // data for population model
  // int<lower = 1> nSubjects;
  // int<lower = 1> start[nSubjects];
  // int<lower = 1> end[nSubjects];
  // real<lower = 0> weight[nSubjects];

  // data for priors
  real<lower = 0> CLHatPrior;
  real<lower = 0> QHatPrior;
  real<lower = 0> V1HatPrior;
  real<lower = 0> V2HatPrior;
  real<lower = 0> kaHatPrior;
  real<lower = 0> CLHatPriorCV;
  real<lower = 0> QHatPriorCV;
  real<lower = 0> V1HatPriorCV;
  real<lower = 0> V2HatPriorCV;
  real<lower = 0> kaHatPriorCV;
  real<lower = 0> circ0HatPrior;
  real<lower = 0> circ0HatPriorCV;
  real<lower = 0> mttHatPrior;
  real<lower = 0> mttHatPriorCV;
  real<lower = 0> gammaPrior;
  real<lower = 0> gammaPriorCV;
  real<lower = 0> alphaHatPrior;
  real<lower = 0> alphaHatPriorCV;
}

transformed data{
  vector[nObsPK] logCObs = log(cObs);
  vector[nObsPD] logNeutObs = log(neutObs);
  int nTheta = 9;  // number of ODE parameters
  int nIIV = 7;  // parameters with IIV
  int nCmt = 8;  // number of compartments
  real biovar[nCmt];
  real tlag[nCmt];
  
  for (i in 1:nCmt) {
    biovar[i] = 1;
    tlag[i] = 0;
  }
}

parameters{
  real<lower = 0> CLHat;
  real<lower = 0> QHat;
  real<lower = 0> V1Hat;
  real<lower = 0> V2Hat;
  real<lower = 0> kaHat;
  real<lower = 0> mttHat;
  real<lower = 0> circ0Hat;
  real<lower = 0> alphaHat;
  real<lower = 0> gamma;
  real<lower = 0> sigma;
  real<lower = 0> sigmaNeut;
  
  // IIV parameters
  cholesky_factor_corr[nIIV] L;
  vector<lower = 0>[nIIV] omega;
  // matrix[nIIV, nSubjects] etaStd;
  
}

transformed parameters{
  row_vector[nt] cHat;
  row_vector[nObsPK] cHatObs;
  row_vector[nt] neutHat;
  row_vector[nObsPD] neutHatObs;
  matrix[nCmt, nt] x;
  real<lower = 0> parms[nTheta]; // The [1] indicates the parameters are constant
  
  // variables for Matt's trick
  vector<lower = 0>[nIIV] theta;

  // Matt's trick to use unit scale
  thetaHat[1] = CLHat; 
  thetaHat[2] = QHat;
  thetaHat[3] = V1Hat;
  thetaHat[4] = V2Hat;
  thetaHat[5] = mttHat;
  thetaHat[6] = circ0Hat;
  thetaHat[7] = alphaHat;

  // for(i in 1:nSubjects) {

    parms[1] = theta[i, 1] * (weight[i] / 70)^0.75; // CL
    parms[2] = theta[i, 2] * (weight[i] / 70)^0.75; // Q
    parms[3] = theta[i, 3] * (weight[i] / 70); // V1
    parms[4] = theta[i, 4] * (weight[i] / 70); // V2
    parms[5] = kaHat; // ka
    parms[6] = theta[i, 5]; // mtt
    parms[7] = theta[i, 6]; // circ0
    parms[8] = gamma;
    parms[9] = theta[i, 7]; // alpha
  
    x[start[i]:end[i]] = pmx_solve_rk45(
      twoCptNeutModelODE, 
      nCmt,
      time[start[i]:end[i]], 
      amt[start[i]:end[i]], 
      rate[start[i]:end[i]], 
      ii[start[i]:end[i]], 
      evid[start[i]:end[i]], 
      cmt[start[i]:end[i]], 
      addl[start[i]:end[i]], 
      ss[start[i]:end[i]],
      parms, 
      biovar, 
      tlag,
      1e-6, 
      1e-6, 
      1e6
    );
    
    cHat[start[i]:end[i]] = x[2, start[i]:end[i]] / parms[3];  // divide by V1
    neutHat[start[i]:end[i]] = x[8, start[i]:end[i]] + parms[7];  // Add baseline
    
  }
  
  cHatObs = cHat[iObsPK];
  neutHatObs = neutHat[iObsPD];

}

model{
  // Priors
  CLHat ~ lognormal(log(CLHatPrior), CLHatPriorCV);
  QHat ~ lognormal(log(QHatPrior), QHatPriorCV);
  V1Hat ~ lognormal(log(V1HatPrior), V1HatPriorCV);
  V2Hat ~ lognormal(log(V2HatPrior), V2HatPriorCV);
  kaHat ~ lognormal(log(kaHatPrior), kaHatPriorCV);
  sigma ~ cauchy(0, 1);
  
  mttHat ~ lognormal(log(mttHatPrior), mttHatPriorCV);
  circ0Hat ~ lognormal(log(circ0HatPrior), circ0HatPriorCV);
  alphaHat ~ lognormal(log(alphaHatPrior), alphaHatPriorCV);
  gamma ~ lognormal(log(gammaPrior), gammaPriorCV);
  sigmaNeut ~ cauchy(0, 1);

  // observed data likelihood
  logCObs ~ normal(log(cHatObs), sigma);
  logNeutObs ~ normal(log(neutHatObs), sigmaNeut);
}

generated quantities{
  matrix[nCmt, nt] xPred;
  real<lower = 0> parmsPred[nTheta];
  row_vector[nt] cHatPred;
  row_vector[nt] neutHatPred;
  vector<lower = 0>[nObsPK] cHatObsCond;
  row_vector<lower = 0>[nObsPK] cHatObsPred;
  vector<lower = 0>[nObsPD] neutHatObsCond;
  row_vector<lower = 0>[nObsPD] neutHatObsPred;
  
  // Variables for IIV  
  matrix[nIIV, nSubjects] etaStdPred;
  matrix<lower = 0>[nSubjects, nIIV] thetaPredM;
  corr_matrix[nIIV] rho;
  
  // for(i in 1:nSubjects) {
    parmsPred[1] = theta[1] * (weight[i] / 70)^0.75; // CL
    parmsPred[2] = theta[2] * (weight[i] / 70)^0.75; // Q
    parmsPred[3] = theta[3] * (weight[i] / 70); // V1
    parmsPred[4] = theta[4] * (weight[i] / 70); // V2
    parmsPred[5] = kaHat; // ka
    parmsPred[6] = theta[i, 5]; // mtt
    parmsPred[7] = theta[i, 6]; // circ0
    parmsPred[8] = gamma; // gamma
    parmsPred[9] = theta[i, 7]; // alpha
    
    xPred[start[i]:end[i]] = pmx_solve_rk45(
      twoCptNeutModelODE, 
      nCmt,
      time[start[i]:end[i]], 
      amt[start[i]:end[i]],
      rate[start[i]:end[i]],
      ii[start[i]:end[i]],
      evid[start[i]:end[i]],
      cmt[start[i]:end[i]],
      addl[start[i]:end[i]],
      ss[start[i]:end[i]],
      parmsPred, 
      biovar, 
      tlag,
      1e-6, 1e-6, 1e6
    );
    
    cHatPred[start[i]:end[i]] = xPred[2, start[i]:end[i]] / parmsPred[3]; // divide by V1
    neutHatPred[start[i]:end[i]] = xPred[8, start[i]:end[i]] + parmsPred[7]; // Add baseline
  // }
  
  // predictions for observed data records
  cHatObsPred = cHatPred[iObsPK];
  neutHatObsPred = neutHatPred[iObsPD];
  
  for(i in 1:nObsPK) {
    cHatObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObs[i])), sigma));
    cHatObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), cHatObsPred[i])), sigma));
  }
  
  for(i in 1:nObsPD) {
    neutHatObsCond[i] = exp(normal_rng(log(fmax(machine_precision(), neutHatObs[i])), sigmaNeut));
    neutHatObsPred[i] = exp(normal_rng(log(fmax(machine_precision(), neutHatObsPred[i])), sigmaNeut));
  }
}