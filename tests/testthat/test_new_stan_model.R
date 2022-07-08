test_that("Neutropenia model correctly generates", {
  parameters <- list(CL = 5, V = 50,
                     SLOPE = 0.1, MTT = 100, 
                     CIRC0 = 5, GAMMA = 0.2)
  iiv <- list(
    CL = 0.2, V = 0.5, 
    SLOPE = 1.0, MTT = 0.2, 
    CIRC0 = 0.5, GAMMA = 0.2
  )
  ruv <- list(
    pk = list(add = 0.2),
    pd = list(add = 0.3)
  )
  parameter_definitions <- list(
    "CL" = "CL",
    "Q" = 0,
    "V" = "V",
    "V2" = 1,
    "KA" = 0,
    "MTT" = "MTT",
    "CIRC0" = "CIRC0",
    "GAMMA" = "GAMMA",
    "ALPHA" = "SLOPE"
  )
  ode <- "
  real ka  = 0;
  real k10 = CL / V;
  real k12 = 0;
  real k21 = 0;
  real ktr = 4 / MTT;
  
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
  
  EDrug = ALPHA * conc; // slope model, not Emax
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
  "
  code <- new_stan_model(
    parameters = parameters,
    parameter_definitions = parameter_definitions,
    ode = ode,
    covariate_definitions = NULL,
    obs_types = c("pk", "pd"),
    custom_ipred = list(
      "pk" = "A[2, ] ./ V;",
      "pd" = "A[8, ] + theta[7];"
    ),
    solver = "pmx_solve_rk45",
    verbose = F
  )
  code_ref <- readRDS(test_path("data", "neutropenia_model_code.rds"))
  # saveRDS(code, test_path("data", "neutropenia_model_code.rds"))
  expect_equal(code, code_ref)
})

test_that("variable_definitions argument works", {
  parameters <- list(
    CL = 3.96,
    V = 10.8
  )
  iiv <- list(CL = 0.24, V = 0.17)
  ruv <- list(prop = 0.106, add = 22.2)
  model <- new_stan_model(
    parameters = parameters,
    variable_definitions = list(
      "BMI" = "WT[j]/(HT[j]*HT[j]/10000)",
      "FFM" = "(SEX[j] * (0.88 + ((1-0.88)/(1+(AGE[j]/13.4)^-12.7))) * ((9270 * WT[j])/(6680 + (216 * BMI)))) + ((1-SEX[j])*(1.11 + ((1-1.11)/(1+(AGE[j]/7.1)^-1.1))) * ((9270 * WT[j])/(8780 + (244 * BMI))))"
    ),
    parameter_definitions = list(
      "CL" = "CL * (FFM/12)^0.75 * (0.451 + (1-0.451)*(1-exp(-AGE[j]*1.37))) * (1 + (-0.2 * REGI[j])) * (1 -0.135*(time[j]>24))",
      "V"  = "V * FFM/12",
      "KA" = "0" 
    ),
    covariate_definitions = list(
      "SEX" = "int",
      "WT" = "real",
      "HT" = "real",
      "AGE" = "real",
      "REGI" = "int"
    ),
    solver = "pmx_solve_onecpt",
    scale = "((V * FFM/12)/1000)",
    verbose = T
  )
  expect_true(length(grep("BMI = ", model)) == 1)
  expect_true(length(grep("FFM = ", model)) == 1)
})
