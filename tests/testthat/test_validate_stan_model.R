test_that("Stan output and PKPDsim outputs are compared", {
  # set up inputs
  regimen <- PKPDsim::new_regimen(
    amt = 1500, 
    n = 4, 
    times = c(0, 12, 24, 36), 
    type = 'infusion',
    t_inf = 2
  )
  covariates <- list(
    WT = PKPDsim::new_covariate(value = 70, unit = "kg"),
    CRCL = PKPDsim::new_covariate(value = 5, unit = "l/hr"),
    CL_HEMO = PKPDsim::new_covariate(value = 0, unit = "l/hr")
  )
  tdm_data <- data.frame(
    t = c(2.5, 11.5), 
    dv = c(40, 14)
  )
  data <- new_stan_data(
    regimen,
    covariates, 
    tdm_data,
    dose_cmt = 2,
    parameters = list(
      CL = 2.99, 
      TH_CRCL = 0.0154, 
      Q = 2.28, 
      V2 = 0.732, 
      TDM_INIT = 0, 
      V1 = 0.675
    ),
    iiv = list(
      CL = 0.27, Q = 0.49, V1 = 0.15, V2 = 1.3, TH_CRCL = 0, TDM_INIT = 0
    ),
    ruv = list(prop = 0.15,add = 1.6)
  )
  
  # mock simulation outputs
  stan_output <- readRDS(test_path("data", "posterior_vanco_thomson_2.rds"))
  preds <- c(
    29.8335132689186, 9.94635377665206, 31.4276459500772, 11.7479800220603, 
    31.9493018602996, 12.758149364427, 24.4845567701294, 14.0371333175407,
    29.5437110187346, 8.96539061793017, 33.3060066029347, 14.0423845504041,
    30.4192041359278, 12.6715619628023, 34.4128280770078, 10.6007329934919,
    30.6510336756912, 12.011021375196, 36.6609207928251, 13.5621855882144,
    29.4546564175147, 11.282732387831, 32.1560512398151, 14.5954867052115,
    29.2213144250503, 11.5761673313184, 23.8614348913975, 16.0298889917733,
    32.3847310987911, 9.63537408249287, 27.7897754172766, 12.3518951580974,
    34.4516563449611, 11.3931831855777, 29.2311995888074, 7.59511626545619,
    26.7158836838708, 14.4587081171436, 31.2396517906757, 12.7226407220033,
    30.3030666116049, 12.8154881625673, 29.160007489036, 12.560474530974,
    32.6484503211349, 10.1530530585419, 31.268808108099, 10.732647949853,
    29.0043049225806, 9.17923076127185, 30.5585470409264, 11.5102036327435,
    33.5821299326876, 13.9811072249464, 25.6870394486081, 11.2761670211494,
    22.1823522660096, 11.6544954592759, 29.5784297926934, 16.5207723849649,
    29.5746466575885, 11.1607142145561, 29.5746466575885, 11.1607142145561,
    33.4486163870817, 14.9604334348592, 26.5512741571895, 12.2943185020551,
    33.3636966016467, 11.9615305267572, 26.6055465685213, 15.3448957563464,
    33.1759814157704, 19.6204190532324, 33.0569341486395, 13.385769623477,
    29.8532188831928, 10.604976791072, 32.0101357658655, 13.1206085242913,
    27.8441408589727, 10.1969991422083, 33.9463477883596, 14.9689800445401,
    25.8012322565308, 11.8959466574679, 28.8726843263567, 14.4398820724473,
    31.3798533832855, 9.6180936269985, 28.5539102213508, 15.5238742899571,
    32.6236667531091, 10.9621952651619, 35.4026959587258, 15.8778580527398,
    26.4624387414905, 11.481483470146, 29.2387486249668, 8.27583935050552
  )
  pkpdsim_output <- data.frame(
    id = 1,
    t = rep(c(2.5, 11.5), 50),
    comp = "obs",
    y = preds,
    obs_type = 1
  )
  
  # SIMULATIONS PASS COMPARISON
  mockery::stub(
    validate_stan_model,
    "get_mcmc_posterior",
    stan_output
  )
  mockery::stub(
    validate_stan_model,
    "purrr::map_dfr",
    pkpdsim_output
  )
  # check output for expected values: test passing
  message_outputs <- capture.output(
    validate_stan_model(
      stan_model = mod,
      pkpdsim_model = pkvancothomson::model(),
      data = data,
      mapping = list("V1" = "V")
    ),
    type = "message"
  )
  expect_true(any(grepl("Validation successful.", message_outputs)))
  
  
  # SIMULATIONS FAIL COMPARISON
  pkpdsim_output$y[1:3] <- 100
  mockery::stub(
    validate_stan_model,
    "purrr::map_dfr",
    pkpdsim_output
  )
  # check output for expected values: test failing
  expect_warning(
    message_fail <- capture.output(
        validate_stan_model(
        stan_model = mod,
        pkpdsim_model = pkvancothomson::model(),
        data = data,
        mapping = list("V1" = "V")
      ),
      type = "message"
    ),
    "Validation not succesful"
  )
  expect_true(any(grepl("70.1665", message_fail)))
  
})