## Define regimen, covariates, and TDM data
regimen <- PKPDsim::new_regimen(
  amt = 1500, 
  n = 4, 
  times = c(0, 12, 24, 36), 
  type = 'infusion',
  t_inf = 2
)
covariates <- list(
  WT = PKPDsim::new_covariate(value = 70, unit = "kg"),
  CRCL = PKPDsim::new_covariate(value = 5, unit = "l/hr")
)
tdm_data <- data.frame(
  t = c(2.5, 11.5), 
  dv = c(40, 14)
)
prior <- list(CL = 2.99, Q = 2.28, V2 = 0.732, V1 = 0.675)

test_that("fixing model parameter works", {

  ## Create combined dataset for Torsten/Stan to read:
  expect_message(
    res <- new_stan_data(
      regimen,
      covariates, 
      tdm_data,
      dose_cmt = 2,
      parameters = prior,
      iiv = list(CL = 0.27, Q = 0.49, V1 = 0.15, V2 = 1.3),
      fix = "V2",
      ruv = list(
        prop = 0.15,
        add = 1.6
      )
    ),
    "Fixing V2 parameter. Please note that this slows down HMC sampling and should be used only for testing. For production models, remove the parameter from the likelihood definition and hardcode it in the model definition."
  )
  
  expect_equal( # fixes V2
    res$stan_data$omega_V2, 
    0.001
  )
  expect_equal( # but not V1
    res$stan_data$omega_V1, 
    0.15
  )
  
})