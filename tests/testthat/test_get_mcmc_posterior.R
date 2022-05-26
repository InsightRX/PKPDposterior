# Set up environment for tests; re-use objects to reduce test time
mod <- load_model("pk_vanco_thomson")

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
  ruv = list(prop = 0.15, add = 1.6),
  ltbs = FALSE
)


test_that("Unsupported MCMC method produces error", {
  expect_error(get_mcmc_posterior(mod,data, method = "imaginary_method"))
})

test_that("Invalid STAN model produces error", {
  expect_error(get_mcmc_posterior(list(), data))
})

test_that("Get posterior estimate: hmc method", {
  # Underlying functions from dependencies use `print` which doesn't play nicely
  # with `expect_message` and similar functions. Use capture.output to catch.

  tmp <- capture.output(
    post1 <- get_mcmc_posterior(
      mod = mod,
      data = data,
      method = "hmc",
      iter_warmup = 5,  # not realistic value, just to speed up testing
      iter_sampling = 5, # not realistic value, just to speed up testing
      adapt_delta = 0.95,
      verbose = FALSE
    ),
    post2 <- get_mcmc_posterior(
      mod = mod,
      data = data,
      method = "vi",
      verbose = FALSE
    )
  )
  expected_names <- c(
    "raw", "draws_df", "settings", "data", "map", 
    "observed_post", "sampler_diagnostics"
  )
  
  expect_named(post1, expected_names, ignore.order = TRUE)
  expect_true(inherits(post1, "PKPDposterior"))
  expect_equal(
    post1$map,
    list(
      CL = 2.21950367643213, 
      Q = 2.60915639369181, 
      V1 = 0.577034739551428, 
      V2 = 0.43590902324682,
      TH_CRCL = 0.0153410247780514, 
      TDM_INIT = -0.00382955986679012
    )
  )
  expect_equal(post1$observed_post$median, c(30.1667, 11.3313))
  expect_true(inherits(post1$draws_df, "data.frame"))
  expect_true(inherits(post1$draws_df, "draws_df"))
  
  expect_named(post2, expected_names[1:6], ignore.order = TRUE)
  expect_true(inherits(post2, "PKPDposterior"))
  expect_equal(
    post2$map,
    list(
      lp_approx__ = -1.16811510661346, 
      CL = 2.68096798893404, 
      Q = 1.57646993508194, 
      V1 = 0.564520072805394, 
      V2 = 0.205024981316967,
      TH_CRCL = 0.0153795608145162, 
      TDM_INIT = -0.00132721983660977
    )
  )
  expect_equal(post2$observed_post$mean, c(30.2991531, 11.60416205))
  expect_true(inherits(post2$draws_df, "data.frame"))
  expect_true(inherits(post2$draws_df, "draws_df"))
})

test_that("Skip processing = TRUE returns unprocessed object", {
  # Underlying functions from dependencies use `print` which doesn't play nicely
  # with `expect_message` and similar functions. Use capture.output to catch.
  
  tmp <- capture.output(
    post2 <- get_mcmc_posterior(
      mod = mod,
      data = data,
      skip_processing = TRUE,
      method = "hmc",
      iter_warmup = 5,  # not realistic value, just to speed up testing
      iter_sampling = 5, # not realistic value, just to speed up testing
      adapt_delta = 0.95,
      verbose = FALSE
    )
  )
  expected_names <- c("raw", "draws_df", "settings", "data")
  
  expect_named(post2, expected_names, ignore.order = TRUE)
  expect_false(inherits(post2, "PKPDposterior"))
  expect_true(inherits(post2$draws_df, "data.frame"))
  expect_true(inherits(post2$draws_df, "draws_df"))
})
