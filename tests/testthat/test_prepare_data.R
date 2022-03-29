test_that("prepare_data returns correct format", {
  regimen <- PKPDsim::new_regimen(
    amt = 1500,
    n = 4,
    times = c(0, 12, 24, 36),
    type = 'infusion'
  )
  covariates <- list(
    WT = PKPDsim::new_covariate(
      value = c(150, 149.5),
      times = c(0, 30),
      unit = "kg"
    ),
    CRCL = PKPDsim::new_covariate(
      value = c(6.5, 6.7),
      times = c(0, 12),
      unit = "l/hr"
    )
  )
  tdm_data <- data.frame(
    t = c(1, 2),
    dv = c(900, 800),
    cmt = c(2, 2)
  )
  res <- prepare_data(
    regimen, 
    covariates, 
    tdm_data,
    parameters = list(CL = 5, V = 50),
    iiv = list(CL = 0.1, V = 0.2),
    ruv = list(prop = 0.1, add = 1)
  )
  
  ## Structure looks as expected
  expect_type(res, "list")
  expect_equal(res$n_t, 9)
  expect_equal(res$n_obs_pk, 2)

  ## Covariates are present and filled in
  expect_true(all(names(covariates) %in% names(res)))
  expect_equal(res$WT, c(rep(150, 7), rep(149.5, 2)))
  expect_equal(res$CRCL, c(rep(6.5, 5), rep(6.7, 4)))

  ## NONMEM columns are present
  nm_cols <- c("cmt", "evid", "addl", "ss", "amt", "time", "rate", "ii")
  expect_true(all(nm_cols %in% names(res)))
})

test_that("covariates_to_nm produces expected result", {
  covariates <- list(
    AGE = PKPDsim::new_covariate(value = 40),
    WT = PKPDsim::new_covariate(value = 70, unit = "kg"),
    CRCL = PKPDsim::new_covariate(
      value = c(6.5, 6.4, 6.6),
      times = c(0, 24, 48),
      unit = "l/hr"
    )
  )
  res <- covariates_to_nm(covariates)

  expect_named(
    res,
    c("TIME", "AGE", "WT", "CRCL", "ID", "EVID", "MDV", "AMT", "RATE", "DV", "CMT")
  )
  expect_equal(res$TIME, c(0, 24, 48))
})

test_that("tdm_to_nm produces expected result", {
  tdm_data <- data.frame(
    t = c(1, 2),
    dv = c(900, 800),
    cmt = c(2, 2)
  )
  res <- tdm_to_nm(tdm_data)

  expect_named(res, c("TIME", "DV", "CMT", "EVID", "ID", "MDV", "AMT", "RATE"))
  expect_equal(res$TIME, tdm_data$t)
  expect_equal(res$DV, tdm_data$dv)
  expect_equal(res$CMT, tdm_data$cmt)
})
