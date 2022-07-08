nm <- data.frame(
  ID = 1,
  OCC = 1,
  TIME = c(0, 12.2, 25.2, 37.4, 42.9, 49.7, 51.4, 70.1, 78.7, 91.8, 102.25),
  EVID = c(1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1),
  DV = c(0, 0, 0, 0, 28.8, 20.2, 0, 0, 0, 0, 0),
  AMT = c(1000, 1000, 1000, 1000, 0, 0, 1000, 1000, 1000, 1000, 1000),
  LENGTH = c(2, 2, 2, 2, 0, 0, 2, 1.5, 1.5, 1.5, 1.5),
  RATE = c(500, 500, 500, 500, 0, 0, 500, 666.7, 666.7, 666.7, 666.7),
  SEX = 0,
  AGE = 62,
  WT = 138,
  CR = c(1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.02, 1.1, 1.4),
  HT = 200,
  CRCL = c(8.8, 8.8, 8.8, 8.8, 8.8, 8.8, 8.8, 8.8, 8.8, 8.15, 6.404),
  IBW = 93.1,
  LBW = 83.6,
  FFM = 90.5
)


test_that("Missing columns are caught", {
  expect_error(
    nonmem_to_stan_data(
      nm,
      c("WT", "required_covariate"),
      parameters = list(CL = 5, V = 50),
      iiv = list(CL = 0.1, V = 0.2),
      ruv = list(prop = 0.1, add = 1)
    ),
    "all(covariate_cols %in% colnames(nm)) is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    nonmem_to_stan_data(
      nm[,setdiff(colnames(nm), "EVID")],
      c("WT", "CRCL"),
      parameters = list(CL = 5, V = 50),
      iiv = list(CL = 0.1, V = 0.2),
      ruv = list(prop = 0.1, add = 1)
    ),
    "all(c(\"EVID\", \"TIME\", \"DV\") %in% colnames(nm)) is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    nonmem_to_stan_data(
      nm[,setdiff(colnames(nm), "ID")],
      c("WT", "CRCL"),
      parameters = list(CL = 5, V = 50),
      iiv = list(CL = 0.1, V = 0.2),
      ruv = list(prop = 0.1, add = 1)
    ),
    NA
  )
})

test_that("More than one patient throws error", {
  nm$ID[1] <- 2
  expect_error(
    nonmem_to_stan_data(
      nm,
      c("WT", "required_covariate"),
      parameters = list(CL = 5, V = 50),
      iiv = list(CL = 0.1, V = 0.2),
      ruv = list(prop = 0.1, add = 1)
    ),
    "Please supply NONMEM data for one patient at a time.",
    fixed = TRUE
  )
})

test_that("Missing TDMs or doses throws correct errors", {
  expect_error(
    nonmem_to_stan_data(
      nm[nm$EVID == 0, ],
      c("WT"),
      parameters = list(CL = 5, V = 50),
      iiv = list(CL = 0.1, V = 0.2),
      ruv = list(prop = 0.1, add = 1)
    ),
    "any(nm$EVID == 1) is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    nonmem_to_stan_data(
      nm[nm$EVID == 1, ],
      c("WT"),
      parameters = list(CL = 5, V = 50),
      iiv = list(CL = 0.1, V = 0.2),
      ruv = list(prop = 0.1, add = 1)
    ),
    "any(nm$EVID == 0) is not TRUE",
    fixed = TRUE
  )
})

test_that("returns PKPDposterior data object", {
  obj <- nonmem_to_stan_data(
    nm,
    c("WT", "CRCL"),
    parameters = list(CL = 5, V = 50),
    iiv = list(CL = 0.1, V = 0.2),
    ruv = list(prop = 0.1, add = 1)
  )
  # Overall structure is correct
  expect_true(inherits(obj, "list"))
  expect_named(
    obj,
    c("parameters", "fixed", "regimen", "covariates", "data", "iiv", "ruv", "stan_data"),
    ignore.order = TRUE
  )
  
  # regimen and covariate objects correctly parsed
  expect_true(inherits(obj$regimen, "regimen"))
  expect_true(inherits(obj$covariates[["WT"]], "covariate"))
  expect_true(inherits(obj$covariates[["CRCL"]], "covariate"))
  expect_equal(obj$covariates[["CRCL"]]$value, nm$CRCL)
  
  # observations object correctly parsed
  expect_equal(obj$data$t, c(42.9, 49.7))
  expect_equal(obj$data$dv, c(28.8, 20.2))
  expect_equal(obj$data$cmt, c(2, 2))
})

test_that("'No covariates' does not produce error", {
  expect_error(
    nonmem_to_stan_data(
      nm,
      NULL,
      parameters = list(CL = 5, V = 50),
      iiv = list(CL = 0.1, V = 0.2),
      ruv = list(prop = 0.1, add = 1)
    ),
    NA
  )
  expect_error(
    nonmem_to_stan_data(
      nm,
      c(),
      parameters = list(CL = 5, V = 50),
      iiv = list(CL = 0.1, V = 0.2),
      ruv = list(prop = 0.1, add = 1)
    ),
    NA
  )
})

test_that("Multiple obs cmts handled", {
  nm$CMT <- 2
  nm$CMT[nm$EVID == 0] <- c(1, 2)
  obj <- nonmem_to_stan_data(
    nm,
    c("WT", "FFM"),
    parameters = list(CL = 5, V = 50),
    iiv = list(CL = 0.1, V = 0.2),
    ruv = list(prop = 0.1, add = 1)
  )
  expect_equal(obj$data$cmt, c(1, 2))
  expect_equal(obj$stan_data$cmt[obj$stan_data$DV != 0], c(1, 2))
})

test_that("covariate interpolation configurable", {
  obj <- nonmem_to_stan_data(
    nm,
    c("CR", "CRCL"),
    parameters = list(CL = 5, V = 50),
    iiv = list(CL = 0.1, V = 0.2),
    ruv = list(prop = 0.1, add = 1),
    covariates_implementation = list(
      CR = "interpolate",
      CRCL = "locf"
    )
  )
  expect_equal(obj$covariates[["CRCL"]]$implementation, "locf")
  expect_equal(obj$covariates[["CR"]]$implementation, "interpolate")
})

