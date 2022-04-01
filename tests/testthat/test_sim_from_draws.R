post <- readRDS(test_path("data", "posterior_vanco_thomson_1.rds"))
model <- PKPDsim::new_ode_model(
  code = "dAdt[0] = -(CL/V)*A[0] - (Q/V)*A[0] + (Q/V2)*A[1]; dAdt[1] = (Q/V)*A[0] - (Q/V2)*A[1];",
  parameters = list(CL = 1, Q =1 , V=1, V2=1),
  obs = list(scale="V", cmt = 1), 
  install = FALSE
)

test_that("sim_from_draws correctly simulates from draws in posterior object", {
  res <- sim_from_draws(
    post,
    model,
    map = list("V1" = "V"),
    regimen = PKPDsim::new_regimen(
      amt = 1000,
      n = 4, 
      interval = 12,
      t_inf = 1,
      type = "infusion"
    ),
    seed = 12345,
    prior = FALSE,
    n = NULL,
    summarize = FALSE
  )
  expect_equal(class(res), c("PKPDsim_data", "data.frame"))
  expect_equal(nrow(res), 75000)
  expect_true(!any(is.na(res$y)))
})


test_that("sim_from_draws draws less patients when asked", {
  res <- sim_from_draws(
    post,
    model,
    map = list("V1" = "V"),
    regimen = PKPDsim::new_regimen(
      amt = 1000,
      n = 4, 
      interval = 12,
      t_inf = 1,
      type = "infusion"
    ),
    seed = 12345,
    prior = FALSE,
    n = 10,
    summarize = FALSE
  )
  expect_equal(class(res), c("PKPDsim_data", "data.frame"))
  expect_equal(nrow(res), 1500)
  expect_equal(unique(res$id), 1:10)
})

test_that("sim_from_draws draws from prior instead of posterior when asked", {
  res <- sim_from_draws(
    post,
    model,
    map = list("V1" = "V"),
    regimen = PKPDsim::new_regimen(
      amt = 1000,
      n = 4, 
      interval = 12,
      t_inf = 1,
      type = "infusion"
    ),
    seed = 12345,
    prior = TRUE,
    n = 10,
    summarize = FALSE
  )
  expect_equal(class(res), c("PKPDsim_data", "data.frame"))
  expect_equal(round(res$y[1:5],1), c(0.0, 211.3, 31.7, 8.9, 2.5))
})

test_that("sim_from_draws errors when posterior object is malformed", {
  post$draws_df <- NULL
  expect_error(
    res <- sim_from_draws(
      post,
      model,
      map = list("V1" = "V"),
      regimen = PKPDsim::new_regimen(
        amt = 1000,
        n = 4, 
        interval = 12,
        t_inf = 1,
        type = "infusion"
      ),
      seed = 12345,
      prior = TRUE,
      n = 10,
      summarize = FALSE
    ),
    "Provided posterior object does not contain expected info."
  )
})