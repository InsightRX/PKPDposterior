model_pkpd <- readLines(test_path("data", "neutropenia_model_pkpd.stan"))
model_pd   <- readLines(test_path("data", "neutropenia_model_pd.stan"))
data_pkpd  <- readRDS(test_path("data", "neutropenia_data_pkpd.rds"))
data_pd  <- readRDS(test_path("data", "neutropenia_data_pd.rds"))

test_that("throws no error when all obs types present", {
  expect_silent(
    PKPDposterior:::check_input_data(
      model_pkpd,
      data_pkpd
    )
  )
})


test_that("check_input_data() throws error when not all observation types present", {
  expect_error(
    PKPDposterior:::check_input_data(
      model_pkpd,
      data_pd
    )
  )
})


test_that("check_input_data() throws message when not all observation types are defined", {
  expect_message(
    PKPDposterior:::check_input_data(
      model_pd,
      data_pkpd
    )
  )
})

