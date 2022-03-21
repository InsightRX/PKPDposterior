test_that("Output from get_parameter_tables is OK", {
  post <- readRDS(test_path("data", "posterior_vanco_thomson_1.rds"))
  res <- get_parameter_tables(post)
  expect_equal(names(res), c("posterior", "prior"))
  expect_equal(nrow(res$posterior), 2000)
  expect_equal(unique(res$posterior$name), c("CL", "Q", "V2", "V1"))
})

test_that("Output from get_parameter_tables when no prior parameter data should still parse posterior data", {
  post <- readRDS(test_path("data", "posterior_vanco_thomson_1.rds"))
  post$draws_df$prior_CL <- NULL
  post$draws_df$prior_Q <- NULL
  post$draws_df$prior_V1 <- NULL
  post$draws_df$prior_V2 <- NULL
  res <- get_parameter_tables(post)
  expect_equal(names(res), c("posterior"))
  expect_equal(nrow(res$posterior), 2000)
  expect_equal(unique(res$posterior$name), c("CL", "Q", "V2", "V1"))
})

test_that("get_parameter_tables with malformed post object stops", {
  post <- readRDS(test_path("data", "posterior_vanco_thomson_1.rds"))
  post$draws_df <- NULL
  expect_error(get_parameter_tables(post),  "Provided posterior object does not contain expected info.")
})
