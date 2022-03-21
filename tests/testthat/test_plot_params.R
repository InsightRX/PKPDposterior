post <- readRDS(test_path("data", "posterior_vanco_thomson_1.rds"))

test_that("plot_params parses posterior object and returns ggplot2 object", {
  p <- plot_params(post)
  expect_equal(class(p), c("gg", "ggplot"))
})

test_that("plot_params fails when no draws available", {
  post$draws_df <- NULL
  expect_error(
    plot_params(post),
    "Provided posterior object does not contain expected info."
  )
})