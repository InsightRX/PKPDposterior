pred_post <- readRDS(test_path("data", "pred_post_vanco_thomson.rds"))
  
test_that("plot_predictions without TDM data returns ggplot2 object", {
  p <- plot_predictions(pred_post)
  expect_equal(class(p), c("gg", "ggplot"))
  expect_equal(length(p$layers), 2)
})

test_that("plot_predictions *with* TDM data returns ggplot2 object", {
  p <- plot_predictions(
    pred_post, 
    obs = data.frame(t = c(2, 11), c(30, 12))
  )
  expect_equal(class(p), c("gg", "ggplot"))
  expect_equal(length(p$layers), 4)
})

test_that("plot_predictions fails when incorrect data object supplied", {
  pred_post$t <- NULL
  expect_error(
    plot_predictions(pred_post),
    regexp = "Requires a simulation dataset created using sim_from_draws()"
  )
})