test_that("Expected print output", {
  post <- readRDS(test_path("data", "posterior_vanco_thomson_1.rds"))
  expect_output(print(post), "Parameters")
  expect_output(print(post), "Observed data, posterior:")
  expect_output(print(post), "Observed data, prior:")
  expect_output(
    print(post), 
    "posterior_mode[ ]+mean[ ]+median[ ]+sd[ ]+q5[ ]+q95[ ]+rhat"
  )
  expect_output(
    print(post), 
    "time[ ]+dv[ ]+mean[ ]+loc[ ]+pct[ ]+pct5[ ]+pct95"
  )
  expect_output(print(post), "-----|----o")
})
