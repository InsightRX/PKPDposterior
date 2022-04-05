post <- readRDS(test_path("data", "posterior_vanco_thomson_1.rds"))

test_that("Extracting stats on predictions from posterior works", {
  res <- PKPDposterior:::extract_from_draws(post)  
  
  expect_equal(nrow(res), 2)
  expect_equal(ncol(res), 18)
  expect_true(all(!is.na(res)))
})
