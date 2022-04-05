test_that("Extracting stats on predictions from posterior works for PK model", {
  post <- readRDS(test_path("data", "posterior_vanco_thomson_1.rds"))
  suppressMessages(
    res <- extract_from_draws(post)  
  )
  expect_equal(nrow(res), 2)
  expect_equal(ncol(res), 18)
  expect_true(all(!is.na(res)))
})

test_that("Extracting stats on predictions from posterior works for PK-PD model", {
  post <- readRDS(test_path("data", "posterior_neutropenia_1.rds"))
  suppressMessages(
    res <- extract_from_draws(post)
  )
  expect_equal(nrow(res), 6)
  expect_equal(res$type, c(rep("pk", 2), rep("pd", 4)))
  expect_equal(ncol(res), 18)
  expect_true(all(!is.na(res)))
})
