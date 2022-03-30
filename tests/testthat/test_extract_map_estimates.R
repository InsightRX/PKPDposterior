test_that("Calculation of MAP estimates as mode of posterior works", {
  post <- list(
    draws_df = data.frame(
      CL = c(1,2,3,4,5),
      V = c(6,7,8,9,10),
      Q = c(3,4,5,6,7)
    )
  )
  expect_equal(
    extract_map_estimates(post),
    list(CL = 3.00962965532006, V = 7.99037034467994, Q = 5.00962965532006)
  )
})

test_that("Calculation of MAP estimates for specific params", {
  post <- list(
    draws_df = data.frame(
      CL = c(1,2,3,4,5),
      V = c(6,7,8,9,10),
      Q = c(3,4,5,6,7)
    )
  )
  expect_equal(
    extract_map_estimates(post, c("CL", "Q")),
    list(CL = 3.00962965532006, Q = 5.00962965532006)
  )
})

test_that("Calculation of MAP estimates using weighted draws", {
  post <- list(
    draws_df = data.frame(
      CL = c(1,2,3,4,5),
      V = c(6,7,8,9,10),
      Q = c(3,4,5,6,7)
    )
  )
  wts = c(0.4, 0.1, 0.1, 0.2, 0.2)
  expect_equal(
    extract_map_estimates(post, c("CL", "Q"), weights = wts),
    list(CL = 1.2570323870687, Q = 3.2570323870687)
  )
})
