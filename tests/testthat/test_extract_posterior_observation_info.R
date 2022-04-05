post <- readRDS(test_path("data", "posterior_vanco_thomson_1.rds"))

test_that("extract_from_draws works correctly on posterior", {
  res <- PKPDposterior:::extract_from_draws(
    post,
    verbose = FALSE
  )
  
  expect_equal(
    res,
    structure(
      list(
        time = c(2.5, 11.5),
        dv = c(40, 14),
        evid = c(0,
                 0),
        type = c("pk", "pk"),
        variable = c("ipred_obs_pk[1]", "ipred_obs_pk[2]"),
        mean = c(30.8252078, 12.4667503),
        median = c(30.69735, 12.279),
        sd = c(2.9966596159034, 2.20902272172831),
        pct = c(0.994,
                0.766),
        loc = c("-----|----o", "-----|---o-"),
        pct2.5 = c(25.4380725,
                   8.675681),
        pct5 = c(26.183275, 9.223874),
        pct10 = c(27.12148,
                  9.814762),
        pct25 = c(28.75215, 10.951375),
        pct75 = c(32.9174,
                  13.8624),
        pct90 = c(34.77864, 15.46311),
        pct95 = c(36.0434, 16.1921),
        pct97.5 = c(37.0050175, 16.975405)
      ),
      row.names = c(NA, -2L),
      class = "data.frame"
    )
  )
})
