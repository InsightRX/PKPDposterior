post <- readRDS(test_path("data", "posterior_vanco_thomson_1.rds"))
data <- list(
  ID = c(1, 1, 1, 1, 1, 1, 1),
  time = c(0, 0, 2.5, 11.5, 12, 24, 36),
  cmt = c(2, 2, 1, 1, 2, 2, 2),
  DV = c(0, 0, 40, 14, 0, 0, 0),
  amt = c(1500, 0, 0, 0, 1500, 1500, 1500),
  evid = c(1, 2, 0, 0, 1, 1, 1),
  MDV = c(1, 0, 0, 0, 1, 1, 1),
  rate = c(750, 0, 0, 0, 750, 750, 750),
  WT = c(70, 70, 70, 70, 70, 70, 70),
  CRCL = c(5, 5, 5, 5, 5, 5, 5),
  addl = c(0, 0, 0, 0, 0, 0,0),
  ss = c(0, 0, 0, 0, 0, 0, 0),
  ii = c(0, 0, 0, 0, 0, 0, 0),
  theta_CL = 2.99,
  omega_CL = 0.27,
  theta_Q = 2.28,
  omega_Q = 0.49,
  theta_V2 = 0.732,
  omega_V2 = 1.3,
  theta_V1 = 0.675,
  omega_V1 = 0.15,
  dv_pk = c(40, 14),
  i_obs_pk = 3:4,
  n_obs_pk = 2L,
  ruv_prop_pk = 0.15,
  ruv_add_pk = 1.6,
  ltbs_pk = FALSE,
  n_t = 7L
)

test_that("extract_from_draws works correctly on posterior", {
  res <- PKPDposterior:::extract_from_draws(
    post,
    data,
    verbose = FALSE
  )
  
  expect_equal(res,
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
                 class = "data.frame")
  )
})
