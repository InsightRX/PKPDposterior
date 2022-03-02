post <- readRDS(test_path("data", "posterior_vanco_thomson_1.rds"))
data <- list(
  ID = c(1, 1, 1, 1, 1, 1, 1),
  time = c(0, 0, 2.5, 11.5, 12, 24, 36),
  cmt = c(2, 2, 1, 1, 2, 2, 2),
  DV = c(0, 0, 40, 14, 0,  0, 0),
  amt = c(1500, 0, 0, 0, 1500, 1500, 1500),
  evid = c(1, 2, 0, 0, 1, 1, 1),
  MDV = c(1, 0, 0, 0, 1, 1, 1),
  rate = c(750, 0, 0, 0, 750, 750, 750),
  WT = c(70, 70, 70, 70, 70, 70, 70),
  CRCL = c(5, 5, 5, 5, 5, 5, 5),
  addl = c(0, 0, 0, 0, 0, 0, 0),
  ss = c(0, 0, 0, 0, 0, 0, 0),
  ii = c(0, 0, 0, 0, 0, 0, 0),
  cObs = c(40, 14),
  nt = 7L,
  iObs = 3:4,
  nObs = 2L
)

test_that("extract_posterior_observation_info works correctly on prior", {
  res <- PKPDposterior:::extract_posterior_observation_info(
    post, 
    data,
    filter = "cHatObs",
    verbose = FALSE
  )
  expect_equal(
    res,
    structure(
      list(
        time = c(2.5, 11.5),
        dv = c(40, 14),
        evid = c(0, 0),
        mean = c(28.6775198, 11.86244634),
        median = c(28.65635, 11.6916),
        sd = c(3.1717696446851, 2.24080293873372),
        pct = c(1, 0.842),
        loc = c("-----|----o", "-----|---o-"),
        pct2.5 = c(23.0385525, 7.88526725),
        pct5 = c(23.678375, 8.408356),
        pct10 = c(24.39475, 9.106655),
        pct25 = c(26.465875, 10.3014),
        pct75 = c(30.78805, 13.150575),
        pct90 = c(32.91628, 14.93326),
        pct95 = c(33.883885, 15.726845),
        pct97.5 = c(35.15159, 16.50444)
      ),
      row.names = c(NA,-2L),
      class = "data.frame"
    )                                                                                          )
})

test_that("extract_posterior_observation_info works correctly on posterior", {
  res <- PKPDposterior:::extract_posterior_observation_info(
    post,
    data,
    filter = "cObsPred",
    verbose = FALSE
  )
  expect_equal(
    res,
    structure(
      list(
        time = c(2.5, 11.5),
        dv = c(40, 14),
        evid = c(0, 0),
        mean = c(29.4704698, 12.50724176),
        median = c(28.64595, 11.76965),
        sd = c(8.84087797989877, 4.39025981823268),
        pct = c(0.878, 0.676),
        loc = c("-----|----o", "-----|--o--"),
        pct2.5 = c(15.0919375, 5.834745),
        pct5 = c(16.39365, 6.6675555),
        pct10 = c(18.67898, 7.605619),
        pct25 = c(23.38215, 9.331605),
        pct75 = c(33.9654, 15.06235),
        pct90 = c(41.87212, 18.2667),
        pct95 = c(46.05685, 19.89446),
        pct97.5 = c(49.2915375, 22.161765)
      ),
      row.names = c(NA,-2L),
      class = "data.frame"
    )
  )
})
