parameters <- list(KA = 1.19, VMAX1 = 114, 
                   KM = 1.15, CL = 6.16, Q = 15.5, 
                   V2 = 103, V1 = 79)

ode <- c(
  "real k10 = CL / V1;",
  "real k12 = Q / V1;",
  "real k21 = Q / V2;",
  "real VMAXINH = exp(1.5) / (1 + exp(1.5));",
  "real Vmax = VMAX1 * (1 - VMAXINH * (t-1) / ((t-1) + (2.41 - 1)));",
  "dAdt[1] = -KA*A[1];",
  "dAdt[2] =  KA*A[1] - (k10 + k12)*A[2] + k21*A[3] - (Vmax * A[2]/V1)/(KM + A[2]/V1);",
  "dAdt[3] = k12*A[2] - k21*A[3];"
)

test_that("ode funnction generator works", {
  res <- new_ode_function(
    ode = ode,
    parameters = parameters,
    n_cmt = 3
  )
  expect_equal(length(res), 8)
  expect_equal(grep("dAdt", res), c(3,5,6))
})