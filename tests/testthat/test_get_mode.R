test_that("Mode calculation works for numeric vector", {
  expect_equal(round(get_mode(c(1,2,3,4,5,6,7,8,9)),2), 5.02) 
})
test_that("Mode calculation fails when NULL", {
  expect_error(get_mode(NULL))
})
