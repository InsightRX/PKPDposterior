test_that("remap() renames correctly, removes old names", {
  lst1 <- list("A" = 1, "B" = 2, "C" = 3)
  expect_equal(
    remap(lst1, list("D" = "A", "E" = "B")),
    list(C = 3, D = 1, E = 2)
  )
})

test_that("remap() doesn't rename/remove if not included in map", {
  lst2 <- list(X = 1, Y = 1)
  expect_equal(
    remap(lst2, list("D" = "A", "E" = "B")),
    lst2
  )
})