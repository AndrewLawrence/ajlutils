test_that("not in works", {
  x <- 1:10
  y <- 1:5
  expect_equal(x %!in% y, ! x %in% y)
})
