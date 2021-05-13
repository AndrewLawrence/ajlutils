test_that("selfname works", {

  x <- 1:3
  y <- letters[1:3]
  xn <- setNames(x, y)

  expect_equal(selfname(x), setNames(x, as.character(x)))
  expect_equal(selfname(y), structure(y, names = y))
  expect_equal(selfname(x, prefix = "v"),
               structure(x, names = paste0("v", x)))
  expect_equal(selfname(xn, append = "before"),
               structure(x, names = paste0(x, y)))
})
