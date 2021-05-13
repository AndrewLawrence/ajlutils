test_that("get_levels works", {

  expect_equal(letters[1:10],
               get_levels(letters[1:10]))
  expect_equal(letters[1:10],
               get_levels(letters[10:1], locale_order = T))
  expect_equal(letters[1:5],
               get_levels(as.factor(letters[1:10])[1:5], drop = T))

  expect_error(get_levels(data.frame(letters[1:10])))

})
