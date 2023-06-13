
make_example_data <- function(nk = 5,
                              n = 100,
                              make_labels_arbitrary = TRUE) {
  dat <- matrix("1", nrow = n, ncol = nk)
  n_conv <- floor(n / nk)

  for ( i in seq(from = 2, to = nk ) ) {
    # add the new value (i) for n/nk subjects chosen at random from 1.
    sel <- sample(which(dat[, i - 1] == 1), size = n_conv, replace = FALSE)
    dat[sel, seq(from = i, to = nk)] <- as.character(i)
  }
  dat <- apply(dat, 2, as.character)
  if ( isTRUE(make_labels_arbitrary) ) {
    dat <- apply(dat, 2, function(x) {
      tab <- unique(x)
      tab <- setNames(sample(LETTERS, size = length(tab), replace = FALSE),
                      tab)
      tab[match(x, names(tab))]
    })
  }
  dat
}

dat <- make_example_data()
dfdat <- as.data.frame(dat)

bad_dat1 <- dat[, 5:1]

test_that("construction works", {

  # direct construction:
  expect_no_error(new_hcluslabels(dat))
  # Using as_hcluslabels:
  expect_no_error(as_hcluslabels(dat))
  expect_no_error(as_hcluslabels(dfdat))
  # with validation:
  expect_no_error(validate_hcluslabels(as_hcluslabels(dat)))
  expect_no_error(validate_hcluslabels(as_hcluslabels(dfdat)))
  # no error when invalid data is directly constructed:
  expect_no_error(new_hcluslabels(bad_dat1))

  # However bad data will throw an error when constructed with helper:
  expect_error(as_hcluslabels(bad_dat1))

  # expect errors when we don't construct correctly:
  expect_error(validate_hcluslabels(dat1))

})

test_that("relabelling works", {
  expect_null(check_relabelling(new = relabel_hcluslabels(dfdat), old = dfdat))
})

# non-hierarchical (k1==k2) example:
set.seed(42)
old <- rep(LETTERS[1:4], times = c(10, 20, 30, 40))
new <- rep(letters[c(3, 2, 1, 4)], times = c(10, 20, 30, 40))
# shuffle:
new[sample(seq.int(length(new)), size = 20)] <-
  new[sample(seq.int(length(new)), size = 20)]

test_that("non-heirarchical relabelling works", {
  # runs without error:
  expect_no_error((res <- table(old, relabel_cluster_pair(old, new))))
  # result is maximised along the diagonal:
  expect_true(all(apply(res, 1, which.max) == 1:4))
})
