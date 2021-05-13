test_that("make_stars works", {

    exponents <- seq(from = -5, to = 0, length.out = 10)
    p <- 10^exponents

    expect_equal(make_stars(p),
                 c("***", "***", "***", "***", "**", "**", "*", NA, NA, NA))

    # try a custom scheme where you get increasing stars the more zeros you have
    #     (with a single "*" for 0.01 < p <= 0.05)
    custom_cuts <- 10^c(seq(from = -10, to = -2, by = 1), log10(0.05))
    custom_symbols <- sapply(10:1,
                             function(x) paste0(rep("*", x), collapse = ""))

    custom_stars <- make_stars(p, cuts = custom_cuts, symbols = custom_symbols)

    expect_equal(max(sapply(custom_stars, nchar), na.rm = T), 5L)
    expect_equal(which.max(sapply(custom_stars, nchar)), 1L, ignore_attr = TRUE)
})
