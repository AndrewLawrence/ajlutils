
#' @title not in
#' @name not_in
#'
#' @description utility function, infix operator to negate \code{\link{\%in\%}}.
#'     This simplifies the common "\code{! x \%in\% y}" construction to
#'     "\code{x \%!in\% y}"
#'
#' @inheritParams base::`%in%`
#' @return A logical vector, indicating if each element of x was NOT found in
#'     table, values are TRUE or FALSE and never NA.
#' @export
`%!in%` <- function(x, table) { #nolint!
  match(x, table, nomatch = 0L) == 0L
}


#' Compare two sets
#'
#' Will convert input to sets of unique elements.
#'
#' @param x a vector of elements to compare with y
#' @param y a vector of elements to compare with x
#' @param list boolean should matching and non-matching elements be returned?
#'     Or just a count?
#'
#' @return a list of set sizes,
#'     or a list containing both the set sizes and lists of elements
#'
#' @export
compare_sets <- function(x, y, list = FALSE) {
  x <- unique(as.vector(x))
  y <- unique(as.vector(y))

  nx <- length(x)
  ny <- length(y)
  ni <- sum(x %in% y)
  tab <- list(
    nx,
    ny,
    ni,
    nx - ni,
    ny - ni
  )
  names(tab) <- c(
    "x", "y", "intersection", "only_x", "only_y"
  )
  if (list) {
    dat <- list(
      x[x %in% y],
      x[x %!in% y],
      y[y %!in% x]
    )

    names(dat) <- c(
      "both",
      "only_x",
      "only_y"
    )
    return(list(table = tab,
                data = dat))
  } else {
    return(tab)
  }
}


#' selfname
#'
#' Use the elements of an object as their own names. This code extends
#'     \code{\link[stats]{setNames}}. Mostly useful for making iterables
#'     on the fly.
#'
#' @param x an object (vector, list) with elements coercible to character (and
#'     sensible to use as names).
#' @param append "none", "before", "after" options for appending contents of
#'     x to existing names.
#' @param prefix add this before elements in names
#' @param suffix add this after elements in names
#' @param sep separator to use for constructed names
#'
#' @return a named object x
#'
#' @examples{
#'     x <- 1:3
#'     y <- letters[1:3]
#'     xn <- stats::setNames(x,y)
#'
#'     selfname(x)
#'
#'     selfname(y)
#'
#'     selfname(x, prefix = "v")
#'
#'     selfname(x, prefix = y)
#'
#'     selfname(xn, append = "after")
#' }
#'
#'
#' @export
selfname <- function(x,
                     append = "none",
                     prefix = "",
                     suffix = "",
                     sep = "") {
  at <- attributes(x)
  nm <- as.character(x)

  if ("names" %in% names(at)) {
    append <- match.arg(append, c("none", "before", "after"))
    if (append == "before") {
      nm <- paste(nm, names(x), sep = sep)
    }
    if (append == "after") {
      nm <- paste(names(x), nm, sep = sep)
    }
  }

  if (!identical(prefix, "")) {
    nm <- paste(prefix, nm, sep = sep)
  }
  if (!identical(suffix, "")) {
    nm <- paste(nm, suffix, sep = sep)
  }
  stats::setNames(x, nm)
}

#' get_levels
#'
#' Will extract factor levels - or, if input does not have levels,
#'     appropriate ones will be created.
#'
#' Default R behaviour (as.factor) for a character vector would be to produce
#' alphabetical levels. This can vary by locale. Alternatively we can use
#' the data input order which will not vary by locale. This latter approach is
#' used by forcats::as_factor.
#'
#' Logical vectors will always be assigned c(FALSE, TRUE).
#'
#' @param x input, a factor or character vector.
#' @param drop boolean - should levels not in the data be discarded?
#' @param locale_order boolean - should created factor levels be alphabetical?
#'     Note the default option value (FALSE) is to order as in data
#'     which differs from base R's as.factor, but agrees with
#'     forcats::as_factor. For base-R as.factor equivalence use TRUE.
#'
#' @export
get_levels <- function(x,
                       drop = FALSE,
                       locale_order = FALSE) {
  # Check input:
  if (!inherits(x, what = c("factor", "numeric",
                            "character", "logical",
                            "integer"))) {
    stop("Input should be a factor, character or numeric vector")
  }
  # Objects with levels:
  if (is.factor(x)) {
    if (drop) {
      lvl <- levels(droplevels(x))
    } else {
      lvl <- levels(x)
    }
  } else {
    # don't allow logical data to have factor levels c(TRUE, FALSE).
    #     - matches forcats::as_factor.
    if (is.logical(x)) locale_order <- TRUE

    if (locale_order) {
      lvl <- sort(unique(x))
    } else {
      lvl <- unique(x)
    }
  }
  return(lvl)
}

#' make_stars
#'
#' Produce "significance" stars for p-values.
#'     Defaults correspond to R's stats::summary.lm table but without the '.'
#'     for p-values 0.05 < p < 0.1.
#'
#' Note: because the value of the last element of \code{cuts} is the largest
#'     value that will be assigned a symbol this can be set to 1 to assign
#'     a symbol to non-significant p-values.
#'
#' @param pvals vector of p-values.
#' @param symbols symbols to use
#' @param cuts upper (right) bounds for each symbol (must be in ascending order)
#' @param nonsymbol return value for p-values not assigned a symbol in cuts.
#'     default: NA.
#' @param nablank boolean - should NA values be returned as an empty string:
#'     ""(TRUE), or be passed as an NA (FALSE; default)?
#' @param check boolean - should pvals be checked for range?
#'
#' @seealso stats::symnum
#' @examples {
#'   make_stars(10^c(seq(from = -7, to = 0, by = 0.25)))
#'   # try a custom scheme where you get increasing stars the more zeros
#'   #      you have (with a single "*" for 0.01 < p <= 0.05)
#'   custom_cuts <- 10^c(seq(from = -10, to = -2, by = 1), log10(0.05))
#'   custom_symbols <- sapply(10:1,
#'                            function(x) paste0(rep("*", x), collapse = ""))
#'   make_stars(10^c(seq(from = -7, to = 0, by = 0.25)),
#'              cuts = custom_cuts, symbols = custom_symbols)
#' }
#' @export
make_stars <- function(pvals,
                       symbols = c("***", "**", "*"),
                       cuts = c(0.001, 0.01, 0.05),
                       nonsymbol = NA_character_,
                       nablank = FALSE,
                       check = TRUE) {

  pvals <- as.numeric(pvals)
  if (check && any(pvals > 1, na.rm = TRUE) || any(pvals < 0, na.rm = TRUE)) {
    stop("input must be p-values, between 0 and 1 inclusive")
  }

  cuts <- as.numeric(cuts)
  if (any(is.na(cuts))) stop("cuts must be non-missing.")
  if (!identical(cuts, sort(cuts))) stop("cuts must be sorted ascending")
  if (!identical(cuts, unique(cuts))) stop("cuts must be unique")

  cuts <- unique(c(0, cuts, 1))

  unclass(stats::symnum(pvals,
                        corr = FALSE,
                        na = !nablank,
                        cutpoints = cuts,
                        symbols = c(symbols, nonsymbol), legend = FALSE))
}
