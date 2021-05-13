
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
compare_sets <- function(x, y, list=FALSE) {
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
#' @importFrom stats setNames
#'
#' @examples{
#'     x <- 1:3
#'     y <- letters[1:3]
#'     xn <- setNames(x,y)
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
  setNames(x, nm)
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
