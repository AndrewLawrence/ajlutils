#' extended_coefficients
#'
#' Model coefficients can be hard to interpret if there is collinearity in the
#' predictors.
#' This function provides an expanded coefficients table for a typical lm
#'     regression model. Inspired by SPSS this presents:
#'     \itemize{
#'      \item standardised model coefficients (standardised betas)
#'      \item structure coefficients (correlation with model predictions,
#'          \code{cor(y-hat, x)})
#'      \item zero-order correlations (i.e. \code{cor(y, x)})
#'      \item part (semi-partial) correlations (i.e. \code{cor(y, x|z)})
#'      \item partial correlations (i.e. \code{cor(y|z, x|z)})
#'     }
#'     where \code{y} is the outcome, \code{x} is the predictor of interest,
#'     and \code{z} is the set of all other predictors in the model
#'
#'@param dat a data.frame or tibble
#'@param yvar character name of the y-variable in \code{dat} used as the
#'    outcome in the lm model
#'@param xvars character vector of the column names in \code{dat} to
#'    use as x-variables in in the lm model
#'@return A named list object containing:
#'    \itemize{
#'    \item *corrmat* a correlation matrix for the predictors
#'    \item *coefficients* an extended coefficients table
#'    }
#'@examples
#'extended_coefficients(mtcars, "mpg", c("cyl", "disp"))
#'@importFrom stats cor lm coef predict residuals
#'@export
extended_coefficients <- function(dat,
                                  yvar = "y",
                                  xvars = c("X1", "X2")) {
  checkmate::assert_data_frame(dat)
  checkmate::assert_character(yvar, len = 1L)
  checkmate::assert_character(xvars, min.len = 1L)
  checkmate::assert_true(all(c(yvar, xvars) %in% colnames(dat)))

  # fit the model:
  m0_f <- paste0(
    "scale(", yvar, ") ~ ",
    paste0("scale(", xvars, ")", collapse = " + ")
  )
  m0_f
  mod <- lm(formula = m0_f, data = dat)

  # semi-partial data prep:
  sp_m <- as.data.frame(
    lapply(
      setNames(seq_along(xvars), xvars),
      \(.x) {
        f <- paste0(xvars[.x], " ~ ", paste0(xvars[-c(.x)], collapse = " + "))
        residuals(lm(formula = f, data = dat))
      }
    )
  )

  # model coefficients:
  m0_coef <- setNames(as.vector(coef(mod)[-1]), xvars)
  # structure coefficients:
  sc <- setNames(as.vector(cor(predict(mod), dat[, xvars])), xvars)
  # zero order correlations (with y):
  zeroorder <- setNames(as.vector(cor(dat[, yvar], dat[, xvars])), xvars)
  # semi-partial (part):
  semipartial <- setNames(as.vector(cor(dat[, yvar], sp_m[, xvars])), xvars)
  # partial:
  partial <- psych::partial.r(data = dat[, c(yvar, xvars)])[yvar, -1]

  # assemble:
  res <- list(
    coef = m0_coef,
    sc = sc,
    zeroorder = zeroorder,
    semipartial = semipartial,
    partial = partial
  )
  res <- do.call(cbind, res)
  res <- as.data.frame(res)
  nm <- colnames(res)
  res$term <- rownames(res)
  res <- res[, c("term", nm)]
  rownames(res) <- NULL

  list(corrmat = cor(dat[, xvars]),
       coefficients = res)
}
