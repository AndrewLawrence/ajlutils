
#' lmeList_pooled_r.squaredGLMM
#'
#' This function calculates pooled r2 values for multiply-imputed
#' linear mixed effects models.
#'
#' The \code{\link[MuMIn]{r.squaredGLMM}} function from MuMIn package
#' provides r2 values for GLMM such as those produced by
#' \code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
#' and \code{\link[nlme]{lme}}.
#' However, this is not a default method and r2 is not yet pooled
#' when analysing multiply-imputed mixed effects models with mice
#' i.e. \code{mice::pool.r.squared(my_mira)}
#' This function applies the custom r2 function to each model and uses the
#' Harel (2009) fisher-r-to-z method for pooling r2 values.
#'
#' @references Harel, O. (2009).
#' The estimation of R2 and adjusted R2 in incomplete data sets using
#' multiple imputation. Journal of Applied Statistics, 36(10), 1109-1118.
#'
#' @param mlist A list of linear mixed effects models
#'     resulting from multiple imputation
#'
#' @return A matrix of r2 values from \code{\link[MuMIn]{r.squaredGLMM}}
#'
#' @export
lmelist_pooled_rsquaredglmm <- function(mlist) {
  # based on:
  # https://stats.oarc.ucla.edu/stata/faq/how-can-i-estimate-r-squared-for-a-model-estimated-with-multiply-imputed-data/ #nolint
  r2 <- do.call(rbind, lapply(mlist, MuMIn::r.squaredGLMM))
  r <- sqrt(r2)
  z <- apply(psych::fisherz(r), 2, mean)
  r <- psych::fisherz2r(z)
  r2 <- r^2
  return(r2)
}
