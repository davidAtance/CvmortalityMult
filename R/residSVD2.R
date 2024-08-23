#' Function to extract the resid from SVD
#' @description
#' This function uses the first \code{d} components of the singular value
#' decomposition in order to approximate a vector of model residuals by a
#' sum of \code{d} multiplicative terms, with the multiplicative
#' structure determined by two specified factors follows the SVD function by Turner et al. (2023).
#' For \code{glm} and \code{gnm} models from the \code{gnm} R-packages,
#' the matrix entries are weighted working residuals.  The primary use of \code{residSVD} is to
#' generate good starting values for the parameters in \link[gnm]{Mult} terms
#' in models to be fitted using \link[gnm]{gnm}. In this case, we modified the function
#' in order to obtain good starting values for the multi-population mortality models.
#'
#' @param model object with \link[stats]{na.action}, \link[stats]{residuals}, and \link[stats]{weights} methods, e.g. objects inheriting from class \code{"gnm"}.
#' @param fac1 first factor.
#' @param fac2 second factor.
#' @param d integer, the number of multiplicative terms to use in the approximation.
#'
#' @return If \code{d = 1}, a numeric vector; otherwise a numeric matrix with \code{d} columns.
#'
#' @seealso \code{\link{fitLCmulti}}, \code{\link{forecast.fitLCmulti}}
#' \code{\link{multipopulation_cv}}, \code{\link{multipopulation_loocv}},
#' \code{\link{plot.fitLCmulti}}
#'
#' @references
#' Turner, H., & Firth, D. (2023).
#' Generalized nonlinear models in R: An overview of the gnm package.
#' R package version 1.1-5. https://CRAN.R-project.org/package=gnm
#'
#' @import gnm
#'
residSVD2 <- function(model, fac1, fac2, d = 1) {
    if (!is.null(model$call$data)) {
        Data <- as.data.frame(eval(model$call$data, parent.frame()))
        fac1 <- eval(match.call()$fac1, Data, parent.frame())
        fac2 <- eval(match.call()$fac2, Data, parent.frame())
    }
    if (!inherits(model, "glm") && !inherits(model, "lm")) stop(
                                        "model not of class lm, glm or gnm")
    if (!is.factor(fac1)) stop("fac1 must be a factor")
    if (!is.factor(fac2)) stop("fac2 must be a factor")
    Data <- data.frame(fac1, fac2)
    if (!is.null(model$na.action)) Data <- Data[-model$na.action, ]
    weights <- if (!is.null(model$weights)) as.vector(model$weights) else 1
    X <- data.frame(rw = as.vector(model$residuals) * weights, w = weights)
    X <- lapply(X, tapply, Data, sum, simplify = TRUE)
    X <- X$rw/X$w
    X <- svd(naToZero(X), d, d)
    uPart <- sqrt(X$d[seq(d)]) * t(X$u)
    vPart <- sqrt(X$d[seq(d)]) * t(X$v)
#    uPartNegative <- apply(uPart, 1, function(row) all(row < 0))
#    vPartNegative <- apply(vPart, 1, function(row) all(row < 0))
#    multiplier <- ifelse(uPartNegative + vPartNegative == 1, -1, 1)
    multiplier <-  1
    result <- t(cbind(uPart, vPart) * multiplier)
    rownames(result) <- c(paste("fac1", levels(fac1), sep = "."),
                          paste("fac2", levels(fac2), sep = "."))
    colnames(result) <- 1:d
    drop(result)
}
naToZero <- function(vec){
  vec[is.na(vec)] <- 0
  return(vec)
}
