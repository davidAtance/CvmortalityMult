#' MEAN SQUARED ERRORS (MSE)
#' @description
#' R function to estimate the mean squared errors (MSE) for the mortality rates:
#' \deqn{\frac{1}{n}\sum_{x} \sum_{t} \left( qxt1 - qxt2 \right)^2 = \frac{1}{n} SSE}
#' where qxt1 is the real mortality rates `qxt_re`, and qxt2 is the adjusted mortality rates `qxt_aju`.
#' You only have to provide the real value and the fitted or forecasted value for your mortality rates.
#' However, the function is constructed to provide the real value and the fitted or forecasted value of your independent variable.
#' These variables must have the same dimensions to be compared.
#'
#' @param qxt_re real mortality rates used to check the goodness of fit measure.
#' @param qxt_aju adjusted mortality rates using a specific mode.
#' @param wxt weights of the mortality rates or data provided.
#'
#' @return A value of MSE for the data provided.
#'
#' @seealso \code{\link{fit_additive.LC.multi}}, \code{\link{fit_multiplicative.LC.multi}},
#' \code{\link{for_additive.LC.multi}}, \code{\link{for_multiplicative.LC.multi}},
#' \code{\link{MAE}}, \code{\link{SSE}}, \code{\link{MAPE}},
#' \code{\link{multipopulation_cv}},  \code{\link{multipopulation_loocv}}.
#'
#' @references
#'
#' Atance, D., Debón, A., & Navarro, E. (2020).
#' A comparison of forecasting mortality models using resampling methods.
#' Mathematics, 8(9), 1550.
#'
#' @importFrom StMoMo genWeightMat
#'
#' @examples
#' #The example takes more than 5 seconds because it includes
#' #several fitting and forecasting process and hence all
#' #the process is included in dontrun
#' \dontrun{
#' #To show how the function works, we need to provide fitted or forecasted data and the real data.
#' #In this case, we employ the following data of the library:
#' SpainRegions
#' library(gnm)
#' ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)
#' #In this case, we fit for males providing the lxt
#' multiplicative_Spainmales <- fit_multiplicative.LC.multi(qxt = SpainRegions$qx_male,
#'                               periods = c(1991:2020),
#'                               ages = c(ages),
#'                               nPop = 18,
#'                               lxt = SpainRegions$lx_male)
#'
#' #Once, we have the fitted data, we will obtain the SSE for the first population.
#' #We need to obtain wxt (weight of the mortality rates or data provided) using a
#' library(StMoMo)
#' wxt_1pob <- genWeightMat(ages = ages, years = c(1991:2020), clip = 0)
#' MSE(qxt_re = multiplicative_Spainmales$qxt.real$pob1,
#'     qxt_aju = multiplicative_Spainmales$qxt.fitted$pob1,
#'     wxt = wxt_1pob)
#' }
#' @export
MSE <- function(qxt_re, qxt_aju, wxt){
  ind <- (wxt > 0)
  res <- array(NA, dim = dim(wxt))
  res[ind] <- (qxt_re[ind] - qxt_aju[ind])^2
  res[ind] <- replace(res[ind], res[ind] == "Inf", 0)
  res[ind] <- replace(res[ind], res[ind] == "-Inf", 0)
  res[ind] <- replace(res[ind], res[ind] == "NA", 0)
  res[ind] <- replace(res[ind], res[ind] == "NaN", 0)
  n <- length(qxt_re)
  sum(res)/n
}

