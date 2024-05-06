#' Measures of Accuracy
#' @description
#' R function to estimate different measures of accuracy.
#' 1. the sum of squared errors (SSE) for the mortality rates:
#' \deqn{\sum_{x}^{} \sum_{t} \left( qxt1 - qxt2 \right)^{2}}
#' where qxt1 is the real mortality rates `qxt_re`, and qxt2 is the adjusted mortality rates `qxt_aju`.
#' 2. The mean squared errors (MSE) for the mortality rates:
#' \deqn{\frac{1}{n}\sum_{x} \sum_{t} \left( qxt1 - qxt2 \right)^2 = \frac{1}{n} SSE}
#' where qxt1 is the real mortality rates `qxt_re`, and qxt2 is the adjusted mortality rates `qxt_aju`.
#' 3. The mean absolute errors (MAE) for the mortality rates:
#' \deqn{\frac{1}{n}\sum_{x} \sum_{t} \left| qx1 - qxt2 \right|}.
#' where qxt1 is the real mortality rates `qxt_re`, and qxt2 is the adjusted mortality rates `qxt_aju`.
#' 4. The mean absolute percentage error (MAPE) for the mortality rates:
#' \deqn{\frac{1}{n}\sum_{x} \sum_{t}\left| \frac{\left(qxt1 - qxt2\right) }{qxt2} \right|}
#' where qxt1 is the real mortality rates `qxt_re`, and qxt2 is the adjusted mortality rates `qxt_aju`.
#' You only have to provide the real value, the fitted or forecasted value for your mortality rates and the measure of accuracy chosen.
#' However, the function is constructed to provide the real value and the fitted or forecasted value of your independent variable.
#' These variables must have the same dimensions to be compared.
#'
#' @param measure choose the non-penalized measure of accuracy that you want to use; c("`SSE`", "`MSE`", "`MAE`", "`MAPE`", "`All`"). Check the function. In case you do not provide any value, the function will apply the "`SSE`" as measure of forecasting accuracy.
#' @param qxt_re real mortality rates used to check the goodness of fit measure.
#' @param qxt_aju adjusted mortality rates using a specific mode.
#' @param wxt weights of the mortality rates or data provided.
#'
#' @return An object with class \code{"MoA"} including the value of the measure of accuracy for the data provided.
#'
#' @seealso \code{\link{fitLCmulti}}, \code{\link{forecast.fitLCmulti}},
#' \code{\link{multipopulation_cv}}, \code{\link{multipopulation_loocv}}.
#'
#' @references
#' Atance, D., Debón, A., & Navarro, E. (2020).
#' A comparison of forecasting mortality models using resampling methods.
#' Mathematics, 8(9), 1550.
#'
#' @importFrom StMoMo genWeightMat
#'
#' @examples
#' #The example takes more than 5 seconds because it includes
#' #several fitting and forecasting process and hence all
#' #the process is included in donttest
#' \donttest{
#' #To show how the function works, we need to provide fitted or forecasted data and the real data.
#' #In this case, we employ the following data of the library:
#'
#' SpainRegions
#'
#' library(gnm)
#' library(forecast)
#' ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)
#' #In this case, we fit for males providing the lxt
#' multiplicative_Spainmales <- fitLCmulti(model = "multiplicative",
#'                               qxt = SpainRegions$qx_male,
#'                               periods = c(1991:2020),
#'                               ages = c(ages),
#'                               nPop = 18,
#'                               lxt = SpainRegions$lx_male)
#'
#' multiplicative_Spainmales
#' plot(multiplicative_Spainmales)
#'
#' #Once, we have the fitted data, we will obtain different measures of accuracy
#' #for the first population.
#' #We need to obtain wxt (weight of the mortality rates or data provided) using a
#' library(StMoMo)
#' wxt_1pob <- genWeightMat(ages = ages, years = c(1991:2020), clip = 0)
#'
#' ##########################
#' #SSE#
#' ##########################
#' SSE_multSpmales <- MeasureAccuracy(measure = "SSE",
#'                        qxt_re = multiplicative_Spainmales$qxt.real$pob1,
#'                        qxt_aju = multiplicative_Spainmales$qxt.fitted$pob1,
#'                        wxt = wxt_1pob)
#' SSE_multSpmales
#' ##########################
#' #MSE#
#' ##########################
#' MSE_multSpmales <- MeasureAccuracy(measure = "MSE",
#'                        qxt_re = multiplicative_Spainmales$qxt.real$pob1,
#'                        qxt_aju = multiplicative_Spainmales$qxt.fitted$pob1,
#'                        wxt = wxt_1pob)
#' MSE_multSpmales
#' ##########################
#' #MAE#
#' ##########################
#' MAE_multSpmales <- MeasureAccuracy(measure = "MSE",
#'                        qxt_re = multiplicative_Spainmales$qxt.real$pob1,
#'                        qxt_aju = multiplicative_Spainmales$qxt.fitted$pob1,
#'                        wxt = wxt_1pob)
#' MAE_multSpmales
#' ##########################
#' #MAPE#
#' ##########################
#' MAPE_multSpmales <- MeasureAccuracy(measure = "MSE",
#'                         qxt_re = multiplicative_Spainmales$qxt.real$pob1,
#'                         qxt_aju = multiplicative_Spainmales$qxt.fitted$pob1,
#'                         wxt = wxt_1pob)
#' MAPE_multSpmales
#'
#' }
#'
#' @export
MeasureAccuracy <- function(measure = c("SSE", "MSE", "MAE", "MAPE", "All"),
                qxt_re, qxt_aju, wxt){
  valid_measures <- c("SSE", "MSE", "MAE", "MAPE", "All")
  measures <- match.arg(measure, valid_measures)

  if(measure == "SSE"){
    ind <- (wxt > 0)
    res <- array(NA_real_, dim = dim(wxt))
    res[ind] <- (qxt_re[ind] - qxt_aju[ind])^2
    res[ind] <- replace(res[ind], res[ind] == "Inf", 0)
    res[ind] <- replace(res[ind], res[ind] == "-Inf", 0)
    res[ind] <- replace(res[ind], res[ind] == "NA", 0)
    res[ind] <- replace(res[ind], res[ind] == "NaN", 0)
    result <- sum(res)
  }else if(measure == "MSE"){
    ind <- (wxt > 0)
    res <- array(NA_real_, dim = dim(wxt))
    res[ind] <- (qxt_re[ind] - qxt_aju[ind])^2
    res[ind] <- replace(res[ind], res[ind] == "Inf", 0)
    res[ind] <- replace(res[ind], res[ind] == "-Inf", 0)
    res[ind] <- replace(res[ind], res[ind] == "NA", 0)
    res[ind] <- replace(res[ind], res[ind] == "NaN", 0)
    n <- length(qxt_re)
    result <- sum(res)/n
  }else if(measure == "MAE"){
    ind <- (wxt > 0)
    res <- array(NA_real_, dim = dim(wxt))
    res[ind] <- abs(qxt_re[ind] - qxt_aju[ind])
    res[ind] <- replace(res[ind], res[ind] == "Inf", 0)
    res[ind] <- replace(res[ind], res[ind] == "-Inf", 0)
    res[ind] <- replace(res[ind], res[ind] == "NA", 0)
    res[ind] <- replace(res[ind], res[ind] == "NaN", 0)
    n <- length(qxt_re)
    result <- sum(res)/n
  }else if(measure == "MAPE"){
    ind <- (wxt > 0)
    res <- array(NA_real_, dim = dim(wxt))
    res[ind] <- abs((qxt_re[ind] - qxt_aju[ind])/qxt_re[ind])
    res[ind] <- replace(res[ind], res[ind] == "Inf", 0)
    res[ind] <- replace(res[ind], res[ind] == "-Inf", 0)
    res[ind] <- replace(res[ind], res[ind] == "NA", 0)
    res[ind] <- replace(res[ind], res[ind] == "NaN", 0)
    n <- length(qxt_re)
    result <- sum(res)/n
  }else if(measure == "All"){
    ind1 <- ind2 <- ind3 <- ind4 <- (wxt > 0)
    res1 <- res2 <- res3 <- res4 <- array(NA_real_, dim = dim(wxt))

    res1[ind1] <- (qxt_re[ind1] - qxt_aju[ind1])^2
    res1[ind1] <- replace(res1[ind1], res1[ind1] == "Inf", 0)
    res1[ind1] <- replace(res1[ind1], res1[ind1] == "-Inf", 0)
    res1[ind1] <- replace(res1[ind1], res1[ind1] == "NA", 0)
    res1[ind1] <- replace(res1[ind1], res1[ind1] == "NaN", 0)
    result1 <- sum(res1)

    res2[ind2] <- (qxt_re[ind2] - qxt_aju[ind2])^2
    res2[ind2] <- replace(res[ind2], res[ind2] == "Inf", 0)
    res2[ind2] <- replace(res[ind2], res[ind2] == "-Inf", 0)
    res2[ind2] <- replace(res[ind2], res[ind2] == "NA", 0)
    res2[ind2] <- replace(res[ind2], res[ind2] == "NaN", 0)
    n2 <- length(qxt_re)
    result2 <- sum(res2)/n2

    res3[ind3] <- abs(qxt_re[ind3] - qxt_aju[ind3])
    res3[ind3] <- replace(res3[ind3], res3[ind3] == "Inf", 0)
    res3[ind3] <- replace(res3[ind3], res3[ind3] == "-Inf", 0)
    res3[ind3] <- replace(res3[ind3], res3[ind3] == "NA", 0)
    res3[ind3] <- replace(res3[ind3], res3[ind3] == "NaN", 0)
    n3 <- length(qxt_re)
    result3 <- sum(res3)/n3

    res4[ind4] <- abs((qxt_re[ind4] - qxt_aju[ind4])/qxt_re[ind4])
    res4[ind4] <- replace(res4[ind4], res4[ind4] == "Inf", 0)
    res4[ind4] <- replace(res4[ind4], res4[ind4] == "-Inf", 0)
    res4[ind4] <- replace(res4[ind4], res4[ind4] == "NA", 0)
    res4[ind4] <- replace(res4[ind4], res4[ind4] == "NaN", 0)
    n4 <- length(qxt_re)
    result4 <- sum(res4)/n4

    result <- c(result1, result2, result3, result4)
    result <- matrix(result, nrow = 1, ncol = 4,
                     dimnames = list("value", c("SSE", "MSE", "MAE", "MAPE")))
  }

  return <- list(value = result,
                 measure = measure)
  class(return) <- "MoA"
  return
}
#' @export
print.MoA <- function(x, ...) {
  if(!is.null(x)){
    if(!"MoA" %in% class(x))
      stop("The 'x' does not have the 'MoA' structure of R CvmortalityMult package.")
  }
  if(!is.list(x)){
    stop("The 'x' is not a list. Use 'MoA' function first.")
  }

  cat("Measure of Accuracy employed:", x$measure,"\n")
  cat("The error obtained is:\n")
  print(round(x$value, 6) )
}
