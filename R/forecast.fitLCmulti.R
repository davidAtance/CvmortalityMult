#' Function to forecast multi-population mortality model
#' @description R function for forecasting additive and multiplicative multi-population mortality model developed by: Debon et al (2011) and Russolillo et al. (2011), respectively.
#' This model follows the structure of the well-known Lee-Carter model (Lee and Carter, 1992) but including an additive or multiplicative parameter to capture the behavior of each population considered.
#' This parameter seeks to capture the individual behavior of every population considered.
#' It should be mentioned that this function is developed for fitting several populations.
#' However, in case you only consider one population, the function will fit the single population version of the Lee-Carter model, the classical one.
#'
#' @param object object \code{"fitLCmulti"} developed using function `fitLCmulti()`. With this object the function will determine the multi-population fitted with the function `fitLCmulti()`.
#' @param nahead number of periods ahead to forecast.
#' @param ktmethod method used to forecast the value of `kt` Arima(p,d,q) or ARIMA(0,1,0); c("`Arimapdq`", "`arima010`").
#' @param kt_include.cte if you want that `kt` include constant in the arima process.
#' @param ... other arguments for \code{\link{iarima}}.
#'
#' @return A list with class \code{"forLCmulti"} including different components of the forecasting process:
#' * `ax` parameter that captures the average shape of the mortality curve in all considered populations.
#' * `bx` parameter that explains the age effect x with respect to the general trend `kt` in the mortality rates of all considered populations.
#' * `arimakt` the arima selected for the `kt` time series.
#' * `kt.fitted` obtained values for the tendency behavior captured by `kt`.
#' * `kt.fut` projected values of `kt` for the nahead periods ahead.
#' * `kt.futintervals` arima selected and future values of `kt` with the different intervals, lower and upper, 80\% and 90\%.
#' * `Ii` parameter that captures the differences in the pattern of mortality in any region i with respect to Region 1.
#' * `ktmethod` method selected to forecast the value of `kt` Arima(p,d,q) or ARIMA(0,1,0); c("`Arimapdq`", "`arima010`").
#' * `kt_include.cte` the decision regarding the inclusion of constant in the `kt` arima process.
#' * `formula` additive multi-population mortality formula used to fit the mortality rates.
#' * `model` provided the model selected in every case.
#' * `qxt.real` real mortality rates.
#' * `qxt.fitted` fitted mortality rates using the additive multi-population mortality model.
#' * `logit.qxt.fitted` fitted mortality rates in logit way estimated with the additive multi-population mortality model.
#' * `qxt.future` future mortality rates estimated with the additive multi-population mortality model.
#' * `logit.qxt.future` future mortality rates in logit way estimated with the additive multi-population mortality model.
#' * `nPop` provided number of populations to fit the periods.
#'
#' @seealso \code{\link{fitLCmulti}},
#' \code{\link{plot.fitLCmulti}}, \code{\link{plot.forLCmulti}},
#' \code{\link{multipopulation_cv}}, \code{\link{multipopulation_loocv}}
#'
#' @references
#' Debon, A., Montes, F., & Martinez-Ruiz, F. (2011).
#' Statistical methods to compare mortality for a group with non-divergent populations: an application to Spanish regions.
#' European Actuarial Journal, 1, 291-308.
#'
#' Lee, R.D. & Carter, L.R. (1992).
#' Modeling and forecasting US mortality.
#' Journal of the American Statistical Association, 87(419), 659–671.
#'
#' Russolillo, M., Giordano, G., & Haberman, S. (2011).
#' Extending the Lee–Carter model: a three-way decomposition.
#' Scandinavian Actuarial Journal, 2011(2), 96-117.
#'
#' @importFrom forecast Arima auto.arima forecast
#' @importFrom utils install.packages
#' @importFrom stats plogis qlogis
#'
#' @examples
#' #The example takes more than 5 seconds because it includes
#' #several fitting and forecasting process and hence all
#' #the process is included in donttest
#' \donttest{
#' #First, we present the data that we are going to use
#' SpainRegions
#' ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)
#'
#' library(gnm)
#' library(forecast)
#' #ADDITIVE MULTI-POPULATION MORTALITY MODEL
#' #In the case, the user wants to fit the additive multi-population mortality model
#' additive_Spainmales <- fitLCmulti(model = "additive",
#'                               qxt = SpainRegions$qx_male,
#'                               periods = c(1991:2020),
#'                               ages = c(ages),
#'                               nPop = 18,
#'                               lxt = SpainRegions$lx_male)
#'
#' additive_Spainmales
#'
#' #If the user does not provide the model inside the function fitLCmult()
#' #the multi-population mortality model applied will be additive one.
#'
#' #Once, we have fit the data, it is possible to see the ax, bx, kt, and Ii
#' #provided parameters for the fitting.
#' plot(additive_Spainmales)
#'
#' #Once, we have fit the data, it is possible to forecast the multipopulation
#' #mortality model several years ahead, for example 10, as follows:
#' fut_additive_Spainmales <- forecast(object = additive_Spainmales, nahead = 10,
#'                                     ktmethod = "Arimapdq", kt_include.cte = TRUE)
#'
#' fut_additive_Spainmales
#' #Once the data have been adjusted, it is possible to display the fitted kt and
#' #its out-of-sample forecasting. In addition, the function shows
#' #the logit mortality adjusted in-sample and projected out-of-sample
#' #for the mean age of the data considered in all populations.
#' plot(fut_additive_Spainmales)
#'
#' #MULTIPLICATIVE MULTI-POPULATION MORTALITY MODEL
#' #In the case, the user wants to fit the multiplicative multi-population mortality model
#' multiplicative_Spainmales <- fitLCmulti(model = "multiplicative",
#'                               qxt = SpainRegions$qx_male,
#'                               periods = c(1991:2020),
#'                               ages = c(ages),
#'                               nPop = 18,
#'                               lxt = SpainRegions$lx_male)
#'
#' multiplicative_Spainmales
#'
#' #Once, we have fit the data, it is possible to see the ax, bx, kt, and It
#' #provided parameters for the fitting.
#' plot(multiplicative_Spainmales)
#'
#' #Once, we have fit the data, it is possible to forecast the multipopulation
#' #mortality model several years ahead, for example 10, as follows:
#' fut_multi_Spainmales <- forecast(object = multiplicative_Spainmales, nahead = 10,
#'                                  ktmethod = "Arimapdq", kt_include.cte = TRUE)
#'
#' fut_multi_Spainmales
#' #Once the data have been adjusted, it is possible to display the fitted kt and
#' #its out-of-sample forecasting. In addition, the function shows
#' #the logit mortality adjusted in-sample and projected out-of-sample
#' #for the mean age of the data considered in all populations.
#' plot(fut_multi_Spainmales)
#'
#' #LEE-CARTER FOR SINGLE-POPULATION
#' #As we mentioned in the details of the function, if we only provide the data
#' #from one-population the function fitLCmulti()
#' #will fit the Lee-Carter model for single populations.
#' LC_Spainmales <- fitLCmulti(qxt = SpainNat$qx_male,
#'                               periods = c(1991:2020),
#'                               ages = ages,
#'                               nPop = 1)
#'
#' LC_Spainmales
#'
#' #Once, we have fit the data, it is possible to see the ax, bx, and kt
#' #parameters provided for the single version of the LC.
#' plot(LC_Spainmales)
#'
#' #Once, we have fit the data, it is possible to forecast the multipopulation
#' #mortality model several years ahead, for example 10, as follows:
#' fut_LC_Spainmales <- forecast(object = LC_Spainmales, nahead = 10,
#'                               ktmethod = "Arimapdq", kt_include.cte = TRUE)
#'
#' #Once the data have been adjusted, it is possible to display the fitted kt and
#' #its out-of-sample forecasting. In addition, the function shows
#' #the logit mortality adjusted in-sample and projected out-of-sample
#' #for the mean age of the data considered in all populations.
#' plot(fut_LC_Spainmales)
#'
#' }
#'
#' @export
forecast.fitLCmulti <- function(object, nahead,
                                ktmethod = c("Arimapdq", "arima010"),
                                kt_include.cte = TRUE, ...){

  #First check the structure of object is equal to the previous object created using our function
  if(!is.null(object)){
    if(!"fitLCmulti" %in% class(object))
      stop("The object does not have the 'fitLCmulti' structure of R CvmortalityMult package.")
  }
  if(!is.list(object)){
    stop("The object is not a list. Use 'fitLCmulti' function first.")
  }

  if (!is.numeric(nahead)) {
    stop("nahead has to be a numeric variable.")
    if (nahead <= 0) {
      stop("nahead has to be higher than 0.")
    }
  }

  #We will use plogis() for inverse logit
  #plogis

  if(ktmethod == "Arimapdq"){
    kt.a2 <- auto.arima(object$kt, ic = "bic",
                        allowdrift = kt_include.cte)
  } else if(ktmethod == "arima010"){
    kt.a2 <- Arima(object$kt, order = c(0,1,0),
                   include.drift = kt_include.cte)
  } else{
    stop("the ktmethod provided does not correspond with 'Arimapdq' or 'arima010'.")
  }

  fut.kt <- forecast(kt.a2, h = nahead)

  kt.var <- list(kt.arima = kt.a2, mean.kt = fut.kt$mean,
                 lower = fut.kt$lower, upper = fut.kt$upper)

  lee.fut.logit <- list()
  lee.fut.qxt <- list()
  nages <- length(object$ax)
  ages <- colnames(object$ax)
  lastperiod <- as.numeric(rownames(object$kt)[length(object$kt)])+ nahead
  periods <- c((as.numeric(rownames(object$kt)[length(object$kt)])+1):lastperiod)

  #it <- 1
  if(object$model == "additive"){
    for(it in 1:object$nPop){
      lee.fut.logit[[paste0("pob", it)]] <- matrix(rep(object$ax, nahead), nrow= nages, ncol = nahead)+
        (matrix(object$bx, nrow=nages, ncol=1)%*%matrix(kt.var$mean.kt[1:nahead], nrow=1, ncol=nahead)) +
        matrix(rep(object$Ii[it], nages*nahead), nrow= nages, ncol = nahead)

      lee.fut.qxt[[paste0("pob", it)]] <- plogis(lee.fut.logit[[it]])
      rownames(lee.fut.logit[[it]]) <- rownames(lee.fut.qxt[[it]]) <- ages
      colnames(lee.fut.logit[[it]]) <- colnames(lee.fut.qxt[[it]]) <- periods
    }
  }else if(object$model == "multiplicative"){
    for(it in 1:object$nPop){
      lee.fut.logit[[paste0("pob", it)]] <- matrix(rep(object$ax, nahead), nrow= nages, ncol = nahead)+
        (matrix(object$bx, nrow=nages, ncol=1)%*%matrix(kt.var$mean.kt[1:nahead], nrow=1, ncol=nahead)) *
        matrix(rep(object$Ii[it], nages*nahead), nrow= nages, ncol = nahead)

      lee.fut.qxt[[paste0("pob", it)]] <- plogis(lee.fut.logit[[it]])
      rownames(lee.fut.logit[[it]]) <- rownames(lee.fut.qxt[[it]]) <- ages
      colnames(lee.fut.logit[[it]]) <- colnames(lee.fut.qxt[[it]]) <- periods
    }
  }else if(object$model == "LC-single-pop"){
    for(it in 1:object$nPop){
      lee.fut.logit[[paste0("pob", it)]] <- matrix(rep(object$ax, nahead), nrow= nages, ncol = nahead)+
        (matrix(object$bx, nrow=nages, ncol=1)%*%matrix(kt.var$mean.kt[1:nahead], nrow=1, ncol=nahead))

      lee.fut.qxt[[paste0("pob", it)]] <- plogis(lee.fut.logit[[it]])
      rownames(lee.fut.logit[[it]]) <- rownames(lee.fut.qxt[[it]]) <- ages
      colnames(lee.fut.logit[[it]]) <- colnames(lee.fut.qxt[[it]]) <- periods

  }}

  return <- list(ax = matrix(object$ax, nrow = 1, ncol = nages, dimnames = list("ax", ages)),
                 bx = matrix(object$bx, nrow = 1, ncol = nages, dimnames = list("bx", ages)),
                 arimakt = kt.a2,
                 kt.fitted = matrix(object$kt, nrow = length(object$Periods), ncol = 1, dimnames = list(object$Periods, "kt")),
                 kt.fut = matrix(kt.var$mean.kt[1:nahead], nrow= nahead, ncol=1,
                                 dimnames= list(c(max(object$Periods+1):(max(object$Periods)+nahead)),"kt")),
                 kt.futintervals = kt.var,
                 Ii = matrix(object$Ii, nrow = object$nPop, ncol = 1, dimnames = list(c(1:object$nPop), "Ii")),
                 ktmethod = ktmethod,
                 kt_include.cte = kt_include.cte,
                 formula = object$formula,
                 model = object$model,
                 qxt.real = object$qxt.real,
                 qxt.fitted = object$qxt.fitted.qxt,
                 logit.qxt.fitted = object$logit.qxt.fitted,
                 qxt.future = lee.fut.qxt,
                 logit.qxt.future = lee.fut.logit,
                 Ages = object$Ages,
                 Periods = object$Periods,
                 FutPeriods = c((max(object$Periods)+1):(max(object$Periods)+nahead)),
                 nPop = object$nPop)
  class(return) <- "forLCmulti"
  return

}
#' @export
print.forLCmulti <- function(x, ...) {
  if(!is.null(x)){
    if(!"forLCmulti" %in% class(x))
      stop("The 'x' does not have the 'forLCmulti' structure of R CvmortalityMult package.")
  }
  if(!is.list(x)){
    stop("The 'x' is not a list. Use 'forecast.fitLCmulti' function first.")
  }

  if(x$nPop != 1){
    if(x$model == "additive"){
      cat("Forecasting the additive multi-population mortality model: \n")
    } else if(x$model == "multiplicative"){
      cat("Forecasting the multiplicative multi-population mortality model: \n")
    }
  } else if(x$nPop == 1){
    cat("Forecasting the single-population version of the Lee-Carter model: \n")
  }

  print(x$formula)

  cat(paste("\nFitting years:", min(x$Periods), "-", max(x$Periods)))
  cat(paste("\nForecasting years:", min(x$FutPeriods), "-", max(x$FutPeriods)))
  cat(paste("\nForecasting and Fitting ages:", min(x$Ages), "-", max(x$Ages), "\n"))
  cat(paste("\nForecasting:", x$nPop, "populations \n"))
}
