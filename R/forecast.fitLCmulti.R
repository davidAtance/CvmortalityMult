#' Function to forecast multi-population mortality model
#' @description R function for forecasting additive, multiplicative, common-factor (CFM), augmented-common-factor (ACFM), or joint-k multi-population mortality model developed by: Debon et al. (2011), Russolillo et al. (2011), Carter and Lee (1992), LI and Lee (2005), and Carter and Lee (2011), respectively.
#' These models follow the structure of the well-known Lee-Carter model (Lee and Carter, 1992) but include different parameter(s) to capture the behavior of each population considered in different ways.
#' This parameter seeks to capture the individual behavior of every population considered.
#' In case, you want to understand in depth each model, please see Villegas et al. (2017).
#' It should be mentioned that this function is developed for fitting several populations.
#' However, in case you only consider one population, the function will fit the single population version of the Lee-Carter model, the classical one.
#'
#' @param object object \code{"fitLCmulti"} developed using function `fitLCmulti()`. With this object the function will determine the multi-population fitted with the function `fitLCmulti()`.
#' @param nahead number of periods ahead to forecast.
#' @param ktmethod method used to forecast the value of `kt` Arima(p,d,q) or ARIMA(0,1,0); c("`Arimapdq`", "`arima010`").
#' @param ... other arguments for \code{\link[StMoMo:iarima]{iarima}}.
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
#' * `formula` additive multi-population mortality formula used to fit the mortality rates.
#' * `model` provided the model selected in every case.
#' * `qxt.crude` corresponds to the crude mortality rates. These crude rate are directly obtained as: $$q_{x,t,i}=d_{x,t,i}/E_{x,t,i}^{0}$$, with the number of deaths recorded $$d_{x,t,i}$$, and relative to those initially exposed to risk $$E_{x,t,i}$$ for age x, period t and in each region i.
#' * `qxt.fitted` fitted mortality rates using the additive multi-population mortality model.
#' * `logit.qxt.fitted` fitted mortality rates in logit way estimated with the additive multi-population mortality model.
#' * `qxt.future` future mortality rates estimated with the additive multi-population mortality model.
#' * `logit.qxt.future` future mortality rates in logit way estimated with the additive multi-population mortality model.
#' * `nPop` provided number of populations to fit the periods.
#'
#' @seealso \code{\link{fitLCmulti}},
#' \code{\link{plot.fitLCmulti}}, \code{\link{plot.forLCmulti}},
#' \code{\link{multipopulation_cv}}, \link[StMoMo:iarima]{iarima}
#'
#' @references
#'
#' Carter, L.R. and Lee, R.D. (1992).
#' Modeling and forecasting US sex differentials in mortality.
#' International Journal of Forecasting, 8(3), 393–411.
#'
#' Debon, A., Montes, F., & Martinez-Ruiz, F. (2011).
#' Statistical methods to compare mortality for a group with non-divergent populations: an application to Spanish regions.
#' European Actuarial Journal, 1, 291-308.
#'
#' Lee, R.D. & Carter, L.R. (1992).
#' Modeling and forecasting US mortality.
#' Journal of the American Statistical Association, 87(419), 659–671.
#'
#' Li, N. and Lee, R.D. (2005).
#' Coherent mortality forecasts for a group of populations: An extension of the Lee-Carter method.
#' Demography, 42(3), 575–594.
#'
#' Russolillo, M., Giordano, G., & Haberman, S. (2011).
#' Extending the Lee–Carter model: a three-way decomposition.
#' Scandinavian Actuarial Journal, 2011(2), 96-117.
#'
#' Villegas, A. M., Haberman, S., Kaishev, V. K., & Millossovich, P. (2017).
#' A comparative study of two-population models for the assessment of basis risk in longevity hedges.
#' ASTIN Bulletin, 47(3), 631-679.
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
#' #1. ADDITIVE MULTI-POPULATION MORTALITY MODEL
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
#'                                     ktmethod = "Arimapdq")
#'
#' fut_additive_Spainmales
#' #Once the data have been adjusted, it is possible to display the fitted kt and
#' #its out-of-sample forecasting. In addition, the function shows
#' #the logit mortality adjusted in-sample and projected out-of-sample
#' #for the mean age of the data considered in all populations.
#' plot(fut_additive_Spainmales)
#'
#' #2. MULTIPLICATIVE MULTI-POPULATION MORTALITY MODEL
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
#'                                  ktmethod = "Arimapdq")
#'
#' fut_multi_Spainmales
#' #Once the data have been adjusted, it is possible to display the fitted kt and
#' #its out-of-sample forecasting. In addition, the function shows
#' #the logit mortality adjusted in-sample and projected out-of-sample
#' #for the mean age of the data considered in all populations.
#' plot(fut_multi_Spainmales)
#'
#' #3. COMMON-FACTOR MULTI-POPULATION MORTALITY MODEL
#' #In the case, the user wants to fit the common-factor multi-population mortality model
#' cfm_Spainmales <- fitLCmulti(model = "CFM",
#'                              qxt = SpainRegions$qx_male,
#'                              periods = c(1991:2020),
#'                              ages = c(ages),
#'                              nPop = 18,
#'                              lxt = SpainRegions$lx_male)
#'
#' cfm_Spainmales
#'
#' #Once, we have fit the data, it is possible to see the ax, bx, kt, and It
#' #provided parameters for the fitting.
#' plot(cfm_Spainmales)
#'
#' #Once, we have fit the data, it is possible to forecast the multipopulation
#' #mortality model several years ahead, for example 10, as follows:
#' fut_cfm_Spainmales <- forecast(object = cfm_Spainmales, nahead = 10,
#'                                ktmethod = "Arimapdq")
#'
#' fut_cfm_Spainmales
#' #Once the data have been adjusted, it is possible to display the fitted kt and
#' #its out-of-sample forecasting. In addition, the function shows
#' #the logit mortality adjusted in-sample and projected out-of-sample
#' #for the mean age of the data considered in all populations.
#' plot(fut_cfm_Spainmales)
#'
#' #4. AUGMENTED-COMMON-FACTOR MULTI-POPULATION MORTALITY MODEL
#' #In the case, the user wants to fit the augmented-common-factor multi-population mortality model
#' acfm_Spainmales <- fitLCmulti(model = "ACFM",
#'                               qxt = SpainRegions$qx_male,
#'                               periods = c(1991:2020),
#'                               ages = c(ages),
#'                               nPop = 18,
#'                               lxt = SpainRegions$lx_male)
#'
#' acfm_Spainmales
#'
#' #Once, we have fit the data, it is possible to see the ax, bx, kt, and It
#' #provided parameters for the fitting.
#' plot(acfm_Spainmales)
#'
#' #Once, we have fit the data, it is possible to forecast the multipopulation
#' #mortality model several years ahead, for example 10, as follows:
#' fut_acfm_Spainmales <- forecast(object = acfm_Spainmales, nahead = 10,
#'                                ktmethod = "Arimapdq")
#'
#' fut_acfm_Spainmales
#' #Once the data have been adjusted, it is possible to display the fitted kt and
#' #its out-of-sample forecasting. In addition, the function shows
#' #the logit mortality adjusted in-sample and projected out-of-sample
#' #for the mean age of the data considered in all populations.
#' plot(fut_acfm_Spainmales)
#'
#' #5. JOINT-K MULTI-POPULATION MORTALITY MODEL
#' #In the case, the user wants to fit the joint-K multi-population mortality model
#' jointk_Spainmales <- fitLCmulti(model = "joint-K",
#'                                 qxt = SpainRegions$qx_male,
#'                                 periods = c(1991:2020),
#'                                 ages = c(ages),
#'                                 nPop = 18,
#'                                 lxt = SpainRegions$lx_male)
#'
#' jointk_Spainmales
#'
#' #Once, we have fit the data, it is possible to see the ax, bx, kt, and It
#' #provided parameters for the fitting.
#' plot(jointk_Spainmales)
#'
#' #Once, we have fit the data, it is possible to forecast the multipopulation
#' #mortality model several years ahead, for example 10, as follows:
#' fut_jointk_Spainmales <- forecast(object = jointk_Spainmales, nahead = 10,
#'                                ktmethod = "Arimapdq")
#'
#' fut_jointk_Spainmales
#' #Once the data have been adjusted, it is possible to display the fitted kt and
#' #its out-of-sample forecasting. In addition, the function shows
#' #the logit mortality adjusted in-sample and projected out-of-sample
#' #for the mean age of the data considered in all populations.
#' plot(fut_jointk_Spainmales)
#'
#' #6. LEE-CARTER FOR SINGLE-POPULATION
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
#'                               ktmethod = "Arimapdq")
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
                                ktmethod = c("Arimapdq", "arima010"), ...){

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

  if(object$model != "ACFM"){

    if(object$model == "LC-single-pop"){
      nages <- length(object$ax)
    }else{
      nages <- length(object$ax[1,])
    }

    ages <- object$Ages
    lastperiod <- object$Periods[length(object$kt)]+ nahead
    periods <- c((object$Periods[length(object$kt)]+1):lastperiod)

    if(ktmethod == "Arimapdq"){
      kt.a2 <- auto.arima(object$kt, ... )
      fut.kt <- forecast(kt.a2, h = nahead, ...)

      kt.var <- list(kt.arima = kt.a2, mean.kt = fut.kt$mean,
                     lower = fut.kt$lower, upper = fut.kt$upper)

      kt.fut <- matrix(kt.var$mean.kt[1:nahead], nrow= nahead, ncol=1,
                       dimnames= list(c(max(object$Periods+1):(max(object$Periods)+nahead)),"kt"))

    } else if(ktmethod == "arima010"){
      kt.a2 <- Arima(object$kt, order = c(0,1,0),
                     include.drift = T, ...)
      fut.kt <- forecast(kt.a2, h = nahead, ...)

      kt.var <- list(kt.arima = kt.a2, mean.kt = fut.kt$mean,
                     lower = fut.kt$lower, upper = fut.kt$upper)
      kt.fut <- matrix(kt.var$mean.kt[1:nahead], nrow= nahead, ncol=1,
                       dimnames= list(c(max(object$Periods+1):(max(object$Periods)+nahead)),"kt"))

    } else{
      stop("the ktmethod provided does not correspond with 'Arimapdq' or 'arima010'.")
    }
  }else if(object$model == "ACFM"){
    nages <- length(object$ax[1,])
    ages <- (object$Ages)
    lastperiod <- as.numeric(rownames(object$kt)[length(object$kt[,1])])+ nahead
    periods <- c((as.numeric(rownames(object$kt)[length(object$kt[,1])])+1):lastperiod)

    if(ktmethod == "Arimapdq"){
      kt.a2 <- fut.kt <- kt.var <- list()
      kt.fut<- matrix(NA, nrow= nahead, ncol=object$nPop,
                      dimnames= list(c(max(object$Periods+1):(max(object$Periods)+nahead)),c(1:object$nPop)))

      for(pe in 1:object$nPop){
        kt.a2[[paste0("Pop", pe)]] <- auto.arima(object$kt[,pe], ... )
        fut.kt[[paste0("Pop", pe)]] <- forecast(kt.a2[[pe]], h = nahead)
        kt.var[[paste0("Pop", pe)]] <- list(kt.arima = kt.a2[[pe]],
                                            mean.kt = fut.kt[[pe]]$mean,
                                            lower = fut.kt[[pe]]$lower,
                                            upper = fut.kt[[pe]]$upper)
        kt.fut[,pe] <- kt.var[[pe]]$mean.kt[1:nahead]
      }

    } else if(ktmethod == "arima010"){
      kt.a2 <- fut.kt <- kt.var <- list()
      kt.fut<- matrix(NA, nrow= nahead, ncol=object$nPop,
                      dimnames= list(c(max(object$Periods+1):(max(object$Periods)+nahead)),c(1:object$nPop)))
      for(pe in 1:object$nPop){
        kt.a2[[paste0("Pop", pe)]] <- Arima(object$kt[,pe], order = c(0,1,0),
                                            include.drift = T, ...)
        fut.kt[[paste0("Pop", pe)]] <- forecast(kt.a2[[pe]], h = nahead)
        kt.var[[paste0("Pop", pe)]] <- list(kt.arima = kt.a2[[pe]],
                                            mean.kt = fut.kt[[pe]]$mean,
                                            lower = fut.kt[[pe]]$lower,
                                            upper = fut.kt[[pe]]$upper)
        kt.fut[,pe] <- kt.var[[pe]]$mean.kt[1:nahead]
      }
    } else{
      stop("the ktmethod provided does not correspond with 'Arimapdq' or 'arima010'.")
    }
  }


  lee.fut.logit <- list()
  lee.fut.qxt <- list()

  #it <- 1
  if(object$model == "additive"){
    for(it in 1:object$nPop){
      lee.fut.logit[[paste0("pop", it)]] <- matrix(rep(object$ax, nahead), nrow= nages, ncol = nahead)+
        (matrix(object$bx, nrow=nages, ncol=1)%*%matrix(kt.var$mean.kt[1:nahead], nrow=1, ncol=nahead)) +
        matrix(rep(object$Ii[it], nages*nahead), nrow= nages, ncol = nahead)

      lee.fut.qxt[[paste0("pop", it)]] <- plogis(lee.fut.logit[[it]])
      rownames(lee.fut.logit[[it]]) <- rownames(lee.fut.qxt[[it]]) <- ages
      colnames(lee.fut.logit[[it]]) <- colnames(lee.fut.qxt[[it]]) <- periods
    }
  }else if(object$model == "multiplicative"){
    for(it in 1:object$nPop){
      lee.fut.logit[[paste0("pop", it)]] <- matrix(rep(object$ax, nahead), nrow= nages, ncol = nahead)+
        (matrix(object$bx, nrow=nages, ncol=1)%*%matrix(kt.var$mean.kt[1:nahead], nrow=1, ncol=nahead)) *
        matrix(rep(object$Ii[it], nages*nahead), nrow= nages, ncol = nahead)

      lee.fut.qxt[[paste0("pop", it)]] <- plogis(lee.fut.logit[[it]])
      rownames(lee.fut.logit[[it]]) <- rownames(lee.fut.qxt[[it]]) <- ages
      colnames(lee.fut.logit[[it]]) <- colnames(lee.fut.qxt[[it]]) <- periods
    }
  }else if(object$model == "ACFM"){

    lee.fut.logit[[paste0("pop", 1)]] <- matrix(rep(object$ax[1,], nahead), nrow= nages, ncol = nahead) +
      (matrix(object$bx[1,], nrow=nages, ncol=1)%*%matrix(kt.var[[1]]$mean.kt[1:nahead],nrow=1, ncol=nahead))
    lee.fut.qxt[[paste0("pop", 1)]] <- plogis(lee.fut.logit[[1]])

    for(it in 2:object$nPop){
      lee.fut.logit[[paste0("pop", it)]] <- matrix(rep(object$ax[1,], nahead), nrow= nages, ncol = nahead) +
        (matrix(object$bx[1,], nrow=nages, ncol=1)%*%matrix(kt.var[[1]]$mean.kt[1:nahead],nrow=1, ncol=nahead)) +
        matrix(rep(object$ax[it,], nahead), nrow= nages, ncol = nahead) +
        (matrix(object$bx[it,], nrow=nages, ncol=1)%*%matrix(kt.var[[it]]$mean.kt[1:nahead],nrow=1, ncol=nahead))

      lee.fut.qxt[[paste0("pop", it)]] <- plogis(lee.fut.logit[[it]])
      rownames(lee.fut.logit[[it]]) <- rownames(lee.fut.qxt[[it]]) <- ages
      colnames(lee.fut.logit[[it]]) <- colnames(lee.fut.qxt[[it]]) <- periods
    }

  }else if(object$model == "CFM"){
    for(it in 1:object$nPop){
      lee.fut.logit[[paste0("pop", it)]] <- matrix(rep(object$ax[it,], nahead), nrow= nages, ncol = nahead) +
        (matrix(object$bx[1,], nrow=nages, ncol=1)%*%matrix(kt.var$mean.kt[1:nahead],nrow=1, ncol = nahead))

      lee.fut.qxt[[paste0("pop", it)]] <- plogis(lee.fut.logit[[it]])
      rownames(lee.fut.logit[[it]]) <- rownames(lee.fut.qxt[[it]]) <- ages
      colnames(lee.fut.logit[[it]]) <- colnames(lee.fut.qxt[[it]]) <- periods
    }

  }else if(object$model == "joint-K"){
    for(it in 1:object$nPop){
      lee.fut.logit[[paste0("pop", it)]] <- matrix(rep(object$ax[it,], nahead), nrow= nages, ncol = nahead) +
        (matrix(object$bx[it,], nrow=nages, ncol=1)%*%matrix(kt.var$mean.kt[1:nahead],nrow=1, ncol = nahead))

      lee.fut.qxt[[paste0("pop", it)]] <- plogis(lee.fut.logit[[it]])
      rownames(lee.fut.logit[[it]]) <- rownames(lee.fut.qxt[[it]]) <- ages
      colnames(lee.fut.logit[[it]]) <- colnames(lee.fut.qxt[[it]]) <- periods
    }

  }else if(object$model == "LC-single-pop"){
    for(it in 1:object$nPop){
      lee.fut.logit[[paste0("pop", it)]] <- matrix(rep(object$ax, nahead), nrow= nages, ncol = nahead)+
        (matrix(object$bx, nrow=nages, ncol=1)%*%matrix(kt.var$mean.kt[1:nahead], nrow=1, ncol=nahead))

      lee.fut.qxt[[paste0("pop", it)]] <- plogis(lee.fut.logit[[it]])
      rownames(lee.fut.logit[[it]]) <- rownames(lee.fut.qxt[[it]]) <- ages
      colnames(lee.fut.logit[[it]]) <- colnames(lee.fut.qxt[[it]]) <- periods

  }}

  return <- list(ax = object$ax,
                 bx = object$bx,
                 arimakt = kt.a2,
                 kt.fitted = object$kt,
                 kt.fut = kt.fut,
                 kt.futintervals = kt.var,
                 Ii = object$Ii,
                 ktmethod = ktmethod,
                 formula = object$formula,
                 model = object$model,
                 qxt.crude = object$qxt.crude,
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
    } else if(x$model == "CFM"){
      cat("Forecasting the common-factor multi-population mortality model: \n")
    } else if(x$model == "ACFM"){
      cat("Forecasting the augmented-common-factor multi-population mortality model: \n")
    } else if(x$model == "joint-K"){
      cat("Forecasting the joint-K multi-population mortality model: \n")
    }
  } else if(x$nPop == 1){
    cat("Forecasting the single-population version of the Lee-Carter model: \n")
  }

  print(x$formula)

  cat(paste("\nFitting periods:", min(x$Periods), "-", max(x$Periods)))
  cat(paste("\nForecasting periods:", min(x$FutPeriods), "-", max(x$FutPeriods)))
  cat(paste("\nForecasting and Fitting ages:", min(x$Ages), "-", max(x$Ages), "\n"))
  cat(paste("\nForecasting:", x$nPop, "populations \n"))
}
