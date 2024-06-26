#' Function to plot the parameters of the multi-population mortality models
#' @description
#' Function to plot different results of the forecasting process of multi-population mortality models, the additive (Debon et al., 2011) and the multiplicative (Russolillo et al., 2011), obtained using the `forecast.fitLCmulti` function which are xs of the `forecastLCmulti` class.
#' In fact, the function will show the trend parameter kt fitted for the in-sample periods and its forecast results. Similarly, the behavior of the logit mortality rate for the mean in-sample age and the out-of-sample forecast will be shown for all the populations considered.
#' It should be mentioned that this function is developed for fitting several populations.
#' However, in case you only consider one population, the function will show the single population version of the Lee-Carter model, the classical one.
#'
#' @param x `x` developed using function `forecast.fitLCmulti()` which are objects of the `fortLCmulti` class.
#' @param ... additional arguments to show in the plot appearance.
#'
#' @return plot the trend parameter kt fitted for the in-sample periods and its forecast results for the multi-population mortality models. Similarly, the behavior of the logit mortality rate for the mean in-sample age and the out-of-sample forecast will be shown for all the populations considered.
#'
#' @seealso \code{\link{fitLCmulti}}, \code{\link{forecast.fitLCmulti}},
#' \code{\link{plot.fitLCmulti}},
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
#' Multi-population mortality model developed by:
#' Russolillo, M., Giordano, G., & Haberman, S. (2011).
#' Extending the Lee–Carter model: a three-way decomposition.
#' Scandinavian Actuarial Journal, 2011(2), 96-117.
#'
#' @importFrom graphics par legend lines mtext
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
#' @export
plot.forLCmulti<- function(x, ...){
  if(!is.null(x)){
    if(!"forLCmulti" %in% class(x))
      stop("The x does not have the 'forLCmulti' structure of R CvmortalityMult package.")
  }
  if(!is.list(x)){
    stop("The x is not a list. Use 'forecast.fitLCmulti' function first.")
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  par(mfrow=c(1,2))

  plot(forecast(auto.arima(x$kt.fitted), h = length(x$FutPeriods)),
       xlab = "periods", ylab = "kt")

  showage <- round(length(x$Ages)/2, 0)
  ageshow <- x$Ages[showage]

  ymin <- ymax <- c()
  for(i in 1:x$nPop){
    ymin <- c(ymin, min(x$logit.qxt.fitted[[i]][showage,],
                        x$logit.qxt.future[[i]][showage,]))
    ymax <- c(ymax, max(x$logit.qxt.fitted[[i]][showage,],
                        x$logit.qxt.future[[i]][showage,]))
  }

  xmin <- min(x$Periods)
  xmax <- max(x$FutPeriods)

  plot(c(x$Periods),  x$logit.qxt.fitted$pob1[showage,], type="l",
       xlim=c(xmin, xmax), ylim = c(min(ymin), max(ymax)),
       main = paste0("Age ", ageshow), xlab = "periods", ylab = "logit qxt",
  )
  lines(c(x$FutPeriods), x$logit.qxt.future$pob1[showage,], col="blue")
  if(x$nPop != 1){
    for(i in 2:x$nPop){
      lines(c(x$Periods), x$logit.qxt.fitted[[(i)]][showage,], col="black", lty = (i+1))
      lines(c(x$FutPeriods), x$logit.qxt.future[[(i)]][showage,], col="blue", lty=(i+1))
    }
  }
  legend("bottomleft", names(x$logit.qxt.future), lty = c(1:x$nPop), cex = 0.4)

  if(x$nPop != 1){
    if(x$model == "additive"){
      text1 <- paste0("Forecasting the additive multi-population mortality model")
    } else if(x$model == "multiplicative"){
      text1 <- paste0("Forecasting the multiplicative multi-population mortality model")
    }
  } else if(x$nPop == 1){
    text1 <- paste0("Forecasting the single-population version of the Lee-Carter model")
  }

  mtext(text1, line = -1.25, outer = TRUE, cex = 1.5, col = "black")

}
