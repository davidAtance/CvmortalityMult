#' FUNCTION TO FORECAST ADDITIVE MULTI-POPULATION MORTALITY MODEL
#' @description R function for forecasting additive multi-population mortality model developed by: Debon et al (2011).
#' This model follows the structure of the well-known Lee-Carter model (Lee and Carter, 1992) but including an additive parameter to capture the behavior of each population considered.
#' This parameter seeks to capture the individual behavior of every population considered.
#' It should be mentioned that in case that this function is developed for forecasting several populations.
#' However, in case you only consider one population, the function will forecast the Lee-Carter model for one population.
#'
#' @param fitted.obj object developed using function fit.additive.LC.multi().
#' @param nahead number of periods ahead to forecast.
#' @param ktmethod method used to forecast the value of kt Arima(p,d,q) or ARIMA(0,1,0); c("Arimapdq", "arima010").
#' @param kt_include.cte if you want that kt include constant in the arima process.
#'
#' @return A list with different components of the forecasting process:
#' * `ax` parameter that captures the average shape of the mortality curve in all considered populations.
#' * `bx` parameter that explains the age effect x with respect to the general trend (kt) in the mortality rates of all considered populations.
#' * `arimakt` the arima selected for the \eqn{k_t} time series.
#' * `kt.fitted` obtained \eqn{k_t} values for the tendency behaviour.
#' * `kt.fut` \eqn{k_t} projected values for the nahead periods ahead.
#' * `kt.futintervals` \eqn{k_t} arima selected and future values for this parameter with the different intervals, lower and upper, 80\% and 90\%.
#' * `Ii` paramater that captures the differences in the pattern of mortality in any region i with respect to Region 1.
#' * `formula` additive multi-population mortality formula used to fit the mortality rates.
#' * `qxt.real` real mortality rates.
#' * `qxt.fitted` fitted mortality rates using the additive multi-population mortality model.
#' * `logit.qxt.fitted` fitted mortality rates in logit way estimated with the additive multi-population mortality model.
#' * `qxt.future` future mortality rates estimated with the additive multi-population mortality model.
#' * `logit.qxt.future` future mortality rates in logit way estimated with the additive multi-population mortality model.
#' * `nPop` provided number of populations to fit the periods.
#'
#' @seealso \code{\link{fit.additive.LC.multi}}, \code{\link{fit.multiplicative.LC.multi}},
#' \code{\link{forecast.multiplicative.LC.multi}}, \code{\link{multipopulation_cv}},
#'
#'
#' @references
#' Debon, A., Montes, F., & Martiez-Ruiz, F. (2011).
#' Statistical methods to compare mortality for a group with non-divergent populations: an application to Spanish regions.
#' European Actuarial Journal, 1, 291-308.
#'
#' Lee, R.D. & Carter, L.R. (1992).
#' Modeling and forecasting US mortality.
#' Journal of the American Statistical Association, 87(419), 659–671.
#'
#' @examples
#' #First, we present the data that we are going to use
#' SpainRegions
#' ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)
#'
#' #Before forecast the future value of the mortality rates for different populations
#' #we need the object the fitted values of the additive multi-population mortality model.
#' additive_Spainmales <- fit.additive.LC.multi(qxt = SpainRegions$qx_male,
#'                                              periods = c(1991:2020),
#'                                              ages = c(ages),
#'                                              nPop = 18,
#'                                              lxt = SpainRegions$lx_male)
#' #Once, we have fit the data, it is possible to see the \eqn{\alpha_x}, \eqn{\beta_x}, \eqn{k_t}, and \eqn{I_i} provided parameters for the fitting.
#' plot.fit.additive.LC.multi(additive_Spainmales)
#'
#' #Finally, we forecast 10 years ahead the additive multi-population mortality model
#'
#' fut_additive_Spainmales <- forecast.additive.LC.multi(fitted.obj = additive_Spainmales, nahead = 10,
#'                                                       ktmethod = "Arimapdq", kt_include.cte = TRUE)
#' #As we mentioned in the details of the function, if we only provide the data from one-population the function
#' #\code{\link{fit.additive.LC.multi}} will fit the Lee-Carter model for single populations.
#' LC_Spainmales <- fit.additive.LC.multi(qxt = SpainNat$qx_male,
#'                               periods = c(1991:2020),
#'                               ages = ages,
#'                               nPop = 1)
#' plot.fit.LC.multi(LC_Spainmales)
#' #Again, we can forecast 10 years ahead using the LC mortality model for one-single population.
#' fut_LC_Spainmales <- forecast.additive.LC.multi(fitted.obj = LC_Spainmales,
#'         nahead = 10,ktmethod = "Arimapdq", kt_include.cte = TRUE)
#'
#' @export
forecast.additive.LC.multi <- function(fitted.obj, nahead,
                                       ktmethod = c("Arimapdq", "arima010"),
                                       kt_include.cte = TRUE){
  #Install library taht we need to forecast
  if (!require("forecast", character.only = TRUE)) {
    install.packages("forecast")
    library(forecast)
  } else {library(forecast)}

  #First check the structure of fitted.obj is equal to the previous object created using our function
  if(!identical(names(fitted.obj), c("ax", "bx", "kt", "Ii", "formula", "data.used",
                                     "qxt.real", "qxt.fitted", "logit.qxt.fitted",
                                     "Ages", "Periods","nPop"))){
    stop("The fitted.obj does not have the structure of R-library")
  }
  if(!is.list(fitted.obj)){
    stop("The fitted.obj is not a list. Use 'fit.additive.LC.multi' function first")
  }

  if (!is.numeric(nahead)) {
    stop("nahead has to be a numeric variable")
    if (nahead <= 0) {
      stop("nahead has to be higher than 0")
    }
  }
  #Construct the inv.logit -- function
  inv.logit <- function(x){exp(x)/(1+exp(x))}

  ktmethod <- match.arg(ktmethod)

  if(ktmethod == "Arimapdq"){
    kt.a2 <- auto.arima(fitted.obj$kt, ic = "bic",
                        allowdrift = kt_include.cte)
  } else if(ktmethod == "arima010"){
    kt.a2 <- Arima(fitted.obj$kt, order = c(0,1,0),
                   include.drift = kt_include.cte)
  }

  fut.kt <- forecast(kt.a2, h = nahead)

  kt.var <- list(kt.arima = kt.a2, mean.kt = fut.kt$mean,
                 lower = fut.kt$lower, upper = fut.kt$upper)

  lee.fut.logit <- list()
  lee.fut.qxt <- list()
  nages <- length(fitted.obj$ax)
  ages <- colnames(fitted.obj$ax)
  lastperiod <- as.numeric(rownames(fitted.obj$kt)[length(fitted.obj$kt)])+ nahead
  periods <- c((as.numeric(rownames(fitted.obj$kt)[length(fitted.obj$kt)])+1):lastperiod)

  #it <- 1
  for(it in 1:fitted.obj$nPop){
    lee.fut.logit[[paste0("pob", it)]] <- matrix(rep(fitted.obj$ax, nahead), nrow= nages, ncol = nahead)+
      (matrix(fitted.obj$bx, nrow=nages, ncol=1)%*%matrix(kt.var$mean.kt[1:nahead], nrow=1, ncol=nahead)) +
      matrix(rep(fitted.obj$Ii[it], nages*nahead), nrow= nages, ncol = nahead)

    lee.fut.qxt[[paste0("pob", it)]] <- inv.logit(lee.fut.logit[[it]])
    rownames(lee.fut.logit[[it]]) <- rownames(lee.fut.qxt[[it]]) <- ages
    colnames(lee.fut.logit[[it]]) <- colnames(lee.fut.qxt[[it]]) <- periods
  }

  return <- list(ax = matrix(fitted.obj$ax, nrow = 1, ncol = nages, dimnames = list("ax", ages)),
                 bx = matrix(fitted.obj$bx, nrow = 1, ncol = nages, dimnames = list("bx", ages)),
                 arimakt = kt.a2,
                 kt.fitted = matrix(fitted.obj$kt, nrow = length(fitted.obj$Periods), ncol = 1, dimnames = list(fitted.obj$Periods, "kt")),
                 kt.fut = matrix(kt.var$mean.kt[1:nahead], nrow= nahead, ncol=1,
                                 dimnames= list(c(max(fitted.obj$Periods+1):(max(fitted.obj$Periods)+nahead)),"kt")),
                 kt.futintervals = kt.var,
                 Ipop = matrix(fitted.obj$Ii, nrow = fitted.obj$nPop, ncol = 1, dimnames = list(c(1:fitted.obj$nPop), "Ipop")),
                 formula = fitted.obj$formula,
                 qxt.real = fitted.obj$qxt.real,
                 qxt.fitted = fitted.obj$qxt.fitted.qxt,
                 logit.qxt.fitted = fitted.obj$logit.qxt.fitted,
                 qxt.future = lee.fut.qxt,
                 logit.qxt.future = lee.fut.logit,
                 nPop = fitted.obj$nPop)
  return

}
