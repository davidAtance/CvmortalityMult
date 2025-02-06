#' Function to apply cross-validation techniques for testing the forecasting accuracy
#' of multi-population mortality models
#' @description
#' R function for testing the accuracy out-of-sample using different cross-validation techniques. The multi-population mortality models used by the package are: additive (Debon et al., 2011), multiplicative (Russolillo et al., 2011), common-factor (CFM) (Carter and Lee, 1992), joint-k (Carter and Lee, 2011), and augmented-common-factor (ACFM) (Li and Lee, 2005).
#' We provide a R function that employ the cross-validation techniques for three-way-array, following the preliminary idea for panel-time series, specifically for testing the forecasting ability of single mortality models (Atance et al. 2020).
#' These techniques consist on split the database in two parts: training set (to run the model) and test set (to check the forecasting accuracy of the model).
#' This procedure is repeated several times trying to check the forecasting accuracy in different ways.
#' With this function, the user can provide its own mortality rates for different populations and apply different cross-validation techniques.
#' The user must specify three main inputs in the function (`nahead`, `trainset1`, and `fixed_train_origin`) to apply a specific cross-validation technique between the different options.
#' Indeed, you can apply the next time-series cross-validation techniques, following the terminology employed by Bergmeir et al. (2012):
#'
#' 1. Fixed-Origin.
#' The technique chronologically splits the data set into two parts, first for training the model, and second for testing the forecasting accuracy.
#' This process predicts only once for different forecast horizons which are evaluated to assess the accuracy of the multi-population model, as can be seen in the next Figure.
#' {\figure{FIXED-ORIGIN.jpg}{options: width="100\%" alt="Figure: mai.png"}}
#'
#' The function "`multipopulation_cv()`" understands FIXED-ORIGIN when `trainset1` + `nahead` = number of provided periods and `fixed_train_origin` = `TRUE` (default value).
#' As an example, data set with periods from 1991 to 2020, `trainset1` = 25 and `nahead`= 5, with a total of 30, equals to length of the periods 1991:2020.
#'
#' 2. Rolling-Origin recalibration (RO-recalibration) evaluation.
#' In this technique, the data set is spitted into 'k' sub-sets of data, keeping chronologically order.
#' The first set of data corresponds to the training set where the model is fitted and the forecast are evaluated with a fixed horizon.
#' In every iteration, the model is enlarged and recalibrated adding the test-set periods (`nahead` in the function) to the training set and forecasting the next fixed horizon.
#' The idea is to keep the origin fixed and move the forecast origin in every iteration, as can be seen in the next Figure
#'
#' {\figure{RO_RECALIBRATION.jpg}{options: width="100\%" alt="Figure: mai.png"}}
#'
#' In the package, to apply this technique the users must provided a value of `trainset1` higher than two (to meet with the minimum time-series size), and `fixed_train_origin` = `TRUE` (default value), independently of the assigned value of `nahead`.
#' There are different resampling techniques that can be applied based on the values of `trainset1` and `nahead`.
#' Indeed, when `nahead` = 1 --- Leave-One-Out-Cross-Validation (LOOCV) with RO-recalibration will be applied.
#' Independently, of the number of periods in the first train set (`trainset1`).
#' When, `nahead` and `trainset1` are equal --- K-Fold-Cross-Validation (LOOCV) with RO-recalibration will be applied.
#' For the rest values of `nahead` and `trainset1` a standard time-series CV technique will be implemented.
#'
#' 3. Rolling-Window (RW) evaluation
#' The approach is very similar to the RO-recalibration, but maintaining the training set size constant at each forecast/iteration.
#' Maintaining the chronological order in each forecast, the training set adds the previous. projected periods of the test set and discards the earliest observations, as can be seen in the next Figure.
#'
#' {\figure{RW_RECALIBRATION}{options: width="100\%" alt="Figure: mai.png"}}
#'
#' To apply this technique, the `multipopulation_cv()` function requires that `fixed_train_origin` = c("`FALSE`", "`1`"), regardless of the values of `nahead` and `trainset1`.
#' Equally as in RO-recalibration, LOOCV, and k-fold can be applied with `nahead` = 1, or `nahead` equals to `trainset1`, respectively, but keeping the training set constant through the iterations.
#' Additionally, the common time-series CV approach can be applied for different values of `nahead` and `trainset1`.
#' When `fixed_train_origin` = "`FALSE`", at each iteration the training set adds the next `nahead` periods and discards the oldest keeping the training set size constant.
#' While `fixed_train_origin` = "`1`", at every iteration the training set only incorporates the next period ahead and discards only the latest period; maintaining the length of the training set constant and allowing to assess the forecasting accuracy of the mortality models in the long and medium term with different periods.
#'
#' It should be mentioned that this function is developed for cross-validation the forecasting accuracy of several populations.
#' However, in case you only consider one population, the function will forecast the Lee-Carter model for one population.
#' To test the forecasting accuracy of the selected model, the function provides five different measures: SSE, MSE, MAE, MAPE or All.
#' This measure of accuracy will be provided in different ways: a total measure, among ages considered, among populations and among projected blocked (periods).
#' Depending on how you want to check the forecasting accuracy of the model you could select one or other.
#' In this case, the measures will be obtained using the mortality rates in the normal scale as recommended by Santolino (2023) against the log scale.
#'
#' @param qxt mortality rates used to fit the multi-population mortality models. This rates can be provided in matrix or in data.frame.
#' @param model multi-population mortality model chosen to fit the mortality rates c("`additive`", "`multiplicative`", "`CFM`", "`joint-K`", "`ACFM`"). In case you do not provide any value, the function will apply the "`additive`" option.
#' @param periods number of years considered in the fitting in a vector way c(`minyear`:`maxyear`).
#' @param ages vector with the ages considered in the fitting. If the mortality rates provide from an abridged life tables, it is necessary to provide a vector with the ages, see the example.
#' @param nPop number of population considered for fitting.
#' @param lxt survivor function considered for every population, not necessary to provide.
#' @param nahead is a vector specifying the number of periods to forecast `nahead` periods ahead. It should be noted that when `nahead` is equal to `trainset1` a k-fold CV will be applied. Whereas when `nahead` is equal to 1, the CV process will be a Leave-One-Out CV.
#' @param trainset1 is a vector with the periods for the first training set.  This value must be greater than 2 to meet the minimum time series size (Hyndman and Khandakar, 2008).
#' @param fixed_train_origin option to select whether the origin in the first train set is fixed or not. The default value is `TRUE` where the origin of the first training sets is fixed. The alternatives are: `FALSE` when the first train set is moved in every iteration according to the provided `nahead` value, and 2. `1` when the train set is moved one period ahead in every repetition keeping constant the amount of data, and incorporating the next period observation, and discarding the last available period.
#' @param ktmethod method used to forecast the value of `kt` Arima(p,d,q) or ARIMA(0,1,0); c("`Arimapdq`", "`arima010`").
#' @param measures choose the non-penalized measure of forecasting accuracy that you want to use; c("`SSE`", "`MSE`", "`MAE`", "`MAPE`", "`All`"). Check the function. In case you do not provide any value, the function will apply the "`SSE`" as measure of forecasting accuracy.
#' @param ... other arguments for \code{\link[StMoMo:iarima]{iarima}}.
#'
#' @return An object of the class \code{"MultiCv"} including a `list()` with different components of the cross-validation process:
#' * `ax` parameter that captures the average shape of the mortality curve in all considered populations.
#' * `bx` parameter that explains the age effect x with respect to the general trend `kt` in the mortality rates of all considered populations.
#' * `kt.fitted` obtained values for the tendency behavior captured by `kt` .
#' * `kt.future` future values of `kt` for every iteration in the cross-validation.
#' * `kt.arima`  the arima selected for each `kt` time series.
#' * `Ii` parameter that captures the differences in the pattern of mortality in any region i with respect to Region 1.
#' * `formula` multi-population mortality formula used to fit the mortality rates.
#' * `model` provided the model selected in every case.
#' * `nPop` provided number of populations to fit the periods.
#' * `qxt.crude` corresponds to the crude mortality rates. These crude rate are directly obtained by dividing the number of registered deaths by the number of those initially exposed to the risk for age x, period t and in each region i.
#' * `qxt.future` future mortality rates estimated with the multi-population mortality model.
#' * `logit.qxt.future` future mortality rates in logit way estimated with the multi-population mortality model.
#' * `meas_ages` measure of forecasting accuracy through the ages of the study.
#' * `meas_periodsfut` measure of forecasting accuracy in every forecasting period(s) of the study.
#' * `meas_pop` measure of forecasting accuracy through the populations considered in the study.
#' * `meas_total` a global measure of forecasting accuracy through the ages, periods and populations of the study.
#' * `warn_msgs ` vector with the populations where the model has not converged.
#'
#' @seealso \code{\link{fitLCmulti}}, \code{\link{forecast.fitLCmulti}},
#' \code{\link{plot.fitLCmulti}}, \code{\link{plot.forLCmulti}},
#' \code{\link{MeasureAccuracy}}.
#'
#' @references
#' Atance, D., Debon, A., and Navarro, E. (2020).
#' A comparison of forecasting mortality models using resampling methods.
#' Mathematics 8(9): 1550.
#'
#' Bergmeir, C. & Benitez, J.M. (2012)
#' On the use of cross-validation for time series predictor evaluation.
#' Information Sciences, 191, 192–
#'
#' Carter, L.R. and Lee, R.D. (1992).
#' Modeling and forecasting US sex differentials in mortality.
#' International Journal of Forecasting, 8(3), 393–411.
#'
#' Debon, A., & Atance, D. (2022).
#' Two multi-population mortality models: A comparison of the forecasting accuracy with resampling methods.
#' in Contributions to Risk Analysis: Risk 2022. Fundacion Mapfre
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
#' Scandinavian Actuarial Journal, 96-117.
#'
#' Santolino, M. (2023).
#' Should Selection of the Optimum Stochastic Mortality Model Be Based on the Original or the Logarithmic Scale of the Mortality Rate?.
#' Risks, 11(10), 170.
#'
#' @importFrom gnm gnm residSVD Mult
#' @importFrom forecast Arima auto.arima forecast
#' @importFrom StMoMo genWeightMat
#' @importFrom utils install.packages
#' @importFrom stats plogis qlogis
#'
#' @examples
#'
#' #The example takes more than 5 seconds because they include
#' #several cross-validation methods and hence all the processes are included in "donttest".
#' \donttest{
#' #We present a cross-validation method for spanish male regions using:
# 'SpainRegions
#'
#' ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40,
#'          45, 50, 55, 60, 65, 70, 75, 80, 85, 90)
#' library(gnm)
#' library(forecast)
#' library(StMoMo)
#'
#' #1. FIXED-ORIGIN -- using the ACFM nahead + trainset1 = periods;
#' #fixed_train_origin = TRUE (defualt value)
#' ho_Spainmales_addit <- multipopulation_cv(qxt = SpainRegions$qx_male,
#'                                          model = c("ACFM"),
#'                                          periods =  c(1991:2020), ages = c(ages),
#'                                          nPop = 18, lxt = SpainRegions$lx_male,
#'                                          nahead = 5,
#'                                          trainset1 = 25,
#'                                          ktmethod = c("Arimapdq"),
#'                                          measures = c("SSE"))
#' ho_Spainmales_addit
#'
#' #Once, we have run the function we can check the result in different ways:
#' ho_Spainmales_addit$meas_ages
#' ho_Spainmales_addit$meas_periodsfut
#' ho_Spainmales_addit$meas_pop
#' ho_Spainmales_addit$meas_total
#'
#' #2. Let's continue with a RO-recalibration,
#' #(fixed_train_origin = TRUE (defualt value))
#' #where we have implemented three main CV techniques:
#' #2.1. Leave-One-Out-Cross-Validation (LOOCV) RO-recalibration when nahead = 1;
#' #(independently the number of periods blocked for the first train set; trainset1"
#' loocv_Spainmales_addit <- multipopulation_cv(qxt = SpainRegions$qx_male,
#'                                          model = c("additive"),
#'                                          periods =  c(1991:2020), ages = c(ages),
#'                                          nPop = 18, lxt = SpainRegions$lx_male,
#'                                          nahead = 1, trainset1 = 10,
#'                                          ktmethod = c("Arimapdq"),
#'                                          measures = c("SSE"))
#' loocv_Spainmales_addit
#'
#' #Once, we have run the function we can check the result in different ways:
#' loocv_Spainmales_addit$meas_ages
#' loocv_Spainmales_addit$meas_periodsfut
#' loocv_Spainmales_addit$meas_pop
#' loocv_Spainmales_addit$meas_total
#'
#' #2.2. K-Fold-CV RO-recalibration when nahead = trainset1
#' kfoldcv_Spainmales_addit <- multipopulation_cv(qxt = SpainRegions$qx_male,
#'                                          model = c("additive"),
#'                                          periods =  c(1991:2020), ages = c(ages),
#'                                          nPop = 18, lxt = SpainRegions$lx_male,
#'                                          nahead = 5, trainset1 = 5,
#'                                          ktmethod = c("Arimapdq"),
#'                                          measures = c("SSE"))
#' kfoldcv_Spainmales_addit
#'
#' #Once, we have run the function we can check the result in different ways:
#' kfoldcv_Spainmales_addit$meas_ages
#' kfoldcv_Spainmales_addit$meas_periodsfut
#' kfoldcv_Spainmales_addit$meas_pop
#' kfoldcv_Spainmales_addit$meas_total
#'
#' #2.3. standard time-series CV
#' cv_Spainmales_addit <- multipopulation_cv(qxt = SpainRegions$qx_male,
#'                                          model = c("additive"),
#'                                          periods =  c(1991:2020), ages = c(ages),
#'                                          nPop = 18, lxt = SpainRegions$lx_male,
#'                                          nahead = 5, trainset1 = 10,
#'                                          fixed_train_origin = "TRUE",
#'                                          ktmethod = c("Arimapdq"),
#'                                          measures = c("SSE"))
#'
#' cv_Spainmales_addit
#' #Once, we have run the function we can check the result in different ways:
#' cv_Spainmales_addit$meas_ages
#' cv_Spainmales_addit$meas_periodsfut
#' cv_Spainmales_addit$meas_pop
#' cv_Spainmales_addit$meas_total
#'
#' #3. RW-evaluation (fixed_train_origin = c("FALSE", "1"))
#' #3.1. fixed_train_origin == "TRUE" (The default value)
#' #In this case, the previous processes (Fixed-Origin or RO-recalibration)
#' #3.2. fixed_train_origin == "FALSE"
#' #where the origin in the training set is moved "nahead" period ahead in every iteration.
#' #This process allows to test the forecasting accuracy of "nahead" periods ahead
#' #keeping constant the size of the training and test set. As an example, we present
#' #three methods
#' #3.2.1. LOOCV
#' loocv_Spainmales_addit_rw <- multipopulation_cv(qxt = SpainRegions$qx_male,
#'                                          model = c("additive"),
#'                                          periods =  c(1991:2020), ages = c(ages),
#'                                          nPop = 18, lxt = SpainRegions$lx_male,
#'                                          nahead = 1, trainset1 = 10,
#'                                          fixed_train_origin = "FALSE",
#'                                          ktmethod = c("Arimapdq"),
#'                                          measures = c("SSE"))
#'
#' loocv_Spainmales_addit_rw
#'
#' #Once, we have run the function we can check the result in different ways:
#' loocv_Spainmales_addit_rw$meas_ages
#' loocv_Spainmales_addit_rw$meas_periodsfut
#' loocv_Spainmales_addit_rw$meas_pop
#' loocv_Spainmales_addit_rw$meas_total
#'
#' #3.2.2. K-Fold-CV
#' kfoldcv_Spainmales_addit_rw <- multipopulation_cv(qxt = SpainRegions$qx_male,
#'                                          model = c("additive"),
#'                                          periods =  c(1991:2020), ages = c(ages),
#'                                          nPop = 18, lxt = SpainRegions$lx_male,
#'                                          nahead = 5, trainset1 = 5,
#'                                          fixed_train_origin = "FALSE",
#'                                          ktmethod = c("Arimapdq"),
#'                                          measures = c("SSE"))
#'
#' kfoldcv_Spainmales_addit_rw
#'
#' #Once, we have run the function we can check the result in different ways:
#' kfoldcv_Spainmales_addit$meas_ages
#' kfoldcv_Spainmales_addit$meas_periodsfut
#' kfoldcv_Spainmales_addit$meas_pop
#' kfoldcv_Spainmales_addit$meas_total
#'
#' #3.2.3. standard time-series CV
#' cv_Spainmales_addit_rw <- multipopulation_cv(qxt = SpainRegions$qx_male,
#'                                          model = c("additive"),
#'                                          periods =  c(1991:2020), ages = c(ages),
#'                                          nPop = 18, lxt = SpainRegions$lx_male,
#'                                          nahead = 5, trainset1 = 10,
#'                                          fixed_train_origin = "FALSE",
#'                                          ktmethod = c("Arimapdq"),
#'                                          measures = c("SSE"))
#'
#' cv_Spainmales_addit_rw
#'
#' #Once, we have run the function we can check the result in different ways:
#' cv_Spainmales_addit_rw$meas_ages
#' cv_Spainmales_addit_rw$meas_periodsfut
#' cv_Spainmales_addit_rw$meas_pop
#' cv_Spainmales_addit_rw$meas_total
#'
#' #3.3  RW-evaluation (fixed_train_origin = c("1"))
#' #where the origin in the training set is moved 1 period ahead in every iteration.
#' #This process allows to test the forecasting accuracy of "nahead" periods ahead
#' #modifying the origin in the training set by 1.
#' #When "nahead" = 1 --- we will have a loocv equally as in the previous process,
#' #while using a different value of 1 for "nahead" we will test the forecasting
#' #accuracy of the model in "nahead" periods:
#' cv_Spainmales_addit_rw1 <- multipopulation_cv(qxt = SpainRegions$qx_male,
#'                                          model = c("additive"),
#'                                          periods =  c(1991:2020), ages = c(ages),
#'                                          nPop = 18, lxt = SpainRegions$lx_male,
#'                                          nahead = 5, trainset1 = 15,
#'                                          fixed_train_origin = "1",
#'                                          ktmethod = c("Arimapdq"),
#'                                          measures = c("SSE"))
#' cv_Spainmales_addit_rw1
#'
#' #Once, we have run the function we can check the result in different ways:
#' cv_Spainmales_addit_rw1$meas_ages
#' cv_Spainmales_addit_rw1$meas_periodsfut
#' cv_Spainmales_addit_rw1$meas_pop
#' cv_Spainmales_addit_rw1$meas_total
#'
#' }
#' @export
multipopulation_cv <- function(qxt, model = c("additive", "multiplicative", "CFM", "joint-K", "ACFM"),
                               periods, ages, nPop,
                               lxt=NULL,
                               ktmethod = c("Arimapdq", "arima010"),
                               nahead,
                               trainset1,
                               fixed_train_origin = TRUE,
                               measures = c("SSE", "MSE", "MAE", "MAPE", "All"),
                               ...){
  #Check several things before start
  if(is.null(qxt) || is.null(model) || is.null(periods) || is.null(ages) ||
     is.null(nPop) || is.null(ktmethod) || is.null(nahead) || is.null(trainset1) ||
     is.null(measures) ){
    stop("Arguments qxt, model, periods, ages, nPop, model, ktmethod, nahead, and trainset1, need to be provided.")
  }

  #2. Check that periods and ages are two vectors
  if (!is.vector(periods) || !is.vector(ages)) {
    stop("Period or Ages are not a vector, need to be provided.")
  }

  if(!is.numeric(nahead)){
    stop("nahead must be numeric variable.")
  }
  if(!is.numeric(trainset1)){
    stop("trainset1 must be numeric variable.")
  }

  valid_model <- c("additive", "multiplicative", "CFM","ACFM", "joint-K")
  model <- match.arg(model, valid_model)
  #In case, you donot provide any value of model, the function applies the additive one.

  valid_measures <- c("SSE", "MSE", "MAE", "MAPE", "All")
  measures <- match.arg(measures, valid_measures)

  nperiods <- length(periods)
  nages <- length(ages)

  #3. Check if lxt is provided
  if(is.null(lxt)){
    message("Argument lxt will be obtained as the number of individuals at age
    x and period t, starting with l0=100,000.")
  }

  #3. Check the structure of qxt

  #Starting with the values of qxt
  #1. Check if it is a list of matrix
  if(is.list(qxt)){
    message("Your qxt and lxt data are in list of matrix form.\n")

    #In the case is a list of matrix check that all matrix have the same length
    for(i in 2:length(qxt)){
      dim.actual <- dim(qxt[[i]])
      if(!identical(dim(qxt[[1]]), dim.actual)){
        stop("population ", i, " has a different dim regarding the rest of the chosen populations.")
      }
    }
    #Check the size of the qxt is equal to number of ages, periods and populations provided
    nper.age.pop <- nperiods*nages*nPop
    if(length(qxt[[i]])*length(qxt) != nper.age.pop){
      stop("Number of qxt is different from the period, ages and countries provided.")
    }

    #Check the size of list(qxt) is the same the npop
    if(length(qxt) != nPop){
      stop("Number of n Pop is different regarding the length of the matrix qxt.")
    }

    #Now, we transform the matrix to data.frame to fit the mortality data
    #create an empty data.frame
    df_qxtdata <- data.frame(matrix(NA, nrow = 0, ncol = 4))
    #Combine the matrixs into a data.frame
    for (i in 1:length(qxt)) {
      sec.per <- c()
      for(j in min(periods):max(periods)){
        sec.per <- c(sec.per, rep(j, nages))}
      message("Confirm the periods correspond to columns and ages to rows.\n")
      #Check how the matrix is formed (by ages/columns or by period/columns)

      #Estimate lx for every population, age and period
      if(is.null(lxt)){

        lxt.vector <- c()
        for(j in 1:nperiods){
          lx.prev <- vector("numeric", length(nages))
          lx.prev[1] <- 100000 #contiene lx de la edad cero
          for(nag in 1: (nages - 1)) {lx.prev[nag+1]<-lx.prev[nag]*(1-qxt.pop[(nag+(j-1)*nages)])}
          lx.vector <- c(lx.vector, lx.prev)
        }
      } else {
        if(dim(lxt[[i]])[1] == nages){
          lxt.vector <- as.vector(lxt[[i]])
        } else if (dim(qxt[[i]])[2] == nages) {
          lxt.vector <- as.vector(t(lxt[[i]]))
        }
      }
      #Create data.frame for estimate the model
      df_actual <- data.frame(rep((i-1), nages*periods),
                              sec.per,
                              rep(ages, nperiods),
                              qxt,
                              lxt.vector)

      df_qxtdata <- rbind(df_qxtdata, df_actual)
    }

    #for(i in 1:length(qxt)){
    # qxt.new <-
    #}
  }else if (is.vector(qxt)){

    #Check the size of the qxt is equal to number of ages, periods and populations provided
    nper.age.pop <- nperiods*nages*nPop
    if(length(qxt) != nper.age.pop){
      stop("Number of qxt vector is different from the period, ages and countries provided.")
    }
    message("Your qxt and lxt data are in vector form.\n")
    message("So, please ensure that qxt are provided by age, period and population, as follows:\n")
    message("age0-period0-pop1, age1-period0-pop1, ..., age.n-period0-pop1\n")
    message("age0-period1-pop1, age1-period1-pop1, ..., age.n-period1-pop1\n")
    message("...\n")
    message("age0-period0-pop2, age1-period0-pop2, ..., age.n-period0-pop2\n")

    df_qxtdata <- data.frame(matrix(NA, nrow = 0, ncol = 5))
    for (i in 1:nPop) {
      sec.per <- c()
      for(j in min(periods):max(periods)){
        sec.per <- c(sec.per, rep(j, nages))

      }
      qxt.pop <- qxt[(1+nages*nperiods*(i-1)):(nages*nperiods*(i))]

      #Estimate lx for every population, age and period
      if(is.null(lxt)){
        lxt.pop <- c()
        for(j in 1:nperiods){
          lx.prev <- vector("numeric", length(nages))
          lx.prev[1] <- 100000 #contiene lx de la edad cero
          for(nag in 1: (nages - 1)) {lx.prev[nag+1]<-lx.prev[nag]*(1-qxt.pop[(nag+(j-1)*nages)])}
          lxt.pop <- c(lxt.pop, lx.prev)
        }} else{lxt.pop <- lxt[(1+nages*nperiods*(i-1)):(nages*nperiods*(i))] }

      #Create data.frame for estimate the model
      df_actual <- data.frame(rep(c((i-1)), nages*nperiods),
                              sec.per,
                              rep(ages, nperiods),
                              qxt.pop,
                              lxt.pop)


      df_qxtdata <- rbind(df_qxtdata, df_actual)
    }

  } else{ message("qxt has to be provided as vector or list of matrix.")}

  #Give the correct name to the columns in the data.frame
  colnames(df_qxtdata) <- c("pop", "period", "age", "qxt", "lx")

  #Also, I am going to transform qxt into a matrix to provide as result
  mat_qxt <- list()
  for(i in 1:nPop){
    nperiods
    nages
    qxt1 <- df_qxtdata[df_qxtdata$pop == (i-1),]$qxt

    mat_qxt[[paste0("pop", i)]] <- matrix(qxt1, nrow = nages, ncol = nperiods,
                                          dimnames = list(ages, periods))
  }

  if(trainset1 > length(periods)){
    stop("number of periods in the trainset1 is higher than the total number of periods provided. Modify periods in the trainset1.")
  }

  if(nahead > length(periods)){
    stop("number of periods in the nahead is higher than the total number of periods provided. Modify periods in the nahead.")
  }

  if((trainset1+nahead) > length(periods)){
    stop("number of periods in the trainset1 and nahead is higher than the total number of periods provided. Modify periods in the nahead, or in the trainset1.")
  }

  if(trainset1 < 3){
    stop("number of trainset1 is lower than 3, minimum number of periods required to construct a time series.")
  }

  #We create several items to fill them
  fitting.data <- forecasting.data <- list()
  ax <- bx <- kt <- Ii <- kt.fut <- kt.arima <- warn <- list()
  CV <- NULL

  if(fixed_train_origin == TRUE){
    #Fixed-Origin --- Hold-Out
    if((nahead + trainset1) == nperiods){
      CV <- "Fixed-Origin (Hold-Out) procedure without rolling-window"
      nrep <- 1
      test1 <- nahead
      periods2 <- periods[1] + trainset1
      df_qxtdata2 <- df_qxtdata[df_qxtdata$period < periods2,]

      qxt.forecast <- list()
      for(i in 1:nPop){
        qxt.forecast[[paste0("pop", i)]] <- matrix(NA, nrow = nages, ncol = nahead,
                                                   dimnames = list(ages, c((periods[1] + trainset1):max(df_qxtdata$period))))

      }
      logit.qxt.forecast <- qxt.forecast

      if(model == "additive"){
        object <- fitLCmulti(model = "additive",
                             qxt = df_qxtdata2$qxt,
                             periods = c(min(periods):(periods2-1)),
                             ages = ages,
                             nPop = nPop, lxt = df_qxtdata2$lxt)
        forecast.obj <- forecast(object,
                                 nahead = nahead,
                                 ktmethod = ktmethod, ...)

      } else if(model == "multiplicative"){
        object <- fitLCmulti(model = "multiplicative",
                             qxt = df_qxtdata2$qxt,
                             periods = c(min(periods):(periods2-1)),
                             ages = ages,
                             nPop = nPop, lxt = df_qxtdata2$lxt)
        forecast.obj <- forecast(object,
                                 nahead = nahead,
                                 ktmethod = ktmethod, ...)

      } else if(model == "CFM"){
        object <- fitLCmulti(model = "CFM",
                             qxt = df_qxtdata2$qxt,
                             periods = c(min(periods):(periods2-1)),
                             ages = ages,
                             nPop = nPop, lxt = df_qxtdata2$lxt)
        forecast.obj <- forecast(object,
                                 nahead = nahead,
                                 ktmethod = ktmethod, ...)
      } else if(model == "ACFM"){
        object <- fitLCmulti(model = "ACFM",
                             qxt = df_qxtdata2$qxt,
                             periods = c(min(periods):(periods2-1)),
                             ages = ages,
                             nPop = nPop, lxt = df_qxtdata2$lxt)
        forecast.obj <- forecast(object,
                                 nahead = nahead,
                                 ktmethod = ktmethod, ...)

      } else if(model == "joint-K"){
        object <- fitLCmulti(model = "joint-K",
                             qxt = df_qxtdata2$qxt,
                             periods = c(min(periods):(periods2-1)),
                             ages = ages,
                             nPop = nPop, lxt = df_qxtdata2$lxt)
        forecast.obj <- forecast(object,
                                 nahead = nahead,
                                 ktmethod = ktmethod, ...)

      }

      for(i in 1:nPop){
        qxt.forecast[[i]][,1:nahead] <- forecast.obj$qxt.future[[i]]
      }

      fitting.data[[paste0("loop-1 from ", periods[1], " to ", periods2-1)]] <- object
      forecasting.data[[paste0("loop-1 from ", periods2, " to ", (periods2 + nahead - 1) )]] <- forecast.obj

      ax[[paste0("loop-1 from ", periods[1], " to ", periods2-1)]] <- object$ax
      bx[[paste0("loop-1 from ", periods[1], " to ", periods2-1)]] <- object$bx
      kt[[paste0("loop-1 from ", periods[1], " to ", periods2-1)]] <- object$kt
      Ii[[paste0("loop-1 from ", periods[1], " to ", periods2-1)]] <- object$Ii

      kt.fut[[paste0("loop-1 from ", periods[1], " to ", (periods2 - 1))]] <- forecast.obj$kt.fut
      kt.arima[[paste0("loop-1 from ", periods[1], " to ", (periods2 - 1))]] <- forecast.obj$arimakt

      if(model == "ACFM"){
        warn[[paste0("loop-1 from ", periods[1], " to ", periods2-1)]] <- object$warn_msgs
      }

    }else{
      if(nahead == 1){
        CV <- "Rolling-origin-recalibration evaluation using the LOOCV procedure keeping fixed the origin in the first train set."
      } else if(nahead == trainset1){
        CV <- "Rolling-origin-recalibration evaluation using the nahead-Fold-CV procedure keeping fixed the origin in the first train set."
      } else {
        CV <- "Rolling-origin-recalibration evaluation using the common-CV procedure keeping fixed the origin in the first train set."
      }
      #rolling-origin-recalibration
      #with different sizes for the first trainset.
      #indeed, loocv => nahead = 1, while k-fold => nahead = trainset1
      #Estimate the first set of fitting-sample
      fitting.years1 <- trainset1
      nrep <- (nperiods - fitting.years1)/nahead
      prue1 <- all.equal(nrep, as.integer(nrep))

      #we estimate the number of times (loops) that the model has to be fitted
      #and hence, the number of times to forecast the model
      if(isTRUE(prue1)){
        nrep <- as.integer(nrep)
        test1 <- c(rep(nahead, nrep))
      } else{
        nrep <- floor(nrep) + 1
        test1 <- c(rep(nahead, nrep - 1))
        test2 <- length(periods) - fitting.years1 - (nrep-1)*nahead
        test1 <- c(test1, test2)
      }
      qxt.forecast <- logit.qxt.forecast <- list()
      for(i in 1:nPop){
        qxt.forecast[[paste0("pop", i)]] <- matrix(NA, nrow = nages, ncol = sum(test1),
                                                   dimnames = list(ages, c((periods[1] + fitting.years1):max(df_qxtdata$period))))

      }


      for(reps in 1:nrep){
        if(reps == 1){
          fitting.years1 <- fitting.years1
        } else{
          fitting.years1 <- fitting.years1 + nahead
        }

        periods2 <- periods[1] + fitting.years1
        df_qxtdata2 <- df_qxtdata[df_qxtdata$period < periods2,]

        if(model == "additive"){
          object <- fitLCmulti(model = "additive",
                               qxt = df_qxtdata2$qxt,
                               periods = c(min(periods):(periods2-1)),
                               ages = ages,
                               nPop = nPop, lxt = df_qxtdata2$lxt)
          forecast.obj <- forecast(object,
                                   nahead = test1[reps],
                                   ktmethod = ktmethod, ...)

        } else if(model == "multiplicative"){
          object <- fitLCmulti(model = "multiplicative",
                               qxt = df_qxtdata2$qxt,
                               periods = c(min(periods):(periods2-1)),
                               ages = ages,
                               nPop = nPop, lxt = df_qxtdata2$lxt)
          forecast.obj <- forecast(object,
                                   nahead = test1[reps],
                                   ktmethod = ktmethod, ...)

        } else if(model == "CFM"){
          object <- fitLCmulti(model = "CFM",
                               qxt = df_qxtdata2$qxt,
                               periods = c(min(periods):(periods2-1)),
                               ages = ages,
                               nPop = nPop, lxt = df_qxtdata2$lxt)
          forecast.obj <- forecast(object,
                                   nahead = test1[reps],
                                   ktmethod = ktmethod, ...)

        } else if(model == "ACFM"){
          object <- fitLCmulti(model = "ACFM",
                               qxt = df_qxtdata2$qxt,
                               periods = c(min(periods):(periods2-1)),
                               ages = ages,
                               nPop = nPop, lxt = df_qxtdata2$lxt)
          forecast.obj <- forecast(object,
                                   nahead = test1[reps],
                                   ktmethod = ktmethod, ...)

        } else if(model == "joint-K"){
          object <- fitLCmulti(model = "joint-K",
                               qxt = df_qxtdata2$qxt,
                               periods = c(min(periods):(periods2-1)),
                               ages = ages,
                               nPop = nPop, lxt = df_qxtdata2$lxt)
          forecast.obj <- forecast(object,
                                   nahead = test1[reps],
                                   ktmethod = ktmethod, ...)

        }

        if(reps == 1){
          npos <- 1
          for(i in 1:nPop){
            qxt.forecast[[i]][,1:test1[reps]] <- forecast.obj$qxt.future[[i]]
          }

        } else {
          npos <- sum(test1[1:(reps-1)]) + 1
          npos2 <- sum(test1[1:(reps)])
          for(i in 1:nPop){
            qxt.forecast[[i]][,npos:npos2] <- forecast.obj$qxt.future[[i]]
          }}

        fitting.data[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- object
        forecasting.data[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- forecast.obj

        ax[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- object$ax
        bx[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- object$bx
        kt[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- object$kt
        Ii[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- object$Ii

        kt.fut[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- forecast.obj$kt.fut
        kt.arima[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- forecast.obj$arimakt

        if(model == "ACFM"){
          warn[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- object$warn_msgs
        }
      }
    }

    }else if(fixed_train_origin == "FALSE"){
      #when the user wants a rolling-window-evaluation and every
      if((nahead + trainset1) == nperiods){
        stop("Fixed-Origin or Hold-Out method can not apply this kind of cross-validation, only one training set and one test set.")
      }else{
        #rolling-origin-recalibration
        #with different sizes for the first trainset.
        #indeed, loocv => nahead = 1, while k-fold => nahead = trainset1
        if(nahead == 1){
          CV <- "Rolling-window evaluation using the LOOCV procedure, deleting and adding 1 period ahead to the train set."
        } else if(nahead == trainset1){
          CV <- "Rolling-window evaluation using the nahead-Fold-CV procedure, deleting and adding nahead-periods to the train set."
        } else {
          CV <- "Rolling-window evaluation using the common-CV procedure, deleting and adding nahead-periods to the train set."
        }

        #Estimate the first set of fitting-sample
        fitting.years1 <- trainset1
        nrep <- (nperiods - fitting.years1)/nahead
        prue1 <- all.equal(nrep, as.integer(nrep))

        #we estimate the number of times (loops) that the model has to be fitted
        #and hence, the number of times to forecast the model
        if(isTRUE(prue1)){
          nrep <- as.integer(nrep)
          test1 <- c(rep(nahead, nrep))
        } else{
          nrep <- floor(nrep) + 1
          test1 <- c(rep(nahead, nrep - 1))
          test2 <- length(periods) - fitting.years1 - (nrep-1)*nahead
          test1 <- c(test1, test2)
        }

        qxt.forecast <- logit.qxt.forecast <- list()
        for(i in 1:nPop){
          qxt.forecast[[paste0("pop", i)]] <- matrix(NA, nrow = nages, ncol = sum(test1),
                                                     dimnames = list(ages, c((periods[1] + fitting.years1):max(df_qxtdata$period))))

        }

        for(reps in 1:nrep){
          if(reps == 1){
            fitting.years1 <- fitting.years1
            periods3 <- periods[1] - 1
          } else{
            fitting.years1 <- fitting.years1 + nahead
            periods3 <- periods[1] - 1 + nahead*(reps-1)
          }

          periods2 <- periods[1] + fitting.years1
          df_qxtdata2 <- df_qxtdata[df_qxtdata$period>periods3 & df_qxtdata$period < periods2,]

          if(model == "additive"){
            object <- fitLCmulti(model = "additive",
                                 qxt = df_qxtdata2$qxt,
                                 periods = c(min(df_qxtdata2$period):(periods2-1)),
                                 ages = ages,
                                 nPop = nPop, lxt = df_qxtdata2$lxt)
            forecast.obj <- forecast(object,
                                     nahead = test1[reps],
                                     ktmethod = ktmethod, ...)

          } else if(model == "multiplicative"){
            object <- fitLCmulti(model = "multiplicative",
                                 qxt = df_qxtdata2$qxt,
                                 periods = c(min(df_qxtdata2$period):(periods2-1)),
                                 ages = ages,
                                 nPop = nPop, lxt = df_qxtdata2$lxt)
            forecast.obj <- forecast(object,
                                     nahead = test1[reps],
                                     ktmethod = ktmethod, ...)

          } else if(model == "CFM"){
            object <- fitLCmulti(model = "CFM",
                                 qxt = df_qxtdata2$qxt,
                                 periods = c(min(df_qxtdata2$period):(periods2-1)),
                                 ages = ages,
                                 nPop = nPop, lxt = df_qxtdata2$lxt)
            forecast.obj <- forecast(object,
                                     nahead = test1[reps],
                                     ktmethod = ktmethod, ...)

          } else if(model == "ACFM"){
            object <- fitLCmulti(model = "ACFM",
                                 qxt = df_qxtdata2$qxt,
                                 periods = c(min(df_qxtdata2$period):(periods2-1)),
                                 ages = ages,
                                 nPop = nPop, lxt = df_qxtdata2$lxt)
            forecast.obj <- forecast(object,
                                     nahead = test1[reps],
                                     ktmethod = ktmethod, ...)

          } else if(model == "joint-K"){
            object <- fitLCmulti(model = "joint-K",
                                 qxt = df_qxtdata2$qxt,
                                 periods = c(min(df_qxtdata2$period):(periods2-1)),
                                 ages = ages,
                                 nPop = nPop, lxt = df_qxtdata2$lxt)
            forecast.obj <- forecast(object,
                                     nahead = test1[reps],
                                     ktmethod = ktmethod, ...)

          }
          if(reps == 1){
            npos <- 1
            for(i in 1:nPop){
              qxt.forecast[[i]][,1:test1[reps]] <- forecast.obj$qxt.future[[i]]
            }

          } else {
            npos <- sum(test1[1:(reps-1)]) + 1
            npos2 <- sum(test1[1:(reps)])
            for(i in 1:nPop){
              qxt.forecast[[i]][,npos:npos2] <- forecast.obj$qxt.future[[i]]
            }}

          fitting.data[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- object
          forecasting.data[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- forecast.obj

          ax[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- object$ax
          bx[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- object$bx
          kt[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- object$kt
          Ii[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- object$Ii

          kt.fut[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- forecast.obj$kt.fut
          kt.arima[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- forecast.obj$arimakt

          if(model == "ACFM"){
            warn[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- object$warn_msgs
          }
      }

    }
    }else if(fixed_train_origin == "1"){
      CV <- "Rolling-window evaluation using the nahead-step-ahead forecast procedure, keeping the size of the test and train sets while add and discard one period in every iteration."
      #When the user wants to test the forecast accuracy of
      #nahead-step-ahead forecast with nahead != 1 keeping the size of the
      #test and train sets.
      fitting.years1 <- trainset1
      nrep <- (nperiods - fitting.years1 - nahead + 1)
      prue1 <- TRUE

      #we estimate the number of times (loops) that the model has to be fitted
      #and hence, the number of times to forecast the model
      if(isTRUE(prue1)){
        nrep <- as.integer(nrep)
        test1 <- c(rep(nahead, nrep))
      }

      qxt.forecast <- logit.qxt.forecast <- list()

      for(reps in 1:nrep){
        if(reps == 1){
          fitting.years1 <- fitting.years1
          periods3 <- periods[1] - 1
        } else{
          fitting.years1 <- fitting.years1 + 1
          periods3 <- periods[1] + reps - 2
        }

        periods2 <- periods[1] + fitting.years1
        df_qxtdata2 <- df_qxtdata[df_qxtdata$period>periods3 & df_qxtdata$period < periods2,]
        for(i in 1:nPop){
          qxt.forecast[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]][[paste0("pop", i)]] <-
            matrix(NA, nrow = nages, ncol = test1[reps],
                   dimnames = list(ages, c((periods2):(periods2 + nahead - 1))))
        }

        if(model == "additive"){
          object <- fitLCmulti(model = "additive",
                               qxt = df_qxtdata2$qxt,
                               periods = c(min(df_qxtdata2$period):(periods2-1)),
                               ages = ages,
                               nPop = nPop, lxt = df_qxtdata2$lxt)
          forecast.obj <- forecast(object,
                                   nahead = test1[reps],
                                   ktmethod = ktmethod, ...)

        } else if(model == "multiplicative"){
          object <- fitLCmulti(model = "multiplicative",
                               qxt = df_qxtdata2$qxt,
                               periods = c(min(df_qxtdata2$period):(periods2-1)),
                               ages = ages,
                               nPop = nPop, lxt = df_qxtdata2$lxt)
          forecast.obj <- forecast(object,
                                   nahead = test1[reps],
                                   ktmethod = ktmethod, ...)

        } else if(model == "CFM"){
          object <- fitLCmulti(model = "CFM",
                               qxt = df_qxtdata2$qxt,
                               periods = c(min(df_qxtdata2$period):(periods2-1)),
                               ages = ages,
                               nPop = nPop, lxt = df_qxtdata2$lxt)
          forecast.obj <- forecast(object,
                                   nahead = test1[reps],
                                   ktmethod = ktmethod, ...)

        } else if(model == "ACFM"){
          object <- fitLCmulti(model = "ACFM",
                               qxt = df_qxtdata2$qxt,
                               periods = c(min(df_qxtdata2$period):(periods2-1)),
                               ages = ages,
                               nPop = nPop, lxt = df_qxtdata2$lxt)
          forecast.obj <- forecast(object,
                                   nahead = test1[reps],
                                   ktmethod = ktmethod, ...)

        } else if(model == "joint-K"){
          object <- fitLCmulti(model = "joint-K",
                               qxt = df_qxtdata2$qxt,
                               periods = c(min(df_qxtdata2$period):(periods2-1)),
                               ages = ages,
                               nPop = nPop, lxt = df_qxtdata2$lxt)
          forecast.obj <- forecast(object,
                                   nahead = test1[reps],
                                   ktmethod = ktmethod, ...)

        }
        logit.qxt.forecast[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- list()
        for(i in 1:nPop){
            qxt.forecast[[reps]][[i]][,1:test1[reps]] <- forecast.obj$qxt.future[[i]]
            logit.qxt.forecast[[reps]][[paste0("pop", i)]] <-
              matrix(NA, nrow = nages, ncol = test1[reps],
                     dimnames = list(ages, c((periods2):(periods2 + nahead - 1))))

            logit.qxt.forecast[[reps]][[i]] <- qlogis(qxt.forecast[[reps]][[i]])
        }

        fitting.data[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- object
        forecasting.data[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- forecast.obj

        ax[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- object$ax
        bx[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- object$bx
        kt[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- object$kt
        Ii[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- object$Ii

        kt.fut[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- forecast.obj$kt.fut
        kt.arima[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- forecast.obj$arimakt

        if(model == "ACFM"){
          warn[[paste0("loop-", reps, " from ", (periods3+1), " to ", periods2-1)]] <- object$warn_msgs
        }
      }
    }

  #We estimate the logit for the forecast qxt
  if(fixed_train_origin != '1'){
    for(i in 1:nPop){
    logit.qxt.forecast[[i]] <- qlogis(qxt.forecast[[i]])
    }
    if(measures == "SSE"){
      #We estimate SSE in different options
      meas_prevpops <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("SSE", c(1:nPop)))

      meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), "SSE"))
      meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("SSE", "all ages and periods"))

      meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "SSE"))
      meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("SSE", ages))

      meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_periods <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("SSE", c(1:length(test1))))

      wxt <- genWeightMat(ages = ages, years = c((periods[1] + trainset1):max(df_qxtdata$period)), clip = 0)
      for(i in 1:nPop){
        for(j in 1:(length(test1))){
          if(j == 1){
            npos <- 1
            meas_prevpops[j,i] <- MeasureAccuracy(measure = "SSE",
                                                  qxt_crude = mat_qxt[[i]][,(trainset1+nahead*(j-1)+1):(trainset1+sum(test1[1:j]))],
                                                  qxt_aju = qxt.forecast[[i]][,(1:nahead)],
                                                  wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))$value
            prev <- (mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead[1])])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevtotal[j,i] <- sum(prev)

          } else {
            npos <- sum(test1[1:(j-1)]) + 1
            npos2 <- sum(test1[1:(j)])
            nrw <-
              meas_prevpops[j,i] <- MeasureAccuracy(measure = "SSE",
                                                    qxt_crude = mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)],
                                                    qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                                    wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))$value
            prev <- (mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevtotal[j,i] <- sum(prev)
          }
        }
        meas_pops[,i] <- mean(meas_prevpops[,i])
      }
      for(j in 1:(length(test1))){
        meas_prevtot[j,] <- sum(meas_prevtotal[j,])
      }
      meas_total <- mean(meas_prevtot)

      for(j in 1:nages){
        for(pe in 1:length(test1)){
          for(i in 1:nPop){
            if(pe == 1){
              npos <- 1
              prev <- (mat_qxt[[i]][j,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][j,(1:nahead)])^2
              prev <- replace(prev, prev == "Inf", 0)
              prev <- replace(prev, prev == "-Inf", 0)
              prev <- replace(prev, prev == "NA", 0)
              prev <- replace(prev, prev == "NaN", 0)
              meas_prevag[i,] <- sum(prev)

            } else {
              npos <- sum(test1[1:(pe-1)]) + 1
              npos2 <- sum(test1[1:(pe)])
              prev <- (mat_qxt[[i]][j,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])^2
              prev <- replace(prev, prev == "Inf", 0)
              prev <- replace(prev, prev == "-Inf", 0)
              prev <- replace(prev, prev == "NA", 0)
              prev <- replace(prev, prev == "NaN", 0)
              meas_prevag[i,] <- sum(prev)
            }
          }
          meas_prevages[pe,j] <- sum(meas_prevag)
        }
        meas_ages[,j] <- mean(meas_prevages[,j])
      }

      for(j in 1:length(test1)){
        for(i in 1:nPop){
          if(j == 1){
            npos <- 1
            prev <- (mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)

          } else {
            npos <- sum(test1[1:(j-1)]) + 1
            npos2 <- sum(test1[1:(j)])
            prev <- (mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
          }
        }
        meas_periods[,j] <- sum(meas_prevperiods[,j])
      }
    }else if(measures == "MSE"){
      #We estimate MSE in different options
      meas_prevpops <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("MSE", c(1:nPop)))

      meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), "MSE"))
      meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("MSE", "all ages and periods"))

      meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "MSE"))
      meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("MSE", ages))

      meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_periods <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("MSE", c(1:length(test1))))

      wxt <- genWeightMat(ages = ages, years = c((periods[1] + trainset1):max(df_qxtdata$period)), clip = 0)
      for(i in 1:nPop){
        for(j in 1:(length(test1))){
          if(j == 1){
            npos <- 1
            meas_prevpops[j,i] <- MeasureAccuracy(measure = "MSE",
                                                  qxt_crude = mat_qxt[[i]][,(trainset1+nahead*(j-1)+1):(trainset1+sum(test1[1:j]))],
                                                  qxt_aju = qxt.forecast[[i]][,(1:nahead)],
                                                  wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))$value
            prev <- (mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevtotal[,i] <- sum(prev)
          } else {
            npos <- sum(test1[1:(j-1)]) + 1
            npos2 <- sum(test1[1:(j)])
            meas_prevpops[j,i] <- MeasureAccuracy(measure = "MSE",
                                                  qxt_crude = mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)],
                                                  qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                                  wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))$value
            prev <- (mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevtotal[j,i] <- sum(prev)
          }
        }
        meas_pops[,i] <- mean(meas_prevpops[,i])
      }
      for(j in 1:(length(test1))){
        meas_prevtot[j,] <- sum(meas_prevtotal[j,])/(nages*test1[j]*nPop)
      }
      meas_total <- mean(meas_prevtot)

      for(j in 1:nages){
        for(pe in 1:length(test1)){
          for(i in 1:nPop){
            if(pe == 1){
              npos <- 1
              prev <- (mat_qxt[[i]][j,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][j,(1:nahead)])^2
              prev <- replace(prev, prev == "Inf", 0)
              prev <- replace(prev, prev == "-Inf", 0)
              prev <- replace(prev, prev == "NA", 0)
              prev <- replace(prev, prev == "NaN", 0)
              meas_prevag[i,] <- sum(prev)

            } else {
              npos <- sum(test1[1:(pe-1)]) + 1
              npos2 <- sum(test1[1:(pe)])
              prev <- (mat_qxt[[i]][j,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])^2
              prev <- replace(prev, prev == "Inf", 0)
              prev <- replace(prev, prev == "-Inf", 0)
              prev <- replace(prev, prev == "NA", 0)
              prev <- replace(prev, prev == "NaN", 0)
              meas_prevag[i,] <- sum(prev)
            }
          }
          meas_prevages[pe,j] <- sum(meas_prevag)/(test1[pe]*nPop)
        }
        meas_ages[,j] <- mean(meas_prevages[,j])
      }

      for(j in 1:length(test1)){
        for(i in 1:nPop){
          if(j == 1){
            npos <- 1
            prev <- (mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
          } else {
            npos <- sum(test1[1:(j-1)]) + 1
            npos2 <- sum(test1[1:(j)])
            prev <- (mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
          }
        }
        meas_periods[,j] <- sum(meas_prevperiods[,j])/(nages*nPop*test1[j])
      }
    }else if(measures == "MAE"){
      #We estimate MAE in different options
      meas_prevpops <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("MAE", c(1:nPop)))

      meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), "MAE"))
      meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("MAE", "all ages and periods"))

      meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "MAE"))
      meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("MAE", ages))

      meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_periods <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("MAE", c(1:length(test1))))

      wxt <- genWeightMat(ages = ages, years = c((periods[1] + trainset1):max(df_qxtdata$period)), clip = 0)

      for(i in 1:nPop){
        for(j in 1:(length(test1))){
          if(j == 1){
            npos <- 1
            meas_prevpops[j,i] <- MeasureAccuracy(measure = "MAE",
                                                  qxt_crude = mat_qxt[[i]][,(trainset1+nahead*(j-1)+1):(trainset1 + sum(test1[1:j]))],
                                                  qxt_aju = qxt.forecast[[i]][,(1:nahead)],
                                                  wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))$value
            prev <- abs(mat_qxt[[i]][,(trainset1+1):(trainset1 + nahead)] - qxt.forecast[[i]][,(1:nahead)])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevtotal[j,i] <- sum(prev)
          } else {
            npos <- sum(test1[1:(j-1)]) + 1
            npos2 <- sum(test1[1:(j)])
            meas_prevpops[j,i] <- MeasureAccuracy(measure = "MAE",
                                                  qxt_crude = mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)],
                                                  qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                                  wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))$value
            prev <- abs(mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevtotal[j,i] <- sum(prev)
          }
        }
        meas_pops[,i] <- mean(meas_prevpops[,i])
      }
      for(j in 1:(length(test1))){
        meas_prevtot[j,] <- sum(meas_prevtotal[j,])/(nPop*nages*test1[j])
      }
      meas_total <- mean(meas_prevtot)

      for(j in 1:nages){
        for(pe in 1:length(test1)){
          for(i in 1:nPop){
            if(pe == 1){
              npos <- 1
              prev <- abs(mat_qxt[[i]][j,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][j,(1:nahead)])
              prev <- replace(prev, prev == "Inf", 0)
              prev <- replace(prev, prev == "-Inf", 0)
              prev <- replace(prev, prev == "NA", 0)
              prev <- replace(prev, prev == "NaN", 0)
              meas_prevag[i,] <- sum(prev)
            } else {
              npos <- sum(test1[1:(pe-1)]) + 1
              npos2 <- sum(test1[1:(pe)])
              prev <- abs(mat_qxt[[i]][j,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])
              prev <- replace(prev, prev == "Inf", 0)
              prev <- replace(prev, prev == "-Inf", 0)
              prev <- replace(prev, prev == "NA", 0)
              prev <- replace(prev, prev == "NaN", 0)
              meas_prevag[i,] <- sum(prev)
            }
          }
          meas_prevages[pe,j] <- sum(meas_prevag)/(test1[pe]*nPop)
        }
        meas_ages[,j] <- mean(meas_prevages[,j])
      }

      for(j in 1:length(test1)){
        for(i in 1:nPop){
          if(j == 1){
            npos <- 1
            prev <- abs(mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
          } else {
            npos <- sum(test1[1:(j-1)]) + 1
            npos2 <- sum(test1[1:(j)])
            prev <- abs(mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
          }
        }
        meas_periods[,j] <- sum(meas_prevperiods[,j])/(nages*nPop*test1[j])
      }

    }else if(measures == "MAPE"){
      #We estimate MAPE in different options
      meas_prevpops <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("MAPE", c(1:nPop)))

      meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), " MAPE"))
      meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("MAPE", "all ages and periods"))

      meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "mape"))
      meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("MAPE", ages))

      meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_periods <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("MAPE", c(1:length(test1))))

      wxt <- genWeightMat(ages = ages, years = c((periods[1] + trainset1):max(df_qxtdata$period)), clip = 0)
      for(i in 1:nPop){
        for(j in 1:(length(test1))){
          if(j == 1){
            npos <- 1
            meas_prevpops[j,i] <- MeasureAccuracy(measure = "MAPE",
                                                  qxt_crude = mat_qxt[[i]][,(trainset1+nahead*(j-1)+1):(trainset1+sum(test1[1:j]))],
                                                  qxt_aju = qxt.forecast[[i]][,(1:nahead)],
                                                  wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))$value
            prev <- abs((mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])/mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevtotal[j,i] <- sum(prev)
          } else {
            npos <- sum(test1[1:(j-1)]) + 1
            npos2 <- sum(test1[1:(j)])
            meas_prevpops[j,i] <- MeasureAccuracy(measure = "MAPE",
                                                  qxt_crude = mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)],
                                                  qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                                  wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))$value
            prev <- abs((mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])/mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevtotal[j,i] <- sum(prev)
          }
        }
        meas_pops[,i] <- mean(meas_prevpops[,i])
      }
      for(j in 1:(length(test1))){
        meas_prevtot[j,] <- sum(meas_prevtotal[j,])/(nPop*nages*test1[j])
      }
      meas_total <- mean(meas_prevtot)

      for(j in 1:nages){
        for(pe in 1:length(test1)){
          for(i in 1:nPop){
            if(pe == 1){
              npos <- 1
              prev <- (abs((mat_qxt[[i]][j,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][j,(1:nahead)])/mat_qxt[[i]][j,(trainset1+1):(trainset1+nahead)]))
              prev <- replace(prev, prev == "Inf", 0)
              prev <- replace(prev, prev == "-Inf", 0)
              prev <- replace(prev, prev == "NA", 0)
              prev <- replace(prev, prev == "NaN", 0)
              meas_prevag[i,] <- sum(prev)

            } else {
              npos <- sum(test1[1:(pe-1)]) + 1
              npos2 <- sum(test1[1:(pe)])
              prev <- (abs((mat_qxt[[i]][j,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])/mat_qxt[[i]][j,(trainset1+npos):(trainset1+npos2)]))
              prev <- replace(prev, prev == "Inf", 0)
              prev <- replace(prev, prev == "-Inf", 0)
              prev <- replace(prev, prev == "NA", 0)
              prev <- replace(prev, prev == "NaN", 0)
              meas_prevag[i,] <- sum(prev)
            }
          }
          meas_prevages[pe,j] <- sum(meas_prevag)/(test1[pe]*nPop)
        }
        meas_ages[,j] <- mean(meas_prevages[,j])
      }

      for(j in 1:length(test1)){
        for(i in 1:nPop){
          if(j == 1){
            npos <- 1
            prev <- abs((mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])/mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)

          } else {
            npos <- sum(test1[1:(j-1)]) + 1
            npos2 <- sum(test1[1:(j)])
            prev <- abs((mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])/mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
          }
        }
        meas_periods[,j] <- sum(meas_prevperiods[,j])/(nages*nPop*test1[j])
      }

    }else if(measures == "All"){
      #We estimate ALL measures in different options
      meas_prevpops1 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevpops2 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevpops3 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevpops4 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_pops <- matrix(NA, nrow = 4, ncol=nPop, dimnames = list(c("SSE", "MSE", "MAE", "MAPE"), c(1:nPop)))

      meas_prevtotal1 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtotal2 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtotal3 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtotal4 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))

      meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 4, dimnames = list(c(1:length(test1)),c("SSE", "MSE", "MAE", "MAPE")))
      meas_total <- matrix(NA, nrow = 4, ncol = 1, dimnames = list(c("SSE", "MSE", "MAE", "MAPE"), "all ages and periods"))

      meas_prevages1 <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevages2 <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevages3 <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevages4 <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevag <- matrix(NA, nrow = nPop, ncol = 4, dimnames = list(c(1:nPop), c("SSE", "MSE", "MAE", "MAPE")))
      meas_ages <- matrix(NA, nrow = 4, ncol = nages, dimnames = list(c("SSE", "MSE", "MAE", "MAPE"), ages))

      meas_prevperiods1 <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_prevperiods2 <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_prevperiods3 <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_prevperiods4 <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))

      meas_periods <- matrix(NA, nrow = 4, ncol = length(test1), dimnames = list(c("SSE", "MSE", "MAE", "MAPE"), c(1:length(test1))))
      wxt <- genWeightMat(ages = ages, years = c((periods[1] + trainset1):max(df_qxtdata$period)), clip = 0)
      for(i in 1:nPop){
        for(j in 1:(length(test1))){
          if(j == 1){
            npos <- 1
            meas_prevpops1[j,i] <- MeasureAccuracy(measure = "SSE", qxt_crude = mat_qxt[[i]][,(trainset1+nahead*(j-1)+1):(trainset1+sum(test1[1:j]))],
                                                   qxt_aju = qxt.forecast[[i]][,(1:nahead)],
                                                   wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))$value
            meas_prevpops2[j,i] <- MeasureAccuracy(measure = "MSE", qxt_crude = mat_qxt[[i]][,(trainset1+nahead*(j-1)+1):(trainset1+sum(test1[1:j]))],
                                                   qxt_aju = qxt.forecast[[i]][,(1:nahead)],
                                                   wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))$value
            meas_prevpops3[j,i] <- MeasureAccuracy(measure = "MAE", qxt_crude = mat_qxt[[i]][,(trainset1+nahead*(j-1)+1):(trainset1+sum(test1[1:j]))],
                                                   qxt_aju = qxt.forecast[[i]][,(1:nahead)],
                                                   wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))$value
            meas_prevpops4[j,i] <- MeasureAccuracy(measure = "MAPE", qxt_crude = mat_qxt[[i]][,(trainset1+nahead*(j-1)+1):(trainset1+sum(test1[1:j]))],
                                                   qxt_aju = qxt.forecast[[i]][,(1:nahead)],
                                                   wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))$value
            #SSE
            prev1 <- (mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])^2
            prev1 <- replace(prev1, prev1 == "Inf", 0)
            prev1 <- replace(prev1, prev1 == "-Inf", 0)
            prev1 <- replace(prev1, prev1 == "NA", 0)
            prev1 <- replace(prev1, prev1 == "NaN", 0)
            meas_prevtotal1[j,i] <- sum(prev1)
            #MSE
            prev2 <- (mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])^2
            prev2 <- replace(prev2, prev2 == "Inf", 0)
            prev2 <- replace(prev2, prev2 == "-Inf", 0)
            prev2 <- replace(prev2, prev2 == "NA", 0)
            prev2 <- replace(prev2, prev2 == "NaN", 0)
            meas_prevtotal2[j,i] <- sum(prev2)
            #MAE
            prev3 <- abs(mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])
            prev3 <- replace(prev3, prev3 == "Inf", 0)
            prev3 <- replace(prev3, prev3 == "-Inf", 0)
            prev3 <- replace(prev3, prev3 == "NA", 0)
            prev3 <- replace(prev3, prev3 == "NaN", 0)
            meas_prevtotal3[j,i] <- sum(prev3)
            #MAPE
            prev4 <- abs((mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])/mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)])
            prev4 <- replace(prev4, prev4 == "Inf", 0)
            prev4 <- replace(prev4, prev4 == "-Inf", 0)
            prev4 <- replace(prev4, prev4 == "NA", 0)
            prev4 <- replace(prev4, prev4 == "NaN", 0)
            meas_prevtotal4[j,i] <- sum(prev4)
          } else {
            npos <- sum(test1[1:(j-1)]) + 1
            npos2 <- sum(test1[1:(j)])
            meas_prevpops1[j,i] <- MeasureAccuracy(measure = "SSE", qxt_crude = mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)],
                                                   qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                                   wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))$value
            meas_prevpops2[j,i] <- MeasureAccuracy(measure = "MSE", qxt_crude = mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)],
                                                   qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                                   wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))$value
            meas_prevpops3[j,i] <- MeasureAccuracy(measure = "MAE", qxt_crude = mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)],
                                                   qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                                   wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))$value
            meas_prevpops4[j,i] <- MeasureAccuracy(measure = "MAPE", qxt_crude = mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)],
                                                   qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                                   wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))$value
            #SSE
            prev1 <- (mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
            prev1 <- replace(prev1, prev1 == "Inf", 0)
            prev1 <- replace(prev1, prev1 == "-Inf", 0)
            prev1 <- replace(prev1, prev1 == "NA", 0)
            prev1 <- replace(prev1, prev1 == "NaN", 0)
            meas_prevtotal1[j,i] <- sum(prev1)
            #MSE
            prev2 <- (mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
            prev2 <- replace(prev2, prev2 == "Inf", 0)
            prev2 <- replace(prev2, prev2 == "-Inf", 0)
            prev2 <- replace(prev2, prev2 == "NA", 0)
            prev2 <- replace(prev2, prev2 == "NaN", 0)
            meas_prevtotal2[j,i] <- sum(prev2)
            #MAE
            prev3 <- abs(mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])
            prev3 <- replace(prev3, prev3 == "Inf", 0)
            prev3 <- replace(prev3, prev3 == "-Inf", 0)
            prev3 <- replace(prev3, prev3 == "NA", 0)
            prev3 <- replace(prev3, prev3 == "NaN", 0)
            meas_prevtotal3[j,i] <- sum(prev3)
            #MAPE
            prev4 <- abs((mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])/mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)])
            prev4 <- replace(prev4, prev4 == "Inf", 0)
            prev4 <- replace(prev4, prev4 == "-Inf", 0)
            prev4 <- replace(prev4, prev4 == "NA", 0)
            prev4 <- replace(prev4, prev4 == "NaN", 0)
            meas_prevtotal4[j,i] <- sum(prev4)
          }
        }
        meas_pops[1,i] <- mean(meas_prevpops1[,i])
        meas_pops[2,i] <- mean(meas_prevpops2[,i])
        meas_pops[3,i] <- mean(meas_prevpops3[,i])
        meas_pops[4,i] <- mean(meas_prevpops4[,i])
      }
      for(j in 1:(length(test1))){
        meas_prevtot[j,1] <- sum(meas_prevtotal1[j,])
        meas_prevtot[j,2] <- sum(meas_prevtotal2[j,])/(nPop*nages*test1[j])
        meas_prevtot[j,3] <- sum(meas_prevtotal3[j,])/(nPop*nages*test1[j])
        meas_prevtot[j,4] <- sum(meas_prevtotal4[j,])/(nPop*nages*test1[j])
      }
      meas_total[1,1] <- mean(meas_prevtot[,1])
      meas_total[2,1] <- mean(meas_prevtot[,2])
      meas_total[3,1] <- mean(meas_prevtot[,3])
      meas_total[4,1] <- mean(meas_prevtot[,4])

      for(j in 1:nages){
        for(pe in 1:length(test1)){
          for(i in 1:nPop){
            if(pe == 1){
              npos <- 1
              #SSE
              prev <- (mat_qxt[[i]][j,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][j,(1:nahead)])^2
              prev1 <- replace(prev1, prev1 == "Inf", 0)
              prev1 <- replace(prev1, prev1 == "-Inf", 0)
              prev1 <- replace(prev1, prev1 == "NA", 0)
              prev1 <- replace(prev1, prev1 == "NaN", 0)
              meas_prevag[i,1] <- sum(prev1)

              #MSE
              prev2 <- (mat_qxt[[i]][j,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][j,(1:nahead)])^2
              prev2 <- replace(prev2, prev2 == "Inf", 0)
              prev2 <- replace(prev2, prev2 == "-Inf", 0)
              prev2 <- replace(prev2, prev2 == "NA", 0)
              prev2 <- replace(prev2, prev2 == "NaN", 0)
              meas_prevag[i,2] <- sum(prev2)

              #MAE
              prev3 <- abs(mat_qxt[[i]][j,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][j,(1:nahead)])
              prev3 <- replace(prev3, prev3 == "Inf", 0)
              prev3 <- replace(prev3, prev3 == "-Inf", 0)
              prev3 <- replace(prev3, prev3 == "NA", 0)
              prev3 <- replace(prev3, prev3 == "NaN", 0)
              meas_prevag[i,3] <- sum(prev3)
              #MAPE
              prev4 <- abs((mat_qxt[[i]][j,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][j,(1:nahead)])/mat_qxt[[i]][j,(trainset1+1):(trainset1+nahead)])
              prev4 <- replace(prev4, prev4 == "Inf", 0)
              prev4 <- replace(prev4, prev4 == "-Inf", 0)
              prev4 <- replace(prev4, prev4 == "NA", 0)
              prev4 <- replace(prev4, prev4 == "NaN", 0)
              meas_prevag[i,4] <- sum(prev4)
            } else {
              npos <- sum(test1[1:(pe-1)]) + 1
              npos2 <- sum(test1[1:(pe)])
              prev <- (mat_qxt[[i]][j,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])^2
              #SSE
              prev <- (mat_qxt[[i]][j,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])^2
              prev1 <- replace(prev1, prev1 == "Inf", 0)
              prev1 <- replace(prev1, prev1 == "-Inf", 0)
              prev1 <- replace(prev1, prev1 == "NA", 0)
              prev1 <- replace(prev1, prev1 == "NaN", 0)
              meas_prevag[i,1] <- sum(prev1)

              #MSE
              prev2 <- (mat_qxt[[i]][j,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])^2
              prev2 <- replace(prev2, prev2 == "Inf", 0)
              prev2 <- replace(prev2, prev2 == "-Inf", 0)
              prev2 <- replace(prev2, prev2 == "NA", 0)
              prev2 <- replace(prev2, prev2 == "NaN", 0)
              meas_prevag[i,2] <- sum(prev2)

              #MAE
              prev3 <- abs(mat_qxt[[i]][j,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])
              prev3 <- replace(prev3, prev3 == "Inf", 0)
              prev3 <- replace(prev3, prev3 == "-Inf", 0)
              prev3 <- replace(prev3, prev3 == "NA", 0)
              prev3 <- replace(prev3, prev3 == "NaN", 0)
              meas_prevag[i,3] <- sum(prev3)
              #MAPE
              prev4 <- abs((mat_qxt[[i]][j,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])/mat_qxt[[i]][j,(trainset1+npos):(trainset1+npos2)])
              prev4 <- replace(prev4, prev4 == "Inf", 0)
              prev4 <- replace(prev4, prev4 == "-Inf", 0)
              prev4 <- replace(prev4, prev4 == "NA", 0)
              prev4 <- replace(prev4, prev4 == "NaN", 0)
              meas_prevag[i,4] <- sum(prev4)
            }
          }
          meas_prevages1[pe,j] <- sum(meas_prevag[,1])
          meas_prevages2[pe,j] <- sum(meas_prevag[,2])/(test1[pe]*nPop)
          meas_prevages3[pe,j] <- sum(meas_prevag[,3])/(test1[pe]*nPop)
          meas_prevages4[pe,j] <- sum(meas_prevag[,4])/(test1[pe]*nPop)
        }
        meas_ages[1,j] <- mean(meas_prevages1[,j])
        meas_ages[2,j] <- mean(meas_prevages2[,j])
        meas_ages[3,j] <- mean(meas_prevages3[,j])
        meas_ages[4,j] <- mean(meas_prevages4[,j])
      }

      for(j in 1:length(test1)){
        for(i in 1:nPop){
          if(j == 1){
            npos <- 1
            #SSE
            prev1 <- (mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])^2
            prev1 <- replace(prev1, prev1 == "Inf", 0)
            prev1 <- replace(prev1, prev1 == "-Inf", 0)
            prev1 <- replace(prev1, prev1 == "NA", 0)
            prev1 <- replace(prev1, prev1 == "NaN", 0)
            meas_prevperiods1[i,j] <- sum(prev1)
            #MSE
            prev2 <- (mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])^2
            prev2 <- replace(prev2, prev2 == "Inf", 0)
            prev2 <- replace(prev2, prev2 == "-Inf", 0)
            prev2 <- replace(prev2, prev2 == "NA", 0)
            prev2 <- replace(prev2, prev2 == "NaN", 0)
            meas_prevperiods2[i,j] <- sum(prev2)
            #MAE
            prev3 <- abs(mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])
            prev3 <- replace(prev3, prev3 == "Inf", 0)
            prev3 <- replace(prev3, prev3 == "-Inf", 0)
            prev3 <- replace(prev3, prev3 == "NA", 0)
            prev3 <- replace(prev3, prev3 == "NaN", 0)
            meas_prevperiods3[i,j] <- sum(prev3)
            #MAPE
            prev4 <- abs((mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)] - qxt.forecast[[i]][,(1:nahead)])/mat_qxt[[i]][,(trainset1+1):(trainset1+nahead)])
            prev4 <- replace(prev4, prev4 == "Inf", 0)
            prev4 <- replace(prev4, prev4 == "-Inf", 0)
            prev4 <- replace(prev4, prev4 == "NA", 0)
            prev4 <- replace(prev4, prev4 == "NaN", 0)
            meas_prevperiods4[i,j] <- sum(prev4)
          } else {
            npos <- sum(test1[1:(j-1)]) + 1
            npos2 <- sum(test1[1:(j)])
            #SSE
            prev1 <- (mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
            prev1 <- replace(prev1, prev1 == "Inf", 0)
            prev1 <- replace(prev1, prev1 == "-Inf", 0)
            prev1 <- replace(prev1, prev1 == "NA", 0)
            prev1 <- replace(prev1, prev1 == "NaN", 0)
            meas_prevperiods1[i,j] <- sum(prev1)
            #MSE
            prev2 <- (mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
            prev2 <- replace(prev2, prev2 == "Inf", 0)
            prev2 <- replace(prev2, prev2 == "-Inf", 0)
            prev2 <- replace(prev2, prev2 == "NA", 0)
            prev2 <- replace(prev2, prev2 == "NaN", 0)
            meas_prevperiods2[i,j] <- sum(prev2)
            #MAE
            prev3 <- abs(mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])
            prev3 <- replace(prev3, prev3 == "Inf", 0)
            prev3 <- replace(prev3, prev3 == "-Inf", 0)
            prev3 <- replace(prev3, prev3 == "NA", 0)
            prev3 <- replace(prev3, prev3 == "NaN", 0)
            meas_prevperiods3[i,j] <- sum(prev3)
            #MAPE
            prev4 <- abs((mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)] - qxt.forecast[[i]][,(npos:npos2)])/mat_qxt[[i]][,(trainset1+npos):(trainset1+npos2)])
            prev4 <- replace(prev4, prev4 == "Inf", 0)
            prev4 <- replace(prev4, prev4 == "-Inf", 0)
            prev4 <- replace(prev4, prev4 == "NA", 0)
            prev4 <- replace(prev4, prev4 == "NaN", 0)
            meas_prevperiods4[i,j] <- sum(prev4)
          }
        }
        meas_periods[1,j] <- sum(meas_prevperiods1[,j])
        meas_periods[2,j] <- sum(meas_prevperiods2[,j])/(nages*nPop*test1[j])
        meas_periods[3,j] <- sum(meas_prevperiods3[,j])/(nages*nPop*test1[j])
        meas_periods[4,j] <- sum(meas_prevperiods4[,j])/(nages*nPop*test1[j])
      }
    }
    else(stop("measures must be equal to SSE, MSE, MAE, MAPE or All."))

  }else if(fixed_train_origin == '1'){

    if(measures == "SSE"){
      #We estimate SSE in different options
      meas_prevpops <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("SSE", c(1:nPop)))

      meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), "SSE"))
      meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("SSE", "all ages and periods"))

      meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "SSE"))
      meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("SSE", ages))

      meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_periods2 <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("SSE", c(1:length(test1))))

      meas_periods <- list('block_SSE' = matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("SSE", c(1:length(test1)))),
                           'SSE_AHEAD' = matrix(NA, nrow = 1, ncol = nahead, dimnames = list("SSE-AHEAD", c(1:nahead))))

      wxt <- genWeightMat(ages = ages, years = c((periods[1] + trainset1):max(df_qxtdata$period)), clip = 0)
      for(i in 1:nPop){
        for(j in 1:(length(test1))){
          meas_prevpops[j,i] <- MeasureAccuracy(measure = "SSE",
                                                qxt_crude = mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))],
                                                qxt_aju = qxt.forecast[[j]][[i]],
                                                wxt = genWeightMat(ages = ages, years = c(1:nahead), clip = 0))$value
          prev <- (mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])^2
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevtotal[j,i] <- sum(prev)
        }
        meas_pops[,i] <- mean(meas_prevpops[,i])
      }
      for(j in 1:(length(test1))){
        meas_prevtot[j,] <- sum(meas_prevtotal[j,])
      }
      meas_total <- mean(meas_prevtot)

      for(j in 1:nages){
        for(pe in 1:length(test1)){
          for(i in 1:nPop){
            prev <- (mat_qxt[[i]][j,(trainset1+pe):(trainset1+nahead+(pe-1))] - qxt.forecast[[pe]][[i]][j,])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevag[i,] <- sum(prev)
          }
          meas_prevages[pe,j] <- sum(meas_prevag)
        }
        meas_ages[,j] <- mean(meas_prevages[,j])
      }

      for(j in 1:length(test1)){
        for(i in 1:nPop){
          prev <- (mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])^2
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevperiods[i,j] <- sum(prev)
        }
        meas_periods$block_SSE[,j] <- sum(meas_prevperiods[,j])
      }

      for(pre in 1:nahead){
        for(j in 1:length(test1)){
          for(i in 1:nPop){
            prev <- (mat_qxt[[i]][,(trainset1+j):(trainset1+j+(pre-1))] - qxt.forecast[[j]][[i]][,1:pre])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
          }
          meas_periods2[,j] <- sum(meas_prevperiods[,j])
        }
        meas_periods$SSE_AHEAD[,pre] <- sum(meas_periods2)}

    }else if(measures == "MSE"){
      #We estimate MSE in different options
      meas_prevpops <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("MSE", c(1:nPop)))

      meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), "MSE"))
      meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("MSE", "all ages and periods"))

      meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "MSE"))
      meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("MSE", ages))

      meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_periods2 <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("MSE", c(1:length(test1))))

      meas_periods <- list('block_MSE' = matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("MSE", c(1:length(test1)))),
                           'MSE_AHEAD' = matrix(NA, nrow = 1, ncol = nahead, dimnames = list("MSE-AHEAD", c(1:nahead))))

      wxt <- genWeightMat(ages = ages, years = c((periods[1] + trainset1):max(df_qxtdata$period)), clip = 0)
      for(i in 1:nPop){
        for(j in 1:(length(test1))){
          meas_prevpops[j,i] <- MeasureAccuracy(measure = "MSE",
                                                qxt_crude = mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))],
                                                qxt_aju = qxt.forecast[[j]][[i]],
                                                wxt = genWeightMat(ages = ages, years = c(1:nahead), clip = 0))$value
          prev <- (mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])^2
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevtotal[,i] <- sum(prev)
        }
        meas_pops[,i] <- mean(meas_prevpops[,i])
      }
      for(j in 1:(length(test1))){
        meas_prevtot[j,] <- sum(meas_prevtotal[j,])/(nages*test1[j]*nPop)
      }
      meas_total <- mean(meas_prevtot)

      for(j in 1:nages){
        for(pe in 1:length(test1)){
          for(i in 1:nPop){
            prev <- (mat_qxt[[i]][j,(trainset1+pe):(trainset1+nahead+(pe-1))] - qxt.forecast[[pe]][[i]][j,])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevag[i,] <- sum(prev)
          }
          meas_prevages[pe,j] <- sum(meas_prevag)/(test1[pe]*nPop)
        }
        meas_ages[,j] <- mean(meas_prevages[,j])
      }

      for(j in 1:length(test1)){
        for(i in 1:nPop){
          prev <- (mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
        }
        meas_periods$block_MSE[,j] <- sum(meas_prevperiods[,j])/(nages*nPop*test1[j])
      }
      for(pre in 1:nahead){
        for(j in 1:length(test1)){
          for(i in 1:nPop){
            prev <- (mat_qxt[[i]][,(trainset1+j):(trainset1+j+(pre-1))] - qxt.forecast[[j]][[i]][,1:pre])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
          }
          meas_periods2[,j] <- sum(meas_prevperiods[,j])
        }
        meas_periods$MSE_AHEAD[,pre] <- sum(meas_periods2)/(nages*nPop*pre)}

    }else if(measures == "MAE"){
      #We estimate MAE in different options
      meas_prevpops <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("MAE", c(1:nPop)))

      meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), "MAE"))
      meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("MAE", "all ages and periods"))

      meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "MAE"))
      meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("MAE", ages))

      meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_periods2 <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("MAE", c(1:length(test1))))

      meas_periods <- list('block_SSE' = matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("MAE", c(1:length(test1)))),
                           'SSE_AHEAD' = matrix(NA, nrow = 1, ncol = nahead, dimnames = list("MAE-AHEAD", c(1:nahead))))

      wxt <- genWeightMat(ages = ages, years = c((periods[1] + trainset1):max(df_qxtdata$period)), clip = 0)

      for(i in 1:nPop){
        for(j in 1:(length(test1))){
          meas_prevpops[j,i] <- MeasureAccuracy(measure = "MAE",
                                                qxt_crude = mat_qxt[[i]][,(trainset1+j):(trainset1 + nahead + (j-1))],
                                                qxt_aju = qxt.forecast[[j]][[i]],
                                                wxt = genWeightMat(ages = ages, years = c(1:nahead), clip = 0))$value
          prev <- abs(mat_qxt[[i]][,(trainset1+j):(trainset1 + nahead + (j-1))] - qxt.forecast[[j]][[i]])
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevtotal[j,i] <- sum(prev)
        }
        meas_pops[,i] <- mean(meas_prevpops[,i])
      }
      for(j in 1:(length(test1))){
        meas_prevtot[j,] <- sum(meas_prevtotal[j,])/(nPop*nages*test1[j])
      }
      meas_total <- mean(meas_prevtot)

      for(j in 1:nages){
        for(pe in 1:length(test1)){
          for(i in 1:nPop){
            prev <- abs(mat_qxt[[i]][j,(trainset1+pe):(trainset1+nahead+(pe-1))] - qxt.forecast[[pe]][[i]][j,])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevag[i,] <- sum(prev)
          }
          meas_prevages[pe,j] <- sum(meas_prevag)/(test1[pe]*nPop)
        }
        meas_ages[,j] <- mean(meas_prevages[,j])
      }

      for(j in 1:length(test1)){
        for(i in 1:nPop){
          prev <- abs(mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
        }
        meas_periods$block_MAE[,j] <- sum(meas_prevperiods[,j])/(nages*nPop*test1[j])
      }
      for(pre in 1:nahead){
        for(j in 1:length(test1)){
          for(i in 1:nPop){
            prev <- abs(mat_qxt[[i]][,(trainset1+j):(trainset1+j+(pre-1))] - qxt.forecast[[j]][[i]][,1:pre])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
          }
          meas_periods2[,j] <- sum(meas_prevperiods[,j])
        }
        meas_periods$SSE_AHEAD[,pre] <- sum(meas_periods2)/(nages*nPop*pre)}

    }else if(measures == "MAPE"){
      #We estimate MAPE in different options
      meas_prevpops <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("MAPE", c(1:nPop)))

      meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), " MAPE"))
      meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("MAPE", "all ages and periods"))

      meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "MAPE"))
      meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("MAPE", ages))

      meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_periods2 <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("MAPE", c(1:length(test1))))

      meas_periods <- list('block_MAPE' = matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("MAPE", c(1:length(test1)))),
                           'MAPE_AHEAD' = matrix(NA, nrow = 1, ncol = nahead, dimnames = list("MAPE-AHEAD", c(1:nahead))))

      wxt <- genWeightMat(ages = ages, years = c((periods[1] + trainset1):max(df_qxtdata$period)), clip = 0)
      for(i in 1:nPop){
        for(j in 1:(length(test1))){
          meas_prevpops[j,i] <- MeasureAccuracy(measure = "MAPE",
                                                qxt_crude = mat_qxt[[i]][,(trainset1+j):(trainset1+nahead(j-1))],
                                                qxt_aju = qxt.forecast[[j]][[i]],
                                                wxt = genWeightMat(ages = ages, years = c(1:nahead), clip = 0))$value
          prev <- abs((mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])/mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))])
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevtotal[j,i] <- sum(prev)
          }
        meas_pops[,i] <- mean(meas_prevpops[,i])
      }
      for(j in 1:(length(test1))){
        meas_prevtot[j,] <- sum(meas_prevtotal[j,])/(nPop*nages*test1[j])
      }
      meas_total <- mean(meas_prevtot)

      for(j in 1:nages){
        for(pe in 1:length(test1)){
          for(i in 1:nPop){
            prev <- (abs((mat_qxt[[i]][j,(trainset1+pe):(trainset1+nahead+(pe-1))] - qxt.forecast[[pe]][[i]][j,])/mat_qxt[[i]][j,(trainset1+pe):(trainset1+nahead+(pe-1))]))
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevag[i,] <- sum(prev)
          }
          meas_prevages[pe,j] <- sum(meas_prevag)/(test1[pe]*nPop)
        }
        meas_ages[,j] <- mean(meas_prevages[,j])
      }

      for(j in 1:length(test1)){
        for(i in 1:nPop){
          prev <- abs((mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])/mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))])
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevperiods[i,j] <- sum(prev)
        }
        meas_periods$block_MAPE[,j] <- sum(meas_prevperiods[,j])/(nages*nPop*test1[j])
      }
      for(pre in 1:nahead){
        for(j in 1:length(test1)){
          for(i in 1:nPop){
            prev <- abs((mat_qxt[[i]][,(trainset1+j):(trainset1+j+(pre-1))] - qxt.forecast[[j]][[i]][,1:pre])/mat_qxt[[i]][,(trainset1+j):(trainset1+j+(pre-1))])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevperiods[i,j] <- sum(prev)
          }
          meas_periods2[,j] <- sum(meas_prevperiods[,j])
        }
        meas_periods$MAPE_AHEAD[,pre] <- sum(meas_periods2)/(nages*nPop*pre)}

    }else if(measures == "All"){
      #We estimate ALL measures in different options
      meas_prevpops1 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevpops2 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevpops3 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevpops4 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_pops <- matrix(NA, nrow = 4, ncol=nPop, dimnames = list(c("SSE", "MSE", "MAE", "MAPE"), c(1:nPop)))

      meas_prevtotal1 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtotal2 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtotal3 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
      meas_prevtotal4 <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))

      meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 4, dimnames = list(c(1:length(test1)),c("SSE", "MSE", "MAE", "MAPE")))
      meas_total <- matrix(NA, nrow = 4, ncol = 1, dimnames = list(c("SSE", "MSE", "MAE", "MAPE"), "all ages and periods"))

      meas_prevages1 <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevages2 <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevages3 <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevages4 <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
      meas_prevag <- matrix(NA, nrow = nPop, ncol = 4, dimnames = list(c(1:nPop), c("SSE", "MSE", "MAE", "MAPE")))
      meas_ages <- matrix(NA, nrow = 4, ncol = nages, dimnames = list(c("SSE", "MSE", "MAE", "MAPE"), ages))

      meas_prevperiods1 <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_prevperiods2 <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_prevperiods3 <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
      meas_prevperiods4 <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))

      meas_periods2 <- matrix(NA, nrow = 4, ncol = length(test1), dimnames = list(c("SSE", "MSE", "MAE", "MAPE"), c(1:length(test1))))

      meas_periods <- list('block' = matrix(NA, nrow = 4, ncol = length(test1), dimnames = list(c("SSE", "MSE", "MAE", "MAPE"), c(1:length(test1)))),
                           'ahead' = matrix(NA, nrow = 4, ncol = nahead, dimnames = list(c("SSE", "MSE", "MAE", "MAPE"), c(1:nahead))))

      wxt <- genWeightMat(ages = ages, years = c((periods[1] + trainset1):max(df_qxtdata$period)), clip = 0)
      for(i in 1:nPop){
        for(j in 1:(length(test1))){
            meas_prevpops1[j,i] <- MeasureAccuracy(measure = "SSE", qxt_crude = mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))],
                                                   qxt_aju = qxt.forecast[[j]][[i]],
                                                   wxt = genWeightMat(ages = ages, years = c(1:nahead), clip = 0))$value
            meas_prevpops2[j,i] <- MeasureAccuracy(measure = "MSE", qxt_crude = mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))],
                                                   qxt_aju = qxt.forecast[[j]][[i]],
                                                   wxt = genWeightMat(ages = ages, years = c(1:nahead), clip = 0))$value
            meas_prevpops3[j,i] <- MeasureAccuracy(measure = "MAE", qxt_crude = mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))],
                                                   qxt_aju = qxt.forecast[[j]][[i]],
                                                   wxt = genWeightMat(ages = ages, years = c(1:nahead), clip = 0))$value
            meas_prevpops4[j,i] <- MeasureAccuracy(measure = "MAPE", qxt_crude = mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))],
                                                   qxt_aju = qxt.forecast[[j]][[i]],
                                                   wxt = genWeightMat(ages = ages, years = c(1:nahead), clip = 0))$value
            #SSE
            prev1 <- (mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])^2
            prev1 <- replace(prev1, prev1 == "Inf", 0)
            prev1 <- replace(prev1, prev1 == "-Inf", 0)
            prev1 <- replace(prev1, prev1 == "NA", 0)
            prev1 <- replace(prev1, prev1 == "NaN", 0)
            meas_prevtotal1[j,i] <- sum(prev1)
            #MSE
            prev2 <- (mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])^2
            prev2 <- replace(prev2, prev2 == "Inf", 0)
            prev2 <- replace(prev2, prev2 == "-Inf", 0)
            prev2 <- replace(prev2, prev2 == "NA", 0)
            prev2 <- replace(prev2, prev2 == "NaN", 0)
            meas_prevtotal2[j,i] <- sum(prev2)
            #MAE
            prev3 <- abs(mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])
            prev3 <- replace(prev3, prev3 == "Inf", 0)
            prev3 <- replace(prev3, prev3 == "-Inf", 0)
            prev3 <- replace(prev3, prev3 == "NA", 0)
            prev3 <- replace(prev3, prev3 == "NaN", 0)
            meas_prevtotal3[j,i] <- sum(prev3)
            #MAPE
            prev4 <- abs((mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])/mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))])
            prev4 <- replace(prev4, prev4 == "Inf", 0)
            prev4 <- replace(prev4, prev4 == "-Inf", 0)
            prev4 <- replace(prev4, prev4 == "NA", 0)
            prev4 <- replace(prev4, prev4 == "NaN", 0)
            meas_prevtotal4[j,i] <- sum(prev4)
        }
        meas_pops[1,i] <- mean(meas_prevpops1[,i])
        meas_pops[2,i] <- mean(meas_prevpops2[,i])
        meas_pops[3,i] <- mean(meas_prevpops3[,i])
        meas_pops[4,i] <- mean(meas_prevpops4[,i])
      }
      for(j in 1:(length(test1))){
        meas_prevtot[j,1] <- sum(meas_prevtotal1[j,])
        meas_prevtot[j,2] <- sum(meas_prevtotal2[j,])/(nPop*nages*test1[j])
        meas_prevtot[j,3] <- sum(meas_prevtotal3[j,])/(nPop*nages*test1[j])
        meas_prevtot[j,4] <- sum(meas_prevtotal4[j,])/(nPop*nages*test1[j])
      }
      meas_total[1,1] <- mean(meas_prevtot[,1])
      meas_total[2,1] <- mean(meas_prevtot[,2])
      meas_total[3,1] <- mean(meas_prevtot[,3])
      meas_total[4,1] <- mean(meas_prevtot[,4])

      for(j in 1:nages){
        for(pe in 1:length(test1)){
          for(i in 1:nPop){
              #SSE
              prev <- (mat_qxt[[i]][j,(trainset1+pe):(trainset1+nahead+(pe-1))] - qxt.forecast[[pe]][[i]][j,])^2
              prev1 <- replace(prev1, prev1 == "Inf", 0)
              prev1 <- replace(prev1, prev1 == "-Inf", 0)
              prev1 <- replace(prev1, prev1 == "NA", 0)
              prev1 <- replace(prev1, prev1 == "NaN", 0)
              meas_prevag[i,1] <- sum(prev1)

              #MSE
              prev2 <- (mat_qxt[[i]][j,(trainset1+pe):(trainset1+nahead+(pe-1))] - qxt.forecast[[pe]][[i]][j,])^2
              prev2 <- replace(prev2, prev2 == "Inf", 0)
              prev2 <- replace(prev2, prev2 == "-Inf", 0)
              prev2 <- replace(prev2, prev2 == "NA", 0)
              prev2 <- replace(prev2, prev2 == "NaN", 0)
              meas_prevag[i,2] <- sum(prev2)

              #MAE
              prev3 <- abs(mat_qxt[[i]][j,(trainset1+pe):(trainset1+nahead+(pe-1))] - qxt.forecast[[pe]][[i]][j,])
              prev3 <- replace(prev3, prev3 == "Inf", 0)
              prev3 <- replace(prev3, prev3 == "-Inf", 0)
              prev3 <- replace(prev3, prev3 == "NA", 0)
              prev3 <- replace(prev3, prev3 == "NaN", 0)
              meas_prevag[i,3] <- sum(prev3)
              #MAPE
              prev4 <- abs((mat_qxt[[i]][j,(trainset1+pe):(trainset1+nahead+(pe-1))] - qxt.forecast[[pe]][[i]][j,])/mat_qxt[[i]][j,(trainset1+pe):(trainset1+nahead+(pe-1))])
              prev4 <- replace(prev4, prev4 == "Inf", 0)
              prev4 <- replace(prev4, prev4 == "-Inf", 0)
              prev4 <- replace(prev4, prev4 == "NA", 0)
              prev4 <- replace(prev4, prev4 == "NaN", 0)
              meas_prevag[i,4] <- sum(prev4)
          }
          meas_prevages1[pe,j] <- sum(meas_prevag[,1])
          meas_prevages2[pe,j] <- sum(meas_prevag[,2])/(test1[pe]*nPop)
          meas_prevages3[pe,j] <- sum(meas_prevag[,3])/(test1[pe]*nPop)
          meas_prevages4[pe,j] <- sum(meas_prevag[,4])/(test1[pe]*nPop)
        }
        meas_ages[1,j] <- mean(meas_prevages1[,j])
        meas_ages[2,j] <- mean(meas_prevages2[,j])
        meas_ages[3,j] <- mean(meas_prevages3[,j])
        meas_ages[4,j] <- mean(meas_prevages4[,j])
      }

      for(j in 1:length(test1)){
        for(i in 1:nPop){
            #SSE
            prev1 <- (mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])^2
            prev1 <- replace(prev1, prev1 == "Inf", 0)
            prev1 <- replace(prev1, prev1 == "-Inf", 0)
            prev1 <- replace(prev1, prev1 == "NA", 0)
            prev1 <- replace(prev1, prev1 == "NaN", 0)
            meas_prevperiods1[i,j] <- sum(prev1)
            #MSE
            prev2 <- (mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])^2
            prev2 <- replace(prev2, prev2 == "Inf", 0)
            prev2 <- replace(prev2, prev2 == "-Inf", 0)
            prev2 <- replace(prev2, prev2 == "NA", 0)
            prev2 <- replace(prev2, prev2 == "NaN", 0)
            meas_prevperiods2[i,j] <- sum(prev2)
            #MAE
            prev3 <- abs(mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])
            prev3 <- replace(prev3, prev3 == "Inf", 0)
            prev3 <- replace(prev3, prev3 == "-Inf", 0)
            prev3 <- replace(prev3, prev3 == "NA", 0)
            prev3 <- replace(prev3, prev3 == "NaN", 0)
            meas_prevperiods3[i,j] <- sum(prev3)
            #MAPE
            prev4 <- abs((mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))] - qxt.forecast[[j]][[i]])/mat_qxt[[i]][,(trainset1+j):(trainset1+nahead+(j-1))])
            prev4 <- replace(prev4, prev4 == "Inf", 0)
            prev4 <- replace(prev4, prev4 == "-Inf", 0)
            prev4 <- replace(prev4, prev4 == "NA", 0)
            prev4 <- replace(prev4, prev4 == "NaN", 0)
            meas_prevperiods4[i,j] <- sum(prev4)
        }
        meas_periods$block[1,j] <- sum(meas_prevperiods1[,j])
        meas_periods$block[2,j] <- sum(meas_prevperiods2[,j])/(nages*nPop*test1[j])
        meas_periods$block[3,j] <- sum(meas_prevperiods3[,j])/(nages*nPop*test1[j])
        meas_periods$block[4,j] <- sum(meas_prevperiods4[,j])/(nages*nPop*test1[j])
      }
      for(pre in 1:nahead){
      for(j in 1:length(test1)){
        for(i in 1:nPop){

            #SSE
            prev1 <- (mat_qxt[[i]][,(trainset1+j):(trainset1+j+(pre-1))] - qxt.forecast[[j]][[i]][,1:pre])^2
            prev1 <- replace(prev1, prev1 == "Inf", 0)
            prev1 <- replace(prev1, prev1 == "-Inf", 0)
            prev1 <- replace(prev1, prev1 == "NA", 0)
            prev1 <- replace(prev1, prev1 == "NaN", 0)
            meas_prevperiods1[i,j] <- sum(prev1)
            #MSE
            prev2 <- (mat_qxt[[i]][,(trainset1+j):(trainset1+j+(pre-1))] - qxt.forecast[[j]][[i]][,1:pre])^2
            prev2 <- replace(prev2, prev2 == "Inf", 0)
            prev2 <- replace(prev2, prev2 == "-Inf", 0)
            prev2 <- replace(prev2, prev2 == "NA", 0)
            prev2 <- replace(prev2, prev2 == "NaN", 0)
            meas_prevperiods2[i,j] <- sum(prev2)
            #MAE
            prev3 <- abs(mat_qxt[[i]][,(trainset1+j):(trainset1+j+(pre-1))] - qxt.forecast[[j]][[i]][,1:pre])
            prev3 <- replace(prev3, prev3 == "Inf", 0)
            prev3 <- replace(prev3, prev3 == "-Inf", 0)
            prev3 <- replace(prev3, prev3 == "NA", 0)
            prev3 <- replace(prev3, prev3 == "NaN", 0)
            meas_prevperiods3[i,j] <- sum(prev3)
            #MAPE
            prev4 <- abs((mat_qxt[[i]][,(trainset1+j):(trainset1+j+(pre-1))] - qxt.forecast[[j]][[i]][,1:pre])/mat_qxt[[i]][,(trainset1+j):(trainset1+j+(pre-1))])
            prev4 <- replace(prev4, prev4 == "Inf", 0)
            prev4 <- replace(prev4, prev4 == "-Inf", 0)
            prev4 <- replace(prev4, prev4 == "NA", 0)
            prev4 <- replace(prev4, prev4 == "NaN", 0)
            meas_prevperiods4[i,j] <- sum(prev4)
        }
        meas_periods2[1,j] <- sum(meas_prevperiods1[,j])
        meas_periods2[2,j] <- sum(meas_prevperiods2[,j])
        meas_periods2[3,j] <- sum(meas_prevperiods3[,j])
        meas_periods2[4,j] <- sum(meas_prevperiods4[,j])
      }
        meas_periods$ahead[1,pre] <- sum(meas_periods2[1,])
        meas_periods$ahead[2,pre] <- sum(meas_periods2[2,])/(nages*nPop*pre)
        meas_periods$ahead[3,pre] <- sum(meas_periods2[3,])/(nages*nPop*pre)
        meas_periods$ahead[4,pre] <- sum(meas_periods2[4,])/(nages*nPop*pre)
      }
    }else(stop("measures must be equal to SSE, MSE, MAE, MAPE or All."))

  }

  if(model != "ACFM"){
    warn <- NULL
  }

  #Now, we start estimating the measures of forecasting-accuracy
  return <- list(ax = ax,
                 bx = bx,
                 kt.fitted = kt,
                 kt.future = kt.fut,
                 kt.arima = kt.arima,
                 Ii = Ii,
                 formula = object$formula,
                 model = object$model,
                 nPop = nPop,
                 Ages = ages,
                 Periods = periods,
                 qxt.crude = mat_qxt,
                 qxt.forecast = qxt.forecast,
                 logit.qxt.forecast = logit.qxt.forecast,
                 meas_ages = meas_ages,
                 meas_periodsfut = meas_periods,
                 meas_pop = meas_pops,
                 meas_total = meas_total,
                 fixed_train_origin = fixed_train_origin,
                 CV_method = CV,
                 warn_msgs = warn)
  class(return) <- "MultiCv"
  return

}
#' @export
print.MultiCv <- function(x, ...) {
  if(!is.null(x)){
    if(!"MultiCv" %in% class(x))
      stop("The 'x' does not have the 'MultiCv' structure of R CvmortalityMult package.")
  }
  if(!is.list(x)){
    stop("The 'x' is not a list. Use 'MultiCv' function first.")
  }

  if(x$nPop != 1){
    if(x$model == "additive"){
      cat("Fitting the additive multi-population mortality model: \n")
    } else if(x$model == "multiplicative"){
      cat("Fitting the multiplicative multi-population mortality model: \n")
    } else if(x$model == "CFM"){
      cat("Fitting the common-factor multi-population mortality model: \n")
    } else if(x$model == "joint-K"){
      cat("Fitting the joint-K multi-population mortality model: \n")
    } else if(x$model == "ACFM"){
      cat("Fitting the augmented-common-factor multi-population mortality model: \n")
    }
  } else if(x$nPop == 1){
    cat("Fitting the single-population version of the Lee-Carter model: \n")
  }
  print(x$formula)
  cat(paste0("\nWe employ the"), x$CV_method, "(with fixed train origin =", x$fixed_train_origin, "); and using as:\n")
  cat(paste("\nFitting and Forecasting periods:", min(x$Periods), "-", max(x$Periods), "\n"))
  cat(paste("\nFitting and Forecasting ages:", min(x$Ages), "-", max(x$Ages), "\n"))
  cat(paste("\nFitting populations:", x$nPop, "\n"))
  cat(paste("\nPeriods using as train and test sets:"))
  for(i in 1:length(names(x$kt.future))){
    cat(paste("\n",i, "train set", as.numeric(rownames(x$kt.fitted[[i]]))[1], "-",
              as.numeric(rownames(x$kt.fitted[[i]]))[length(rownames(x$kt.fitted[[i]]))],
              "and test set", as.numeric(rownames(x$kt.future[[i]])[1]), "-",
              as.numeric(rownames(x$kt.future[[i]])[length(rownames(x$kt.future[[i]]))])))
  }
  if(length(rownames(x$meas_ages)) != 1){
    cat(paste("\nMeasure of forecasting accuracy selected: SSE, MSE, MAE and MAPE\n"))

  }else{
    cat(paste("\nMeasure of forecasting accuracy selected:", rownames(x$meas_ages), "\n"))}

  cat(paste("\nGlobal Forecasting accuracy over ages, periods and populations", round(x$meas_total, 6), "using", rownames(x$meas_ages),"\n"))
}
