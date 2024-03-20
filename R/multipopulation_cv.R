#' FUNCTION TO APPLY CROSS-VALIDATION TECHNIQUES FOR TESTING THE FORECASTING ACCURACY
#' OF MULTI-POPULATION MORTALITY MODELS
#' @description
#' R function for testing the accuracy out-of-sample of different multi-population mortality models, Additive (Debon et al., 2011) and Multiplicative (Russolillo et al., 2011).
#' We provide a R function that employ the cross-validation techniques for data panel-time series (Atance et al. 2020) to test the forecasting accuracy.
#' These techniques consist on split the database in two parts: training set (to run the model) and test set (to check the forecasting accuracy of the model).
#' This procedure is repeated several times trying to check the forecasting accuracy in different ways.
#' With this function, the user can provide its own mortality rates for different populations. The function will split the database chronologically (Bergmeir and Benitez, 2012) based on the nahead which consist on the length of the training set.
#' We have include the following Figure 1 to understand how the R function works.
#' {\figure{CV_technique.jpg}{options: width="100\%" alt="Figure: mai.png"}}
#' It should be mentioned that this function is developed for cross-validation the forecasting accuracy of several populations.
#' However, in case you only consider one population, the function will forecast the Lee-Carter model for one population.
#' To test the forecasting accuracy of the selected model, the function provides five different measures: SSE, MSE, MAE, MAPE or All. Depending on how you want to check the forecasting accuracy of the model you could select one or other.
#' In this case, the measures will be obtained using the mortality rates in the normal scale as recommended by Santolino (2023) against the log scale.
#'
#' @param qxt mortality rates used to fit the multi-population mortality models. This rates can be provided in matrix or in data.frame.
#' @param model choose the multi-population mortality model to fit the mortality rates c("`additive`", "`multiplicative`")
#' @param periods periods considered in the fitting in a vector way c(`minyear`:`maxyear`).
#' @param ages vector with the ages considered in the fitting. If the mortality rates provide from an abridged life tables, it is necessary to provide a vector with the ages, see the example.
#' @param nPop number of population considered for fitting.
#' @param lxt survivor function considered for every population, not necessary to provide.
#' @param nahead is a vector specifying the number of periods to block in the blocked CV. The function operates by using the sum of the periods in nahead and three (the minimum number of years required to construct a time series), as the initial training set. This ensures that the first train set has sufficient observations to forecast the initial test set, which will be of length `nahead`.
#' @param ktmethod method used to forecast the value of `kt` Arima(p,d,q) or ARIMA(0,1,0); c("`Arimapdq`", "`arima010`").
#' @param kt_include.cte if you want that `kt` include constant in the arima process.
#' @param measures choose the non-penalized measure of forecasting accuracy that you want to use; c("`SSE`", "`MSE`", "`MAE`", "`MAPE`", "`All`"). Check the function
#'
#' @return A list with different components of the cross-validation process:
#' * `ax` parameter that captures the average shape of the mortality curve in all considered populations.
#' * `bx` parameter that explains the age effect x with respect to the general trend `kt` in the mortality rates of all considered populations.
#' * `kt.fitted` obtained values for the tendency behavior captured by `kt` .
#' * `kt.future` future values of `kt` for every iteration in the cross-validation.
#' * `kt.arima`  the arima selected for each `kt` time series.
#' * `Ii` parameter that captures the differences in the pattern of mortality in any region i with respect to Region 1.
#' * `formula` multi-population mortality formula used to fit the mortality rates.
#' * `nPop` provided number of populations to fit the periods.
#' * `qxt.real` real mortality rates.
#' * `qxt.future` future mortality rates estimated with the multi-population mortality model.
#' * `logit.qxt.future` future mortality rates in logit way estimated with the multi-population mortality model.
#' * `meas_ages` measure of forecasting accuracy through the ages of the study.
#' * `meas_periodsfut` measure of forecasting accuracy in every forecasting period(s) of the study.
#' * `meas_pop` measure of forecasting accuracy through the populations considered in the study.
#' * `meas_total` a global measure of forecasting accuracy through the ages, periods and populations of the study.
#'
#' @seealso \code{\link{multipopulation_loocv}}, \code{\link{fit_additive.LC.multi}}, \code{\link{fit_multiplicative.LC.multi}},
#' \code{\link{for_additive.LC.multi}}, \code{\link{for_multiplicative.LC.multi}},
#' \code{\link{plotLC.multi}}, \code{\link{SSE}}, \code{\link{MAE}}, \code{\link{MAPE}}.
#'
#' @references
#' Atance, D., Debon, A., and Navarro, E. (2020).
#' A comparison of forecasting mortality models using resampling methods.
#' Mathematics 8(9): 1550.
#'
#' Bergmeir, C. & Benitez, J.M. (2012)
#' On the use of cross-validation for time series predictor evaluation.
#' Information Sciences, 191, 192–213.
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
#'
#' @examples
#'
#' #The example takes more than 5 seconds because it includes
#' #several fitting and forecasting process and hence all
#' #the process is included in donttest
#' \donttest{
#' #We present a cross-validation method for spanish male regions
# 'SpainRegions
#'
#' ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40,
#'          45, 50, 55, 60, 65, 70, 75, 80, 85, 90)
#' library(gnm)
#' #Let start with a simple nahead=5 CV method obtaining the SSE forecasting measure of accuracy
#' cv_Spainmales_addit <- multipopulation_cv(qxt = SpainRegions$qx_male,
#'                                          model = c("additive"),
#'                                          periods =  c(1991:2020), ages = c(ages),
#'                                          nPop = 18, lxt = SpainRegions$lx_male,
#'                                          nahead = 5,
#'                                          ktmethod = c("Arimapdq"),
#'                                          kt_include.cte = TRUE,
#'                                          measures = c("SSE"))
#'
#' #Once, we have run the function we can check the result in different ways:
#' cv_Spainmales_addit$meas_ages
#' cv_Spainmales_addit$meas_periodsfut
#' cv_Spainmales_addit$meas_pop
#' cv_Spainmales_addit$meas_total
#' }
#' @export
multipopulation_cv <- function(qxt, model = c("additive", "multiplicative"),
                               periods, ages, nPop, lxt=NULL,
                               nahead,
                               ktmethod = c("Arimapdq", "arima010"),
                               kt_include.cte = TRUE,
                               measures = c("SSE", "MSE", "MAE", "MAPE", "All")){
  #Check several things before start
  if(is.null(qxt) || is.null(periods) || is.null(ages) ||
     is.null(nPop) || is.null(model) || is.null(nahead)){
    stop(warning("Arguments qxt, periods, ages, nPop, and model, need to be provided."))
  }

  #2. Check that periods and ages are two vectors
  if (!is.vector(periods) || !is.vector(ages)) {
    stop(warning("Period or Ages are not a vector, need to be provided."))
  }

  if(!is.numeric(nahead)){
    stop(warning("nahead must be numeric variable."))
  }
  #Construct the inv.logit -- function
  inv.logit <- function(x){exp(x)/(1+exp(x))}
  logit<-function(x){log(x/(1-x))}

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
        stop(warning("population ", i, " has a different dim regarding the rest of the chosen populations."))
      }
    }
    #Check the size of the qxt is equal to number of ages, periods and populations provided
    nper.age.pop <- nperiods*nages*nPop
    if(length(qxt[[i]])*length(qxt) != nper.age.pop){
      stop(warning("Number of qxt is different from the period, ages and countries provided."))
    }

    #Check the size of list(qxt) is the same the npop
    if(length(qxt) != nPop){
      stop(warning("Number of n Pop is different regarding the length of the matrix qxt."))
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
      stop(warning("Number of qxt vector is different from the period, ages and countries provided."))
    }
    message("Your qxt and lxt data are in vector form.\n")
    message("So, please ensure that qxt are provided by age, period and population, as follows:\n")
    message("age0-period0-pob1, age1-period0-pob1, ..., age.n-period0-pob1\n")
    message("age0-period1-pob1, age1-period1-pob1, ..., age.n-period1-pob1\n")
    message("...\n")
    message("age0-period0-pob2, age1-period0-pob2, ..., age.n-period0-pob2\n")

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

  if(nahead > length(periods)){
    stop(warning("number of periods ahead is higher than the total number of periods."))
  }

  nper <- 2*nahead + 3

  if(length(periods) < nper){
    stop(warning("number of periods has to be bigger than 2*nahead+3."))
  }

  #Also, I am going to transform qxt into a matrix to provide as result
  mat_qxt <- list()
  for(i in 1:nPop){
    nperiods
    nages
    qxt1 <- df_qxtdata[df_qxtdata$pop == (i-1),]$qxt

    mat_qxt[[paste0("pob", i)]] <- matrix(qxt1, nrow = nages, ncol = nperiods,
                                          dimnames = list(ages, periods))


  }

  #Estimate the first set of fitting-sample
  fitting.years1 <- 3 + nahead
  nrep <- (nperiods - fitting.years1)/nahead
  prue1 <- all.equal(nrep, as.integer(nrep))

  #we estimate the number of times (loops) that the model has to be fitted
  #and hence, the number of times to forecast the model
  if(prue1 == TRUE){
    nrep <- as.integer(nrep)
    test1 <- c(rep(nahead, nrep))
  } else{
    nrep <- floor(nrep) + 1
    test1 <- c(rep(nahead, nrep - 1))
    test2 <- length(periods) - fitting.years1 - (nrep-1)*nahead
    test1 <- c(test1, test2)
  }

  #we create a list of matrix to save qxt for each forecast population
  qxt.forecast <- list()
  for(i in 1:nPop){
    qxt.forecast[[paste0("pob", i)]] <- matrix(NA, nrow = nages, ncol = sum(test1),
                                               dimnames = list(ages, c((periods[1] + fitting.years1):max(df_qxtdata$period))))

  }
  logit.qxt.forecast <- qxt.forecast

  fitting.data <- forecasting.data <- list()
  ax <- bx <- kt <- Ii <- kt.fut <- kt.arima <- list()
  #reps <- 1
  for(reps in 1:nrep){
    fitting.years1 <- 3 + nahead*reps

    periods2 <- periods[1] + fitting.years1
    df_qxtdata2 <- df_qxtdata[df_qxtdata$period < periods2,]

    if(model == "additive"){
      fitted.obj <- fit_additive.LC.multi(qxt = df_qxtdata2$qxt,
                                          periods = c(min(periods):(periods2-1)),
                                          ages = ages,
                                          nPop = nPop, lxt = df_qxtdata2$lxt)
      forecast.obj <- for_additive.LC.multi(fitted.obj,
                                                 nahead = test1[reps],
                                                 ktmethod = ktmethod,
                                                 kt_include.cte = kt_include.cte)

    } else if(model == "multiplicative"){
      fitted.obj <- fit_multiplicative.LC.multi(qxt = df_qxtdata2$qxt,
                                                periods = c(min(periods):(periods2-1)),
                                                ages = ages,
                                                nPop = nPop, lxt = df_qxtdata2$lxt)
      forecast.obj <- for_multiplicative.LC.multi(fitted.obj,
                                                 nahead = test1[reps],
                                                 ktmethod = ktmethod,
                                                 kt_include.cte = kt_include.cte)


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

    fitting.data[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- fitted.obj
    forecasting.data[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- forecast.obj

    ax[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- fitted.obj$ax
    bx[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- fitted.obj$bx
    kt[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- fitted.obj$kt
    Ii[[paste0("loop-", reps, " from ", periods[1], " to ", periods2-1)]] <- fitted.obj$Ii

    kt.fut[[paste0("loop-", reps, " from ", periods2, " to ", (periods2+test1[reps]-1))]] <- forecast.obj$kt.fut
    kt.arima[[paste0("loop-", reps, " from ", periods2, " to ", (periods2+test1[reps]-1))]] <- forecast.obj$arimakt

  }

  #We estimate the logit for the forecast qxt
  for(i in 1:nPop){
    logit.qxt.forecast[[i]] <- logit(qxt.forecast[[i]])
  }

  if(measures == "SSE"){
    #We estimate SSE in different options
    meas_prevpops <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
    meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("sse", c(1:nPop)))

    meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
    meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), "mse"))
    meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("sse", "all ages and periods"))

    meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
    meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "sse"))
    meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("sse", ages))

    meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
    meas_periods <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("sse", c(1:length(test1))))

    wxt <- genWeightMat(ages = ages, years = c((periods[1] + 3 + nahead):max(df_qxtdata$period)), clip = 0)
    for(i in 1:nPop){
      for(j in 1:(length(test1))){
        if(j == 1){
          npos <- 1
          meas_prevpops[j,i] <- SSE(qxt_re = mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])],
                                    qxt_aju = qxt.forecast[[i]][,(1:nahead[1])],
                                    wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))
          prev <- (mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])^2
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevtotal[j,i] <- sum(prev)

        } else {
          npos <- sum(test1[1:(j-1)]) + 1
          npos2 <- sum(test1[1:(j)])
          meas_prevpops[j,i] <- SSE(qxt_re = mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)],
                                    qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                    wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))
          prev <- (mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
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
            prev <- (mat_qxt[[i]][j,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][j,(1:nahead[1])])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevag[i,] <- sum(prev)

          } else {
            npos <- sum(test1[1:(pe-1)]) + 1
            npos2 <- sum(test1[1:(pe)])
            prev <- (mat_qxt[[i]][j,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])^2
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
          prev <- (mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])^2
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevperiods[i,j] <- sum(prev)

        } else {
          npos <- sum(test1[1:(j-1)]) + 1
          npos2 <- sum(test1[1:(j)])
          prev <- (mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
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
    meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("mse", c(1:nPop)))

    meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
    meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), "mse"))
    meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("mse", "all ages and periods"))

    meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
    meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "mse"))
    meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("mse", ages))

    meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
    meas_periods <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("mse", c(1:length(test1))))

    wxt <- genWeightMat(ages = ages, years = c((periods[1] + 3 + nahead):max(df_qxtdata$period)), clip = 0)
    for(i in 1:nPop){
      for(j in 1:(length(test1))){
        if(j == 1){
          npos <- 1
          meas_prevpops[j,i] <- MSE(qxt_re = mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])],
                                    qxt_aju = qxt.forecast[[i]][,(1:nahead[1])],
                                    wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))
          prev <- (mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])^2
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevtotal[,i] <- sum(prev)
        } else {
          npos <- sum(test1[1:(j-1)]) + 1
          npos2 <- sum(test1[1:(j)])
          meas_prevpops[j,i] <- MSE(qxt_re = mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)],
                                    qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                    wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))
          prev <- (mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
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
            prev <- (mat_qxt[[i]][j,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][j,(1:nahead[1])])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevag[i,] <- sum(prev)

          } else {
            npos <- sum(test1[1:(pe-1)]) + 1
            npos2 <- sum(test1[1:(pe)])
            prev <- (mat_qxt[[i]][j,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])^2
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevag[i,] <- sum(prev)
          }
        }
        meas_prevages[pe,j] <- sum(meas_prevag)/(test1[j]*nPop)
      }
      meas_ages[,j] <- mean(meas_prevages[,j])
    }

    for(j in 1:length(test1)){
      for(i in 1:nPop){
        if(j == 1){
          npos <- 1
          prev <- (mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])^2
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevperiods[i,j] <- sum(prev)
        } else {
          npos <- sum(test1[1:(j-1)]) + 1
          npos2 <- sum(test1[1:(j)])
          prev <- (mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
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
    meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("mae", c(1:nPop)))

    meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
    meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), "mae"))
    meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("mae", "all ages and periods"))

    meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
    meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "mae"))
    meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("mae", ages))

    meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
    meas_periods <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("mae", c(1:length(test1))))

    wxt <- genWeightMat(ages = ages, years = c((periods[1] + 3 + nahead):max(df_qxtdata$period)), clip = 0)

    for(i in 1:nPop){
      for(j in 1:(length(test1))){
        if(j == 1){
          npos <- 1
          meas_prevpops[j,i] <- MAE(qxt_re = mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])],
                                    qxt_aju = qxt.forecast[[i]][,(1:nahead[1])],
                                    wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))
          prev <- abs(mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevtotal[j,i] <- sum(prev)
        } else {
          npos <- sum(test1[1:(j-1)]) + 1
          npos2 <- sum(test1[1:(j)])
          meas_prevpops[j,i] <- MAE(qxt_re = mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)],
                                    qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                    wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))
          prev <- abs(mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])
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
            prev <- abs(mat_qxt[[i]][j,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][j,(1:nahead[1])])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevag[i,] <- sum(prev)
          } else {
            npos <- sum(test1[1:(pe-1)]) + 1
            npos2 <- sum(test1[1:(pe)])
            prev <- abs(mat_qxt[[i]][j,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevag[i,] <- sum(prev)
          }
        }
        meas_prevages[pe,j] <- sum(meas_prevag)/(test1[j]*nPop)
      }
      meas_ages[,j] <- mean(meas_prevages[,j])
    }

    for(j in 1:length(test1)){
      for(i in 1:nPop){
        if(j == 1){
          npos <- 1
          prev <- abs(mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevperiods[i,j] <- sum(prev)
        } else {
          npos <- sum(test1[1:(j-1)]) + 1
          npos2 <- sum(test1[1:(j)])
          prev <- abs(mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])
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
    meas_pops <- matrix(NA, nrow = 1, ncol=nPop, dimnames = list("sse", c(1:nPop)))

    meas_prevtotal <- matrix(NA, nrow = length(test1), ncol=nPop, dimnames = list(c(1:length(test1)), c(1:nPop)))
    meas_prevtot <- matrix(NA, nrow = length(test1), ncol = 1, dimnames = list(c(1:length(test1)), "mse"))
    meas_total <- matrix(NA, nrow = 1, ncol = 1, dimnames = list("sse", "all ages and periods"))

    meas_prevages <- matrix(NA, nrow = length(test1), ncol = nages, dimnames = list(c(1:length(test1)), ages))
    meas_prevag <- matrix(NA, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "sse"))
    meas_ages <- matrix(NA, nrow = 1, ncol = nages, dimnames = list("sse", ages))

    meas_prevperiods <- matrix(NA, nrow = nPop, ncol = length(test1), dimnames = list(c(1:nPop), c(1:length(test1))))
    meas_periods <- matrix(NA, nrow = 1, ncol = length(test1), dimnames = list("sse", c(1:length(test1))))

    wxt <- genWeightMat(ages = ages, years = c((periods[1] + 3 + nahead):max(df_qxtdata$period)), clip = 0)
    for(i in 1:nPop){
      for(j in 1:(length(test1))){
        if(j == 1){
          npos <- 1
          meas_prevpops[j,i] <- MAPE(qxt_re = mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])],
                                     qxt_aju = qxt.forecast[[i]][,(1:nahead[1])],
                                     wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))
          prev <- abs((mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])/mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])])
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevtotal[j,i] <- sum(prev)
        } else {
          npos <- sum(test1[1:(j-1)]) + 1
          npos2 <- sum(test1[1:(j)])
          meas_prevpops[j,i] <- MAPE(qxt_re = mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)],
                                     qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                     wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))
          prev <- abs((mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])/mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)])
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
            prev <- (abs((mat_qxt[[i]][j,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][j,(1:nahead[1])])/mat_qxt[[i]][j,(3+nahead+1):(3+nahead+test1[1])]))
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevag[i,] <- sum(prev)

          } else {
            npos <- sum(test1[1:(pe-1)]) + 1
            npos2 <- sum(test1[1:(pe)])
            prev <- (abs((mat_qxt[[i]][j,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])/mat_qxt[[i]][j,(3+nahead+npos):(3+nahead+npos2)]))
            prev <- replace(prev, prev == "Inf", 0)
            prev <- replace(prev, prev == "-Inf", 0)
            prev <- replace(prev, prev == "NA", 0)
            prev <- replace(prev, prev == "NaN", 0)
            meas_prevag[i,] <- sum(prev)
          }
        }
        meas_prevages[pe,j] <- sum(meas_prevag)/(test1[j]*nPop)
      }
      meas_ages[,j] <- mean(meas_prevages[,j])
    }

    for(j in 1:length(test1)){
      for(i in 1:nPop){
        if(j == 1){
          npos <- 1
          prev <- abs((mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])/mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])])
          prev <- replace(prev, prev == "Inf", 0)
          prev <- replace(prev, prev == "-Inf", 0)
          prev <- replace(prev, prev == "NA", 0)
          prev <- replace(prev, prev == "NaN", 0)
          meas_prevperiods[i,j] <- sum(prev)

        } else {
          npos <- sum(test1[1:(j-1)]) + 1
          npos2 <- sum(test1[1:(j)])
          prev <- abs((mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])/mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)])
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
    wxt <- genWeightMat(ages = ages, years = c((periods[1] + 3 + nahead):max(df_qxtdata$period)), clip = 0)
    for(i in 1:nPop){
      for(j in 1:(length(test1))){
        if(j == 1){
          npos <- 1
          meas_prevpops1[j,i] <- SSE(qxt_re = mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])],
                                     qxt_aju = qxt.forecast[[i]][,(1:nahead[1])],
                                     wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))
          meas_prevpops2[j,i] <- MSE(qxt_re = mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])],
                                     qxt_aju = qxt.forecast[[i]][,(1:nahead[1])],
                                     wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))
          meas_prevpops3[j,i] <- MAE(qxt_re = mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])],
                                     qxt_aju = qxt.forecast[[i]][,(1:nahead[1])],
                                     wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))
          meas_prevpops4[j,i] <- MAPE(qxt_re = mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])],
                                      qxt_aju = qxt.forecast[[i]][,(1:nahead[1])],
                                      wxt = genWeightMat(ages = ages, years = c(1:test1[1]), clip = 0))
         #SSE
          prev1 <- (mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])^2
          prev1 <- replace(prev1, prev1 == "Inf", 0)
          prev1 <- replace(prev1, prev1 == "-Inf", 0)
          prev1 <- replace(prev1, prev1 == "NA", 0)
          prev1 <- replace(prev1, prev1 == "NaN", 0)
          meas_prevtotal1[j,i] <- sum(prev1)
          #MSE
          prev2 <- (mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])^2
          prev2 <- replace(prev2, prev2 == "Inf", 0)
          prev2 <- replace(prev2, prev2 == "-Inf", 0)
          prev2 <- replace(prev2, prev2 == "NA", 0)
          prev2 <- replace(prev2, prev2 == "NaN", 0)
          meas_prevtotal2[j,i] <- sum(prev2)
          #MAE
          prev3 <- abs(mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])
          prev3 <- replace(prev3, prev3 == "Inf", 0)
          prev3 <- replace(prev3, prev3 == "-Inf", 0)
          prev3 <- replace(prev3, prev3 == "NA", 0)
          prev3 <- replace(prev3, prev3 == "NaN", 0)
          meas_prevtotal3[j,i] <- sum(prev3)
          #MAPE
          prev4 <- abs((mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])/mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])])
          prev4 <- replace(prev4, prev4 == "Inf", 0)
          prev4 <- replace(prev4, prev4 == "-Inf", 0)
          prev4 <- replace(prev4, prev4 == "NA", 0)
          prev4 <- replace(prev4, prev4 == "NaN", 0)
          meas_prevtotal4[j,i] <- sum(prev4)
        } else {
          npos <- sum(test1[1:(j-1)]) + 1
          npos2 <- sum(test1[1:(j)])
          meas_prevpops1[j,i] <- SSE(qxt_re = mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)],
                                     qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                     wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))
          meas_prevpops2[j,i] <- MSE(qxt_re = mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)],
                                     qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                     wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))
          meas_prevpops3[j,i] <- MAE(qxt_re = mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)],
                                     qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                     wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))
          meas_prevpops4[j,i] <- MAPE(qxt_re = mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)],
                                      qxt_aju = qxt.forecast[[i]][,(npos:npos2)],
                                      wxt = genWeightMat(ages = ages, years = c(1:test1[j]), clip = 0))
          #SSE
          prev1 <- (mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
          prev1 <- replace(prev1, prev1 == "Inf", 0)
          prev1 <- replace(prev1, prev1 == "-Inf", 0)
          prev1 <- replace(prev1, prev1 == "NA", 0)
          prev1 <- replace(prev1, prev1 == "NaN", 0)
          meas_prevtotal1[j,i] <- sum(prev1)
          #MSE
          prev2 <- (mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
          prev2 <- replace(prev2, prev2 == "Inf", 0)
          prev2 <- replace(prev2, prev2 == "-Inf", 0)
          prev2 <- replace(prev2, prev2 == "NA", 0)
          prev2 <- replace(prev2, prev2 == "NaN", 0)
          meas_prevtotal2[j,i] <- sum(prev2)
          #MAE
          prev3 <- abs(mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])
          prev3 <- replace(prev3, prev3 == "Inf", 0)
          prev3 <- replace(prev3, prev3 == "-Inf", 0)
          prev3 <- replace(prev3, prev3 == "NA", 0)
          prev3 <- replace(prev3, prev3 == "NaN", 0)
          meas_prevtotal3[j,i] <- sum(prev3)
          #MAPE
          prev4 <- abs((mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])/mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)])
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
            prev <- (mat_qxt[[i]][j,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][j,(1:nahead[1])])^2
            prev1 <- replace(prev1, prev1 == "Inf", 0)
            prev1 <- replace(prev1, prev1 == "-Inf", 0)
            prev1 <- replace(prev1, prev1 == "NA", 0)
            prev1 <- replace(prev1, prev1 == "NaN", 0)
            meas_prevag[i,1] <- sum(prev1)

            #MSE
            prev2 <- (mat_qxt[[i]][j,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][j,(1:nahead[1])])^2
            prev2 <- replace(prev2, prev2 == "Inf", 0)
            prev2 <- replace(prev2, prev2 == "-Inf", 0)
            prev2 <- replace(prev2, prev2 == "NA", 0)
            prev2 <- replace(prev2, prev2 == "NaN", 0)
            meas_prevag[i,2] <- sum(prev2)

            #MAE
            prev3 <- abs(mat_qxt[[i]][j,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][j,(1:nahead[1])])
            prev3 <- replace(prev3, prev3 == "Inf", 0)
            prev3 <- replace(prev3, prev3 == "-Inf", 0)
            prev3 <- replace(prev3, prev3 == "NA", 0)
            prev3 <- replace(prev3, prev3 == "NaN", 0)
            meas_prevag[i,3] <- sum(prev3)
            #MAPE
            prev4 <- abs((mat_qxt[[i]][j,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][j,(1:nahead[1])])/mat_qxt[[i]][j,(3+nahead+1):(3+nahead+test1[1])])
            prev4 <- replace(prev4, prev4 == "Inf", 0)
            prev4 <- replace(prev4, prev4 == "-Inf", 0)
            prev4 <- replace(prev4, prev4 == "NA", 0)
            prev4 <- replace(prev4, prev4 == "NaN", 0)
            meas_prevag[i,4] <- sum(prev4)
          } else {
            npos <- sum(test1[1:(pe-1)]) + 1
            npos2 <- sum(test1[1:(pe)])
            prev <- (mat_qxt[[i]][j,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])^2
            #SSE
            prev <- (mat_qxt[[i]][j,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])^2
            prev1 <- replace(prev1, prev1 == "Inf", 0)
            prev1 <- replace(prev1, prev1 == "-Inf", 0)
            prev1 <- replace(prev1, prev1 == "NA", 0)
            prev1 <- replace(prev1, prev1 == "NaN", 0)
            meas_prevag[i,1] <- sum(prev1)

            #MSE
            prev2 <- (mat_qxt[[i]][j,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])^2
            prev2 <- replace(prev2, prev2 == "Inf", 0)
            prev2 <- replace(prev2, prev2 == "-Inf", 0)
            prev2 <- replace(prev2, prev2 == "NA", 0)
            prev2 <- replace(prev2, prev2 == "NaN", 0)
            meas_prevag[i,2] <- sum(prev2)

            #MAE
            prev3 <- abs(mat_qxt[[i]][j,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])
            prev3 <- replace(prev3, prev3 == "Inf", 0)
            prev3 <- replace(prev3, prev3 == "-Inf", 0)
            prev3 <- replace(prev3, prev3 == "NA", 0)
            prev3 <- replace(prev3, prev3 == "NaN", 0)
            meas_prevag[i,3] <- sum(prev3)
            #MAPE
            prev4 <- abs((mat_qxt[[i]][j,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][j,(npos:npos2)])/mat_qxt[[i]][j,(3+nahead+npos):(3+nahead+npos2)])
            prev4 <- replace(prev4, prev4 == "Inf", 0)
            prev4 <- replace(prev4, prev4 == "-Inf", 0)
            prev4 <- replace(prev4, prev4 == "NA", 0)
            prev4 <- replace(prev4, prev4 == "NaN", 0)
            meas_prevag[i,4] <- sum(prev4)
          }
        }
        meas_prevages1[pe,j] <- sum(meas_prevag[,1])
        meas_prevages2[pe,j] <- sum(meas_prevag[,2])/(test1[j]*nPop)
        meas_prevages3[pe,j] <- sum(meas_prevag[,3])/(test1[j]*nPop)
        meas_prevages4[pe,j] <- sum(meas_prevag[,4])/(test1[j]*nPop)
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
          prev1 <- (mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])^2
          prev1 <- replace(prev1, prev1 == "Inf", 0)
          prev1 <- replace(prev1, prev1 == "-Inf", 0)
          prev1 <- replace(prev1, prev1 == "NA", 0)
          prev1 <- replace(prev1, prev1 == "NaN", 0)
          meas_prevperiods1[i,j] <- sum(prev1)
          #MSE
          prev2 <- (mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])^2
          prev2 <- replace(prev2, prev2 == "Inf", 0)
          prev2 <- replace(prev2, prev2 == "-Inf", 0)
          prev2 <- replace(prev2, prev2 == "NA", 0)
          prev2 <- replace(prev2, prev2 == "NaN", 0)
          meas_prevperiods2[i,j] <- sum(prev2)
          #MAE
          prev3 <- abs(mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])
          prev3 <- replace(prev3, prev3 == "Inf", 0)
          prev3 <- replace(prev3, prev3 == "-Inf", 0)
          prev3 <- replace(prev3, prev3 == "NA", 0)
          prev3 <- replace(prev3, prev3 == "NaN", 0)
          meas_prevperiods3[i,j] <- sum(prev3)
          #MAPE
          prev4 <- abs((mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])] - qxt.forecast[[i]][,(1:nahead[1])])/mat_qxt[[i]][,(3+nahead+1):(3+nahead+test1[1])])
          prev4 <- replace(prev4, prev4 == "Inf", 0)
          prev4 <- replace(prev4, prev4 == "-Inf", 0)
          prev4 <- replace(prev4, prev4 == "NA", 0)
          prev4 <- replace(prev4, prev4 == "NaN", 0)
          meas_prevperiods4[i,j] <- sum(prev4)
        } else {
          npos <- sum(test1[1:(j-1)]) + 1
          npos2 <- sum(test1[1:(j)])
          #SSE
          prev1 <- (mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
          prev1 <- replace(prev1, prev1 == "Inf", 0)
          prev1 <- replace(prev1, prev1 == "-Inf", 0)
          prev1 <- replace(prev1, prev1 == "NA", 0)
          prev1 <- replace(prev1, prev1 == "NaN", 0)
          meas_prevperiods1[i,j] <- sum(prev1)
          #MSE
          prev2 <- (mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])^2
          prev2 <- replace(prev2, prev2 == "Inf", 0)
          prev2 <- replace(prev2, prev2 == "-Inf", 0)
          prev2 <- replace(prev2, prev2 == "NA", 0)
          prev2 <- replace(prev2, prev2 == "NaN", 0)
          meas_prevperiods2[i,j] <- sum(prev2)
          #MAE
          prev3 <- abs(mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])
          prev3 <- replace(prev3, prev3 == "Inf", 0)
          prev3 <- replace(prev3, prev3 == "-Inf", 0)
          prev3 <- replace(prev3, prev3 == "NA", 0)
          prev3 <- replace(prev3, prev3 == "NaN", 0)
          meas_prevperiods3[i,j] <- sum(prev3)
          #MAPE
          prev4 <- abs((mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)] - qxt.forecast[[i]][,(npos:npos2)])/mat_qxt[[i]][,(3+nahead+npos):(3+nahead+npos2)])
          prev4 <- replace(prev4, prev4 == "Inf", 0)
          prev4 <- replace(prev4, prev4 == "-Inf", 0)
          prev4 <- replace(prev4, prev4 == "NA", 0)
          prev4 <- replace(prev4, prev4 == "NaN", 0)
          meas_prevperiods4[i,j] <- sum(prev4)
        }
      }
      meas_periods[1,j] <- sum(meas_prevperiods1[,j])
      meas_periods[1,j] <- sum(meas_prevperiods2[,j])/(nages*nPop*test1[j])
      meas_periods[1,j] <- sum(meas_prevperiods3[,j])/(nages*nPop*test1[j])
      meas_periods[1,j] <- sum(meas_prevperiods4[,j])/(nages*nPop*test1[j])
    }
  }else(stop("measures must be equal to SSE, MSE, MAE, MAPE or All"))

  return <- list(ax = ax,
                 bx = bx,
                 kt.fitted = kt,
                 kt.future = kt.fut,
                 kt.arima = kt.arima,
                 Ii = Ii,
                 formula = fitted.obj$formula,
                 nPop = nPop,
                 qxt.real = mat_qxt,
                 qxt.forecast = qxt.forecast,
                 logit.qxt.forecast = logit.qxt.forecast,
                 meas_ages = meas_ages,
                 meas_periodsfut = meas_periods,
                 meas_pop = meas_pops,
                 meas_total = meas_total)

}
