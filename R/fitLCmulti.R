#' Function to fit multi-population mortality models
#' @description
#' R function for fitting additive or multiplicative multi-population mortality model developed by: Debon et al. (2011) and Russolillo et al. (2011), respectively.
#' These model follows the structure of the well-known Lee-Carter model (Lee and Carter, 1992) but including an additive or multiplicative parameter to capture the behavior of each population considered.
#' This parameter seeks to capture the individual behavior of every population considered.
#' It should be mentioned that this function is developed for fitting several populations.
#' However, in case you only consider one population, the function will fit the single population version of the Lee-Carter model, the classical one.
#'
#' @param model multi-population mortality model chosen to fit the mortality rates c("`additive`", "`multiplicative`"). In case you do not provide any value, the function will apply the "`additive`" option.
#' @param qxt mortality rates used to fit the additive multipopulation mortality model. This rates can be provided in matrix or in data.frame.
#' @param periods periods considered in the fitting in a vector way c(minyear:maxyear).
#' @param ages vector with the ages considered in the fitting. If the mortality rates provide from an abridged life tables, it is necessary to provide a vector with the ages, see the example.
#' @param nPop number of population considered for fitting.
#' @param lxt survivor function considered for every population, not necessary to provide.
#'
#' @return A list with class \code{"LCmulti"} including different components of the fitting process:
#' * `ax` parameter that captures the average shape of the mortality curve in all considered populations.
#' * `bx` parameter that explains the age effect x with respect to the general trend `kt` in the mortality rates of all considered populations.
#' * `kt` represent the national tendency of multi-mortality populations during the period.
#' * `Ii` gives an idea of the differences in the pattern of mortality in any region i with respect to Region 1.
#' * `formula` additive multi-population mortality formula used to fit the mortality rates.
#' * `model` provided the model selected in every case.
#' * `data.used` mortality rates used to fit the data.
#' * `qxt.real` real mortality rates.
#' * `qxt.fitted` fitted mortality rates using the additive multi-population mortality model.
#' * `logit.qxt.fitted` fitted mortality rates in logit way.
#' * `Ages` provided ages to fit the data.
#' * `Periods` provided periods to fit the periods.
#' * `nPop` provided number of populations to fit the periods.
#'
#' @seealso \code{\link{forecast.fitLCmulti}},
#' \code{\link{multipopulation_cv}}, \code{\link{multipopulation_loocv}},
#' \code{\link{plot.fitLCmulti}}
#'
#' @references
#' Debon, A., Montes, F., and Martinez-Ruiz, F. (2011).
#' Statistical methods to compare mortality for a group with non-divergent populations: an application to Spanish regions.
#' European Actuarial Journal, 1, 291-308.
#'
#' Lee, R.D. and Carter, L.R. (1992).
#' Modeling and forecasting US mortality.
#' Journal of the American Statistical Association, 87(419), 659–671.
#'
#' Russolillo, M., Giordano, G., & Haberman, S. (2011).
#' Extending the Lee–Carter model: a three-way decomposition.
#' Scandinavian Actuarial Journal, 2011(2), 96-117.
#'
#' @import gnm
#' @importFrom utils install.packages
#' @importFrom stats coef
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
#' #MULTIPLICATIVE MULTI-POPULATION MORTALITY MODEL
#' #In the case, the user wants to fit the multiplicative multi-population mortality model
#' multiplicative_Spainmales <- fitLCmulti(model = "multiplicative",
#'                               qxt = SpainRegions$qx_male,
#'                               periods = c(1991:2020),
#'                               ages = c(ages),
#'                               nPop = 18,
#'                               lxt = SpainRegions$lx_male)
#'
#'
#' multiplicative_Spainmales
#'
#' #Once, we have fit the data, it is possible to see the ax, bx, kt, and It
#' #provided parameters for the fitting.
#' plot(multiplicative_Spainmales)
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
#' }
#' @export
fitLCmulti <- function(model = c("additive", "multiplicative"),
                       qxt, periods, ages, nPop, lxt = NULL){

  #Check several things before start
  if(is.null(model) || is.null(qxt) || is.null(periods) || is.null(ages) || is.null(nPop)){
    warning("Arguments model, qxt, periods, ages, and nPop, need to be provided.")
  }

  #1. Check that the model corresponds to the additive o multiplicative multi-population mortality model
  valid_model <- c("additive", "multiplicative")
  model <- match.arg(model, valid_model)
  #In case, you donot provide any value of model, the function applies the additive one.

  #2. Check that periods and ages are two vectors
  if (!is.vector(periods) || !is.vector(ages)) {
    warning("Period or Ages are not a vector, need to be provided.")
  }

  nperiods <- length(periods)
  nages <- length(ages)

  #3. Check if lxt is provided
  if(is.null(lxt)){
    message("Argument lxt will be obtained as the number of individuals at age
    x and period t, starting with l0=100,000")
  }

  #We will use plogis() for inverse logit
  #plogis

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
      df_actual <- data.frame(rep((i-1), nages*nperiods),
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
    qxt1 <- df_qxtdata[df_qxtdata$pop == (i-1),]$qxt

    mat_qxt[[paste0("pob", i)]] <- matrix(qxt1, nrow = nages, ncol = nperiods,
                                          dimnames = list(ages, periods))
  }

  #Fitting multi-population or single mortality model
  #1. Certify the number of populations:
  #1.1 Number of population different than one (higher than one).
  #The function will apply a multi-population mortality model that can be additive or multiplicative
  if(nPop != 1){
    #Additive multi-population mortality model
    if(model == "additive"){

      emptymodel <- gnm(qxt ~ -1 + as.factor(age), weights=df_qxtdata$lx,
                        family= 'quasibinomial', data=df_qxtdata)

      biplotStart <- residSVD2(emptymodel,
                              factor(df_qxtdata$age),
                              factor(df_qxtdata$period),1)

      #We need to ensure that strating values are equal to b[1]=1 and k[1]=0
      aini <- coef(emptymodel)+biplotStart[1]*biplotStart[(nages+1)]*biplotStart[1:nages]/biplotStart[1]
      bini <- biplotStart[1:nages]/biplotStart[1]
      kini <- biplotStart[1]*biplotStart[(nages+1):(nages+nperiods)]-biplotStart[1]*biplotStart[(nages+1)]

      ###################
      #Adittive model
      lee.geo3 <- gnm(qxt ~ -1 + factor(age) + Mult(factor(period),factor(age)) + factor(pop),
                    weights = df_qxtdata$lx, family = 'quasibinomial',
                    constrain=c((nages+1), (nages+nperiods+1)), constrainTo=c(0,1),
                    data=df_qxtdata, start=c(aini,kini,bini,rep(0,(nPop-1))))

      ax.geo3 <- as.numeric(coef(lee.geo3)[1:nages])
      kt.geo3 <- as.numeric(c(0,coef(lee.geo3)[(nages+2):(nages+nperiods)]))
      bx.geo3 <- as.numeric(c(1,coef(lee.geo3)[(nages+nperiods+2):(nages+nperiods+nages)]))
      Ii.geo3 <- as.numeric(c(0,coef(lee.geo3)[(nages+nperiods+nages+1):length(coef(lee.geo3))]))

      #Multiplicative multi-population mortality model
    } else if(model == "multiplicative") {
      emptymodel <- gnm(qxt ~ -1 + factor(age), weights=df_qxtdata$lx,
                        family= 'quasibinomial', data=df_qxtdata)
      biplotStart <- residSVD2(emptymodel,
                              factor(df_qxtdata$age),
                              factor(df_qxtdata$period),1)

      #We need to ensure that starting values are equal to b[1]=1 and k[1]=0
      aini <- as.numeric(coef(emptymodel)+biplotStart[1]*biplotStart[(nages+1)]*biplotStart[1:nages]/biplotStart[1])
      bini <- as.numeric(biplotStart[1:nages]/biplotStart[1])
      kini <- as.numeric(biplotStart[1]*biplotStart[(nages+1):(nages+nperiods)]-biplotStart[1]*biplotStart[(nages+1)])
      Iini <- as.numeric(rep(0, nPop))

      ###################
      #Multiplicative model
      lee.geo3 <- gnm(qxt ~ -1 + factor(age) + Mult(factor(period), factor(age), factor(pop)),
                      weights = df_qxtdata$lx, family = 'quasibinomial',
                      constrain=c((nages+1), (nages+nperiods+1), (nages+nperiods+nages+1)),
                      constrainTo=c(0, 1, 1),
                      data=df_qxtdata, start=c(aini, bini, kini, exp(Iini)))
      ax.geo3 <- as.numeric(coef(lee.geo3)[1:nages])
      kt.geo3 <- as.numeric(c(0,coef(lee.geo3)[(nages+2):(nages+nperiods)]))
      bx.geo3 <- as.numeric(c(1,coef(lee.geo3)[(nages+nperiods+2):(nages+nperiods+nages)]))
      Ii.geo3 <- as.numeric(c(1,coef(lee.geo3)[(nages+nperiods+nages+2):length(coef(lee.geo3))]))
    }
    #For one-single population the function applies the single-population version of the LC, the classical one.
  } else {
    message("You only provide one country. Thus, we fit LC one-single model population.")
    emptymodel <- gnm(qxt ~ -1 + factor(age), weights = df_qxtdata$lx,
                      family='quasibinomial', data=df_qxtdata)
    biplotStart <- residSVD2(emptymodel,
                            factor(df_qxtdata$age), factor(df_qxtdata$period),1)
    aini <- coef(emptymodel)+biplotStart[1]*biplotStart[(nages+1)]*biplotStart[1:nages]/biplotStart[1]
    bini <- biplotStart[1:nages]/biplotStart[1]
    kini <- biplotStart[1]*biplotStart[(nages+1):(nages+nperiods)]-biplotStart[1]*biplotStart[(nages+1)]

    lee.geo3 <- gnm(qxt ~ -1 + factor(age) + Mult(factor(period), factor(age)),
                    weights = df_qxtdata$lx, family = 'quasibinomial',
                    constrain = c((nages+1), (nages+nperiods+1)), constrainTo=c(0,1),
                    data = df_qxtdata, start = c(aini, bini, kini))
    ax.geo3 <- as.numeric(coef(lee.geo3)[1:nages])
    bx.geo3 <- as.numeric(c(1,coef(lee.geo3)[(nages+nperiods+2):(nages+nperiods+nages)]))
    kt.geo3 <- as.numeric(c(0,coef(lee.geo3)[(nages+2):(nages+nperiods)]))


  }

  if(nPop == 1){
    Ii.geo3 <- matrix(1, nrow = nages, ncol = nperiods)
    formula.used <- paste("Binomial model with predictor: logit q[x,t,i] = a[x] + b[x] k[t]")
    Ii23 <- matrix(1, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "Ipop"))
    model = "LC-single-pop"
  } else{
    if(model == "additive"){
      formula.used <- paste("Binomial model with predictor: logit q[x,t,i] = a[x] + b[x] k[t] + I[i]")
      Ii23 <- matrix(Ii.geo3, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "Ipop"))
    } else if(model == "multiplicative") {
      formula.used <- paste("Binomial model with predictor: logit q[x,t,i] = a[x] + b[x] k[t] I[i]")
      Ii23 <- matrix(Ii.geo3, nrow = nPop, ncol = 1, dimnames = list(c(1:nPop), "Ipop"))
    }
  }

  #Create fitted values of qxt for every population provided
  lee.logit <- list()
  lee.qxt <- list()
  #it <- 1

  if(model == "additive"){
    for(it in 1:nPop){
        lee.logit[[paste0("pob", it)]] <- matrix(rep(ax.geo3, nperiods), nrow= nages, ncol = nperiods)+
          (matrix(bx.geo3, nrow=nages, ncol=1)%*%matrix(kt.geo3,nrow=1, ncol=nperiods)) +
          matrix(rep(Ii.geo3[it], nages*nperiods), nrow= nages, ncol = nperiods)

        lee.qxt[[paste0("pob", it)]] <- plogis(lee.logit[[it]])
        rownames(lee.logit[[it]]) <- rownames(lee.qxt[[it]]) <- ages
        colnames(lee.logit[[it]]) <- colnames(lee.qxt[[it]]) <- periods
  }} else if(model == "multiplicative") {
    for(it in 1:nPop){
      lee.logit[[paste0("pob", it)]] <- matrix(rep(ax.geo3, nperiods), nrow= nages, ncol = nperiods)+
        (matrix(bx.geo3, nrow=nages, ncol=1)%*%matrix(kt.geo3,nrow=1, ncol=nperiods))*
        matrix(rep(Ii.geo3[it], nages*nperiods), nrow= nages, ncol = nperiods)

      lee.qxt[[paste0("pob", it)]] <- plogis(lee.logit[[it]])
      rownames(lee.logit[[it]]) <- rownames(lee.qxt[[it]]) <- ages
      colnames(lee.logit[[it]]) <- colnames(lee.qxt[[it]]) <- periods

    }} else if(model == "LC-single-pop"){
      for(it in 1:nPop){
        lee.logit[[paste0("pob", it)]] <- matrix(rep(ax.geo3, nperiods), nrow= nages, ncol = nperiods)+
          (matrix(bx.geo3, nrow=nages, ncol=1)%*%matrix(kt.geo3,nrow=1, ncol=nperiods))

        lee.qxt[[paste0("pob", it)]] <- plogis(lee.logit[[it]])
        rownames(lee.logit[[it]]) <- rownames(lee.qxt[[it]]) <- ages
        colnames(lee.logit[[it]]) <- colnames(lee.qxt[[it]]) <- periods
  }}


  #Check the length of qxt is acording to the provided data
  return <- list(ax = matrix(ax.geo3, nrow = 1, ncol = nages, dimnames = list("ax", ages)),
                 bx = matrix(bx.geo3, nrow = 1, ncol = nages, dimnames = list("bx", ages)),
                 kt = matrix(kt.geo3, nrow = nperiods, ncol = 1, dimnames = list(c(min(periods):max(periods)), "kt")),
                 Ii = Ii23,
                 formula = formula.used,
                 model = model,
                 data.used = df_qxtdata,
                 qxt.real = mat_qxt,
                 qxt.fitted = lee.qxt,
                 logit.qxt.fitted = lee.logit,
                 Ages = ages,
                 Periods = periods,
                 nPop = nPop)
  class(return) <- "fitLCmulti"
  return

}
#' @export
print.fitLCmulti <- function(x, ...) {
  if(!is.null(x)){
    if(!"fitLCmulti" %in% class(x))
      stop("The 'x' does not have the 'fitLCmulti' structure of R CvmortalityMult package.")
  }
  if(!is.list(x)){
    stop("The 'x' is not a list. Use 'fitLCmulti' function first.")
  }

  if(x$nPop != 1){
    if(x$model == "additive"){
      cat("Fitting the additive multi-population mortality model: \n")
    } else if(x$model == "multiplicative"){
      cat("Fitting the multiplicative multi-population mortality model: \n")
    }
  } else if(x$nPop == 1){
    cat("Fitting the single-population version of the Lee-Carter model: \n")
  }
  print(x$formula)
  cat(paste("\nFitting years:", min(x$Periods), "-", max(x$Periods)))
  cat(paste("\nFitting ages:", min(x$Ages), "-", max(x$Ages), "\n"))
  cat(paste("\nFitting populations:", x$nPop, "\n"))
}
