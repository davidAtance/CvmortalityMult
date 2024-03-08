library(devtools)
library(roxygen2)
library(Rdpack)
library(usethis)
library(testthat)
library(available)
available::available("CvmortalityMult", browse = FALSE)

## easilyr

#Construct help files
#roxygenise()
devtools::load_all()

usethis::use_testthat()
devtools::test_coverage()

?SSE
?plot.fit.LC.multi
?regions
?SpainMap
?multipopulation_loocv

#rm("SpainMap")
#ℹ Loading CvmortalityMult

#The data that we are going to use:
SpainRegions
?SpainRegions
#Mortality Data
#Spain Regions for males and females
#Years 1991 : 2020
#Abbridged Ages 0 : 90

ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40,
          45, 50, 55, 60, 65, 70, 75, 80, 85, 90)

####################################################################
#AJUSTAMOS MODELO ADITIVO
####################################################################
?fit.additive.LC.multi
additive_Spainmales <- fit.additive.LC.multi(qxt = SpainRegions$qx_male,
                      periods = c(1991:2020),
                      ages = c(ages),
                      nPop = 18,
                      lxt = SpainRegions$lx_male)
additive_Spainmales$formula

?plot.fit.LC.multi
plot.fit.LC.multi(additive_Spainmales$Ii)

?SpainMap
Ii_title = expression(I[i])

SpainMap(regionvalue = additive_Spainmales$Ii[2:18],
         main = c("Additive for males"),
         name = c("Ii"))

additive_Spainfemales <- fit.additive.LC.multi(qxt = SpainRegions$qx_female,
                                               periods = c(1991:2020),
                                               ages = c(ages),
                                               nPop = 18)
#' Once, we have fit the data, it is possible to see the ax, bx, kt, and Ii provided parameters for the fitting.
plot.fit.LC.multi(additive_Spainfemales)
SpainMap(regionvalue = additive_Spainfemales$Ii[2:18],
         main = c("Additive for females"),
         name = c(expression(I[i])))

#Now, we fit only the national population of Spain (without the regions) - only one-single population
SpainNat
?SpainNat

LC_Spainmales <- fit.additive.LC.multi(qxt = SpainNat$qx_male,
                                       periods = c(1991:2020),
                                       ages = ages,
                                       nPop = 1)
plot.fit.LC.multi(LC_Spainmales)

####################################################################
#AJUSTAMOS MODELO MULTIPLICATIVO
####################################################################
?fit.multiplicative.LC.multi

multiplicative_Spainmales <- fit.multiplicative.LC.multi(qxt = SpainRegions$qx_male,
                               periods = c(1991:2020),
                               ages = c(ages),
                               nPop = 18,
                               lxt = SpainRegions$lx_male)
SpainMap(regionvalue = multiplicative_Spainmales$Ii[2:18],
         main = c("Multiplicative for males"),
         name = c("Ii"))

# Once, we have fit the data, it is possible to see the \eqn{\alpha_x}, \eqn{\beta_x}, \eqn{k_t}, and \eqn{I_i} provided parameters for the fitting.
plot.fit.LC.multi(multiplicative_Spainmales)
x#'
#' Equal to the previous step but in this case for females and without providing lxt.
multiplicative_Spainfemales <- fit.multiplicative.LC.multi(qxt = SpainRegions$qx_female,
                                 periods = c(1991:2020),
                                 ages = c(ages),
                                 nPop = 18)
#' Once, we have fit the data, it is possible to see the ax, bx, kt, and Ii provided parameters for the fitting.
plot.fit.LC.multi(multiplicative_Spainfemales)
SpainMap(regionvalue = multiplicative_Spainfemales$Ii[2:18],
         main = c("Multiplicative for females"),
         name = c("Ii"))

#'
#' As we mentioned in the details of the function, if we only provide the data from one-population the function
#' \code{\link{fit.multiplicative.LC.multi}} will fit the Lee-Carter model for single populations.
LC_Spainmales <- fit.multiplicative.LC.multi(qxt = SpainNat$qx_male,
                               periods = c(1991:2020),
                               ages = ages,
                               nPop = 1)
plot.fit.LC.multi(LC_Spainmales)

####################################################################
#FORECAST ADDITIVE LC MULTI para una poblacion
####################################################################
SpainRegions

ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40,
          45, 50, 55, 60, 65, 70, 75, 80, 85, 90)

###################################################################
#AJUSTAMOS MODELO ADITIVO
###################################################################
additive_Spainmales <- fit.additive.LC.multi(qxt = SpainRegions$qx_male,
                                             periods = c(1991:2020),
                                             ages = c(ages),
                                             nPop = 18,
                                             lxt = SpainRegions$lx_male)

additive_Spainmales$formula

plot.fit.additive.LC.multi(additive_Spainmales)

?forecast.additive.LC.multi

fut_additive_Spainmales <- forecast.additive.LC.multi(fitted.obj = additive_Spainmales, nahead = 10,
                           ktmethod = "Arimapdq", kt_include.cte = TRUE)

#Igual pero para 1 unica poblacion que ajusta el LC basico
LC_Spainmales <- fit.additive.LC.multi(qxt = SpainNat$qx_male,
                               periods = c(1991:2020),
                               ages = ages,
                               nPop = 1)
plot.fit.additive.LC.multi(LC_Spainmales)

LC_Spainmales$formula

fut_LC_Spainmales <- forecast.additive.LC.multi(fitted.obj = LC_Spainmales, nahead = 10,
                                                ktmethod = "Arimapdq", kt_include.cte = TRUE)
fut_LC_Spainmales$formula

###################################################################
#AJUSTAMOS MODELO MULTIPLICATIVO
###################################################################
multiplicative_Spainmales <- fit.multiplicative.LC.multi(qxt = SpainRegions$qx_male,
                                             periods = c(1991:2020),
                                             ages = c(ages),
                                             nPop = 18,
                                             lxt = SpainRegions$lx_male)

plot.fit.additive.LC.multi(multiplicative_Spainmales)

?fit.multiplicative.LC.multi

fut_multiplicative_Spainmales <- forecast.multiplicative.LC.multi(fitted.obj = multiplicative_Spainmales, nahead = 10,
                                                      ktmethod = "Arimapdq", kt_include.cte = TRUE)

fut_multiplicative_Spainmales$formula

#Igual pero para 1 unica poblacion que ajusta el LC basico
LC_Spainmales <- fit.multiplicative.LC.multi(qxt = SpainNat$qx_male,
                                       periods = c(1991:2020),
                                       ages = ages,
                                       nPop = 1)
plot.fit.additive.LC.multi(LC_Spainmales)

LC_Spainmales$formula

fut_LC_Spainmales <- forecast.multiplicative.LC.multi(fitted.obj = LC_Spainmales, nahead = 10,
                                                ktmethod = "Arimapdq", kt_include.cte = TRUE)
fut_LC_Spainmales$formula

###################################################################
#MULTIPOPULATION_CROSS-VALIDATION
###################################################################
#Vamos hacer la prueba con las tecnicas de validacion cruzada
SpainRegions

ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40,
          45, 50, 55, 60, 65, 70, 75, 80, 85, 90)

?multipopulation_cv
#Tambien se puede ver las diferentes medidas de error que vamos aplicar, con las fórmulas que se aplican y demás
?SSE
?MSE
?MAE
?MAPE

 #Let start with a simple nahead=5 CV method obtaining the SSE forecasting measure of accuracy
cv_Spainmales_addit <- multipopulation_cv(qxt = SpainRegions$qx_male,
                                    model = c("additive"),
                                    periods =  c(1991:2020), ages = c(ages),
                                    nPop = 18, lxt = SpainRegions$lx_male,
                                    nahead = 5,
                                    ktmethod = c("Arimapdq"),
                                    kt_include.cte = TRUE,
                                    measures = c("SSE"))
cv_Spainmales_addit$meas_ages
cv_Spainmales_addit$meas_periodsfut
cv_Spainmales_addit$meas_pop
cv_Spainmales_addit$meas_total

#Equally but with multiplicative
cv_Spainmales_multi <- multipopulation_cv(qxt = SpainRegions$qx_male,
                                          model = c("multiplicative"),
                                          periods =  c(1991:2020), ages = c(ages),
                                          nPop = 18, lxt = SpainRegions$lx_male,
                                          nahead = 5,
                                          ktmethod = c("Arimapdq"),
                                          kt_include.cte = TRUE,
                                          measures = c("SSE"))


cv_Spainmales_multi$meas_ages
cv_Spainmales_multi$meas_periodsfut
cv_Spainmales_multi$meas_pop
cv_Spainmales_multi$meas_total

#Ahora repetimos procedimiento pero para una poblacion:
cv_Spainmales_LC <- multipopulation_cv(qxt = SpainNat$qx_male,
                                        model = c("multiplicative"),
                                        periods =  c(1991:2020), ages = c(ages),
                                       nPop = 1, lxt = SpainNat$lx_male,
                                       nahead = 5,
                                       ktmethod = c("Arimapdq"),
                                       kt_include.cte = TRUE,
                                       measures = c("SSE"))
cv_Spainmales_LC$meas_ages
cv_Spainmales_LC$meas_periodsfut
cv_Spainmales_LC$meas_pop
cv_Spainmales_LC$meas_total

################################################################
#Voy aplicar el LOOCV del articulo del Risk, asi que voy ajustar un data.frame para cuadrarlo
SpReg <- data.frame(ccaa = SpainRegions$ccaa,
                       years = SpainRegions$years,
                       ages = SpainRegions$ages,
                       qx_male = SpainRegions$qx_male,
                       qx_female = SpainRegions$qx_female,
                       lx_male = SpainRegions$lx_male,
                       lx_female = SpainRegions$lx_female,
                       series = SpainRegions$series,
                       label = SpainRegions$label)

SpainReg90.20 <- SpReg[SpReg$years < 2001,]
qxt = SpainReg90.20$qx_male
periods = c(1991:2000)
ages = c(ages)
nPop = 18
lxt = SpainReg90.20$lx_male

multi23 <- fit.multiplicative.LC.multi(qxt = SpainReg90.20$qx_male,
                                       periods = c(1991:2000),
                                       ages = c(ages),
                                       nPop = 18,
                                       lxt = SpainReg90.20$lx_male)
formulti23 <- forecast.multiplicative.LC.multi(fitted.obj = multi23,
                                               nahead = 1,
                                               ktmethod = c("Arimapdq"),
                                               kt_include.cte = TRUE)
multi23
formulti23$qxt.future$pob1

addit23 <- fit.additive.LC.multi(qxt = SpainReg90.20$qx_male,
                                      periods = c(1991:2000),
                                      ages = c(ages),
                                      nPop = 18,
                                      lxt = SpainReg90.20$lx_male)
foraddit23 <- forecast.additive.LC.multi(fitted.obj = addit23,
                                               nahead = 1,
                                               ktmethod = c("Arimapdq"),
                                               kt_include.cte = TRUE)

###########################################################0########
#Voy a intentar replicar los resultados del RISK
loocv_Spainfemales_multi <- multipopulation_loocv(qxt = SpainRegions$qx_female,
                                                model = c("multiplicative"),
                                                periods =  c(1991:2020), ages = c(ages),
                                          nPop = 18, lxt = SpainRegions$lx_female,
                                         ktmethod = c("Arimapdq"),
                                          kt_include.cte = TRUE,
                                         measures = c("SSE"),
                                        trainset1 = 10)
loocv_Spainfemales_multi$kt.future
loocv_Spainfemales_addit <- multipopulation_loocv(qxt = SpainRegions$qx_female,
                                                model = c("additive"),
                                                periods =  c(1991:2020), ages = c(ages),
                                                nPop = 18, lxt = SpainRegions$lx_female,
                                                ktmethod = c("Arimapdq"),
                                                kt_include.cte = TRUE,
                                                measures = c("SSE"),
                                                trainset1 = 10)
loocv_Spainfemales_multi$meas_periodsfut

val <- c(cv_Spainmales_multi$meas_ages, cv_Spainmales_addit$meas_ages,
         loocv_Spainfemales_multi$meas_ages, loocv_Spainfemales_addit$meas_ages)
max1 <- max(val)
min1 <- min(val)
par(mfrow=c(1,2))
plot(ages, cv_Spainmales_multi$meas_ages, type='l', col='black',
     xlab = 'ages', ylab='MAPE', main='Blocked-CV males', lwd = 2, ylim=c(min1, max1))
lines(ages, cv_Spainmales_addit$meas_ages, col ='red', lwd =2)
legend("topleft", lty=c(1), col=c("black", "red"), c("Multiplicative", "Additive"), cex = 1.1, lwd=2)

plot(ages, loocv_Spainfemales_multi$meas_ages, type='l', col='black',
     xlab = 'ages', ylab='MAPE', main='LOOCV females', lwd = 2, ylim=c(min1, max1))
lines(ages, loocv_Spainfemales_addit$meas_ages, col ='red', lwd =2)
legend("topleft", lty=c(1), col=c("black", "red"), c("Multiplicative", "Additive"), cex = 1.1, lwd=2)

val <- c(cv_Spainmales_multi$meas_ages[1:14], cv_Spainmales_addit$meas_ages[1:14],
         loocv_Spainfemales_multi$meas_ages[1:14], loocv_Spainfemales_addit$meas_ages[1:14])
max1 <- max(val)
min1 <- min(val)
par(mfrow = c(1,2))
plot(ages[1:14], cv_Spainmales_multi$meas_ages[1:14], type='l', col='black',
     xlab = 'ages', ylab='MAPE', main='Blocked-CV males', lwd = 2)
lines(ages[1:14], cv_Spainmales_addit$meas_ages[1:14], col ='red', lwd =2)
plot(ages[1:14], loocv_Spainfemales_multi$meas_ages[1:14], type='l', col='black',
     xlab = 'ages', ylab='MAPE', main='LOOCV females', lwd = 2)
lines(ages[1:14], loocv_Spainfemales_addit$meas_ages[1:14], col ='red', lwd =2)



val2 <- c(cv_Spainmales_multi$meas_periodsfut, cv_Spainmales_addit$meas_periodsfut,
           loocv_Spainfemales_multi$meas_periodsfut, loocv_Spainfemales_addit$meas_periodsfut)
periods <- c("1999-2003", "2004-2008", "2009-2013", "2014-2018", "2019-2021")
max2 <- max(val2)
min2 <- min(val2)
par(mfrow=c(1,2))
plot(c(1:5), cv_Spainmales_multi$meas_periodsfut, type='l', col='black', xaxt = "n",
     xlab = 'blocks', ylab='MAPE', main='Blocked-CV males', lwd = 2, ylim = c(min2, max2))
axis(1, 1:5, c("1999-2003", "2004-2008", "2009-2013", "2014-2018", "2019-2021"))
lines(c(1:5), cv_Spainmales_addit$meas_periodsfut, col ='red', lwd =2)
legend("bottomleft", lty=c(1), col=c("black", "red"), c("Multiplicative", "Additive"), cex = 1.1, lwd=2)

plot(c(2001:2020), loocv_Spainfemales_multi$meas_periodsfut, type='l', col='black',
     xlab = 'fut-periods', ylab='MAPE', main='LOOCV females', lwd = 2, ylim = c(min2, max2))
lines(c(2001:2020), loocv_Spainfemales_addit$meas_periodsfut, col ='red', lwd =2)
legend("topleft", lty=c(1), col=c("black", "red"), c("Multiplicative", "Additive"), cex = 1.1, lwd=2)

SpainMap(regionvalue = cv_Spainmales_multi$meas_pop[2:18],
         main = "Blocked-CV males Multiplicative",
         name = "MAPE")
SpainMap(regionvalue = cv_Spainmales_addit$meas_pop[2:18],
         main = "Blocked-CV males Additive",
         name = "MAPE")

SpainMap(regionvalue = loocv_Spainfemales_multi$meas_pop[2:18],
         main = "LOOCV females Multiplicative",
         name = "MAPE")
SpainMap(regionvalue = loocv_Spainfemales_addit$meas_pop[2:18],
         main = "LOOCV females Additive",
         name = "MAPE")

c(cv_Spainmales_multi$meas_total, cv_Spainmales_addit$meas_total,
  loocv_Spainfemales_multi$meas_total, loocv_Spainfemales_addit$meas_total)
