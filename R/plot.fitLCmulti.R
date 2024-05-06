#' Function to plot the parameters of the multi-population mortality models
#' @description
#' R function to plot the parameters for the Additive (Debon et al., 2011) and Multiplicative (Russolillo et al., 2011) Multi-Population mortality model.
#' It should be mentioned that in case that this function is developed for fitting several populations.
#' However, in case you only consider one population, the function will fit the one-population Lee-Carter model (Lee and Carter, 1992).
#'
#' @param fitted.obj object developed using function `fitLCmulti()` which are objects of the `fitLCmulti` class.
#'
#' @return plot the different parameters for the multi-population mortality models `ax`, `bx`, `kt` and `Ii`. This function is valid for both approaches Additive and Multiplicative multi-population mortality models.
#'
#' @seealso \code{\link{fitLCmulti}}, \code{\link{forecast.fitLCmulti}},
#' \code{\link{plot.forLCmult}},
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
#' @importFrom graphics par
#' @importFrom utils install.packages
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
#' multiplicative_Spainmales
#'
#' #Once, we have fit the data, it is possible to see the ax, bx, kt, and It
#' #provided parameters for the fitting.
#' plot(multiplicative_Spainmales)
#'
#' #LEE-CARTER FOR SINGLE-POPULATION
#' #As we mentioned in the details of the function, if we only provide the data
#' #from one-population the function fit_additive.LC.multi()
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
#' }
#' @export
plot.fitLCmulti <- function(fitted.obj){
  pers <- fitted.obj$Periods
  ages <- fitted.obj$Ages
  pops <- fitted.obj$nPop
  ax <- fitted.obj$ax
  bx <- fitted.obj$bx
  kt <- fitted.obj$kt
  Ii <- fitted.obj$Ii

  ax_main <- expression(a[x])
  bx_main <- expression(b[x])
  kt_main <- expression(k[t])
  Ii_main <- expression(I[i])


  if(pops != 1){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    par(mfrow=c(1,4))
    plot(ages, ax, ylab="", xlab="x = age", main =ax_main,
         type="l", lwd=2)
    plot(ages, bx, ylab="", xlab="x = age", main =bx_main,
         type="l", lwd=2)
    plot(pers, kt, ylab="", xlab="t = period", main =kt_main,
         type="l", lwd=2)
    plot(c(1:fitted.obj$nPop), Ii, ylab="", xlab="i = population", main = Ii_main,
         type="l", lwd=2)

  } else if(pops == 1){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    par(mfrow=c(1,3))
    plot(ages, ax, ylab="", xlab="x = age", main =ax_main,
         type="l", lwd=2)
    plot(ages, bx, ylab="", xlab="x = age", main =bx_main,
         type="l", lwd=2)
    plot(pers, kt, ylab="", xlab="t = period", main =kt_main,
         type="l", lwd=2)
    }

}
