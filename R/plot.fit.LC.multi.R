#' FUNCTION TO PLOT MULTI-POPULATION MORTALITY MODEL
#' @description
#' R function to plot the parameters for the Additive (Debon et al., 2011) and Multiplicative (Russolillo et al., 2011) Multi-Populationn mortality model.
#' It should be mentioned that in case that this function is developed for fitting several populations.
#' However, in case you only consider one population, the function will fit the one-population Lee-Carter model (Lee and Carter, 1992).
#'
#' @param fitted.obj object developed using function fit.additive.LC.multi() and fit.multiplicative.LC.multi()
#'
#' @return plot the different parameters for the multi-population mortality models `ax`, `bx`, `kt` and `Ii`. This function is valid for both approaches Additive and Multiplicative multi-population mortality models.
#'
#' @seealso \code{\link{fit.additive.LC.multi}}, \code{\link{fit.multiplicative.LC.multi}},
#' \code{\link{forecast.additive.LC.multi}}, \code{\link{forecast.multiplicative.LC.multi}},
#' \code{\link{multipopulation_cv}}, \code{\link{multipopulation_loocv}}
#'
#' @references'
#' Debon, A., Montes, F., & Martiez-Ruiz, F. (2011).
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
#' @examples
#' SpainRegions
#' ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)
#' #In this case, we fit for males providing the lxt
#' multiplicative_Spainmales <- fit.multiplicative.LC.multi(qxt = SpainRegions$qx_male,
#'                               periods = c(1991:2020),
#'                               ages = c(ages),
#'                               nPop = 18,
#'                               lxt = SpainRegions$lx_male)
#'
#' #Once, we have fit the data, it is possible to see the ax, bx, kt, and Ii
#' #provided parameters for the fitting.
#' plot.fit.LC.multi(multiplicative_Spainmales)
#'
#' #Equal to the previous step but in this case for females and without
#' #providing lxt.
#' multiplicative_Spainfemales <- fit.multiplicative.LC.multi(qxt = SpainRegions$qx_female,
#'                               periods = c(1991:2020),
#'                               ages = c(ages),
#'                               nPop = 18)
#'
#' #Once, we have fit the data, it is possible to see the ax, bx, kt
#' #provided parameters for the fitting.
#' plot.fit.LC.multi(multiplicative_Spainmales)
#'
#' @export
plot.fit.LC.multi <- function(fitted.obj){
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
    par(mfrow=c(1,3))
    plot(ages, ax, ylab="", xlab="x = age", main =ax_main,
         type="l", lwd=2)
    plot(ages, bx, ylab="", xlab="x = age", main =bx_main,
         type="l", lwd=2)
    plot(pers, kt, ylab="", xlab="t = period", main =kt_main,
         type="l", lwd=2)
  }


}
