#' Function to plot the parameters of the multi-population mortality models
#' @description
#' R function to plot the parameters for the Additive (Debon et al., 2011) and Multiplicative (Russolillo et al., 2011) Multi-Population mortality model.
#' It should be mentioned that in case that this function is developed for fitting several populations.
#' However, in case you only consider one population, the function will fit the one-population Lee-Carter model (Lee and Carter, 1992).
#'
#' @param x `x` developed using function `fitLCmulti()` which are objects of the `fitLCmulti` class.
#' @param ... additional arguments to show in the plot appearance.
#'
#' @return plot the different parameters for the multi-population mortality models `ax`, `bx`, `kt` and `Ii`. This function is valid for both approaches Additive and Multiplicative multi-population mortality models.
#'
#' @seealso \code{\link{fitLCmulti}}, \code{\link{forecast.fitLCmulti}},
#' \code{\link{plot.forLCmulti}}, \code{\link{multipopulation_cv}}
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
#' @importFrom stats plogis qlogis
#' @importFrom graphics layout
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
#' #' #3. COMMON-FACTOR MULTI-POPULATION MORTALITY MODEL
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
#' #5. JOINT-K MULTI-POPULATION MORTALITY MODEL
#' #In the case, the user wants to fit the augmented-common-factor multi-population mortality model
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
#' }
#' @export
plot.fitLCmulti <- function(x, ...){
  if(!is.null(x)){
    if(!"fitLCmulti" %in% class(x))
      stop("The 'x' does not have the 'fitLCmulti' structure of R CvmortalityMult package.")
  }
  if(!is.list(x)){
    stop("The 'x' is not a list. Use 'fitLCmulti' function first.")
  }

  pers <- x$Periods
  ages <- x$Ages
  pops <- x$nPop
  ax <- x$ax
  bx <- x$bx
  kt <- x$kt
  Ii <- x$Ii

  ax_main <- expression(a[x])
  bx_main <- expression(b[x])
  kt_main <- expression(k[t])
  Ii_main <- expression(I[i])


  if(pops != 1){
    if(x$model == "additive" | x$model == "multiplicative"){
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))

      par(mfrow=c(1,4))
      plot(ages, ax, ylab="", xlab="x = age", main =ax_main,
           type="l", lwd=2)
      plot(ages, bx, ylab="x", xlab="x = age", main =bx_main,
           type="l", lwd=2)
      plot(pers, kt, ylab="", xlab="t = period", main =kt_main,
           type="l", lwd=2)
      plot(c(1:x$nPop), Ii, ylab="", xlab="i = population", main = Ii_main,
           type="l", lwd=2)
    }else if(x$model == "CFM"){
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))

      layout(matrix(c(1,2,3, 4, 4, 4), ncol = 3, byrow = T),
             heights = c(0.9, 0.1))

      min1 <- min(ax)
      max1 <- max(ax)
      plot(ages, ax[1,], ylab="", xlab="x = age", main =ax_main,
           type="l", lty=1, ylim = c(min1, max1))
      namepop <- c("Pop1")
      for(pe in 2:pops){
        lines(ages, ax[pe,], lty = pe)
        namepop <- c(namepop, paste0("Pop", pe))
      }
      #legend("topleft", c(namepop), lty = c(1:pops), cex = 0.5)
      plot(ages, bx, ylab="x", xlab="x = age", main =bx_main,
           type="l")
      plot(pers, kt, ylab="", xlab="t = period", main =kt_main,
           type="l")
      namepop <- c("Pop1", "Pop2")
      if(pops >=3){
        for(pe in 3:pops){
          namepop <- c(namepop, paste0("Pop", pe))
        }}
      if(pops %% 2 != 0){
        pops2 <- pops + 1
        names_pops <- c(namepop, "NA")
        pop_d <- pops2/2
      } else{
        pops2 <- pops
        names_pops <- c(namepop)
        pop_d <- pops2/2
      }
      legend_order <- matrix(1:pops2, ncol = pop_d, byrow = T)
      par(mar=c(0,0,0,0))
      plot(1, type = "n", axes=F, xlab="", ylab="")
      legend("top", c(names_pops)[legend_order],
             lty = c(1:pops2)[legend_order],
             ncol = pop_d, cex = 0.6)


    }else if(x$model == "ACFM"){
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))

      layout(matrix(c(1,2,3,4,5,6, 7, 7, 7), ncol = 3, byrow = T),
             heights = c(0.45,0.45,0.1))

      plot(ages, ax[1,], ylab="", xlab="x = age", main =ax_main,
           type="l", lty=1)

      plot(ages, bx[1,], ylab="", xlab="x = age", main =bx_main,
           type="l")
      plot(pers, kt[,1], ylab="", xlab="t = period", main =kt_main,
           type="l")

      namepop <- c("Pop1", "Pop2")
      max1 <- max(ax[2:pops,])
      min1 <- min(ax[2:pops,])
      plot(ages, ax[2,], ylab="", xlab="x = age", main =expression(a[x[i]]),
           type="l", ylim = c(min1, max1), lty = 2)
      if(pops >=3){
        for(pe in 3:pops){
          lines(ages, ax[pe,], lty = pe)
          namepop <- c(namepop, paste0("Pop", pe))
        }}

      max2 <- max(bx[2:pops,])
      min2 <- min(bx[2:pops,])
      plot(ages, bx[2,], ylab="", xlab="x = age", main =expression(b[x[i]]),
           type="l", ylim = c(min2, max2), lty = 2)
      if(pops >=3){
        for(pe in 3:pops){
          lines(ages, bx[pe,], lty = pe)
      }}

      max3 <- max(kt[,2:pops])
      min3 <- min(kt[,2:pops])
      plot(pers, kt[,2], ylab="", xlab="t = period", main =expression(k[x[i]]),
           type="l", ylim = c(min3, max3), lty = 2)
      if(pops >=3){
        for(pe in 2:pops){
          lines(pers, kt[,pe], lty = pe)
        }}
      if(pops %% 2 != 0){
        pops2 <- pops + 1
        names_pops <- c(namepop, "NA")
        pop_d <- pops2/2
      } else{
        pops2 <- pops
        names_pops <- c(namepop)
        pop_d <- pops2/2
      }
      legend_order <- matrix(1:pops2, ncol = pop_d, byrow = T)
      par(mar=c(0,0,0,0))
      plot(1, type = "n", axes=F, xlab="", ylab="")
      legend("top", c(names_pops)[legend_order],
             lty = c(1:pops2)[legend_order],
             ncol = pop_d, cex = 0.6)

    }else if(x$model == "joint-K"){
      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar))

      layout(matrix(c(1,2,3, 4, 4, 4), ncol = 3, byrow = T),
             heights = c(0.9, 0.1))

      min1 <- min(ax)
      max1 <- max(ax)
      plot(ages, ax[1,], ylab="", xlab="x = age", main =ax_main,
           type="l", lty=1, ylim = c(min1, max1))
      namepop <- c("Pop1")
      for(pe in 2:pops){
        lines(ages, ax[pe,], lty = pe)
        namepop <- c(namepop, paste0("Pop", pe))
      }

      min2 <- min(bx)
      max2 <- max(bx)
      plot(ages, bx[1,], ylab="", xlab="x = age", main =bx_main,
           type="l", ylim = c(min2, max2))
      for(pe in 2:pops){
        lines(ages, bx[pe,], lty = pe)
      }

      plot(pers, kt, ylab="", xlab="t = period", main =kt_main,
           type="l")
      namepop <- c("Pop1", "Pop2")
      if(pops >=3){
        for(pe in 3:pops){
          namepop <- c(namepop, paste0("Pop", pe))
        }}
      if(pops %% 2 != 0){
        pops2 <- pops + 1
        names_pops <- c(namepop, "NA")
        pop_d <- pops2/2
      } else{
        pops2 <- pops
        names_pops <- c(namepop)
        pop_d <- pops2/2
      }
      legend_order <- matrix(1:pops2, ncol = pop_d, byrow = T)
      par(mar=c(0,0,0,0))
      plot(1, type = "n", axes=F, xlab="", ylab="")
      legend("top", c(names_pops)[legend_order],
             lty = c(1:pops2)[legend_order],
             ncol = pop_d, cex = 0.6)


    }

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
