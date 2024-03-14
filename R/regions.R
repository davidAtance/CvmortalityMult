#' regions
#'
#' Data from the spanish region of Spain which are provided to plot an indicator.
#' This dataset contains a plot with the information of Spain regions (geometry and name of every region).
#'
#' @name regions
#'
#' @format A data frame with 600 rows and 9 columns with the following information
#' * `Codigo` a vector containing the code of every region of Spain.
#' * `Texto` a vector containing the name of every region of Spain.
#' * `Texto_Alt` a vector containing the long name of every region of Spain.
#' * `Ii` a vector containing a possible value of one indicator to be shown.
#' * `geometry` the dimension of every region of Spain. This vector allows to plot the regions of Spain.
#'
#' @docType data
#' @return a plot with the Spain regions colored by the indicator provided.
#'
#' @references
#' Spanish National Institute of Statistics (INE) (2023). Tablas de mortalidad, metodologia.
#' Technical report, Instituto Nacional de Estadistica
#'
#' @examples
#' regions
#' multiplicative_Spainmales <- fit_multiplicative.LC.multi(qxt = SpainRegions$qx_male,
#'                                                 periods = c(1991:2020),
#'                                                 ages = c(ages),
#'                                                 nPop = 18,
#'                                                 lxt = SpainRegions$lx_male)
#'
#' SpainMap(regionvalue = multiplicative_Spainmales$Ii[2:18],
#'          main = c("Multiplicative for males"),
#'          name = c("Ii"))
#'
"regions"
