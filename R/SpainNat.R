#' Spain National Mortality data
#'
#' Data from the spanish national of Spain from the Spanish National Institute of Statistics (INE) for both genders years 1991-2020 and abridged ages from 0 to 90.
#' This dataset contains mortality rates for the total national population of Spain.
#' Additionally, the dataset includes the number of people alive (lxt) for each age and period.
#'
#' @name SpainNat
#'
#' @docType data
#' @usage SpainNat
#'
#' @references
#'
#' Spanish National Institute of Statistics (INE) (2023). Tablas de mortalidad, metodologia.
#' Technical report, Instituto Nacional de Estadistica
#'
#' @examples
#' SpainNat
#' ?SpainNat
#' LC_Spainmales <- fit.additive.LC.multi(qxt = SpainNat$qx_male,
#'                                        periods = c(1991:2020),
#'                                        ages = ages,
#'                                        nPop = 1)
#' plot.fit.additive.LC.multi(LC_Spainmales)
#' @export
#' SpainNat
