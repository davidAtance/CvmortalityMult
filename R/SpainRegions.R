#' Spain Regions Mortality data
#'
#' Data from the spanish region of Spain from the Spanish National Institute of Statistics (INE) for both genders years 1991-2020 and abridged ages from 0 to 90.
#' This dataset contains mortality rates (qxt) from 18 different regions of Spain.
#' Additionally, the dataset includes the number of people alive (lxt) for each age and period.
#'
#' @name SpainRegions
#'
#' @docType data
#' @usage SpainRegions
#'
#' @references
#'
#' Spanish National Institute of Statistics (INE) (2023). Tablas de mortalidad, metodologia.
#' Technical report, Instituto Nacional de Estadistica
#'
#' @examples
#' SpainRegions
#' multiplicative_Spainmales <- fit.multiplicative.LC.multi(qxt = SpainRegions$qx_male,
#'                                                 periods = c(1991:2020),
#'                                                 ages = c(ages),
#'                                                 nPop = 18,
#'                                                 lxt = SpainRegions$lx_male)
#' @export
#' SpainRegions
