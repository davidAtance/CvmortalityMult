#' Spain Regions Mortality data
#'
#' Data from the Spanish region of Spain from the Spanish National Institute of Statistics (INE) for both genders years 1991-2020 and abridged ages from 0 to 90.
#' This dataset contains mortality rates (qxt) from 18 different regions of Spain.
#' Additionally, the dataset includes the number of people alive (lxt) for each age and period.
#'
#' @name SpainRegions
#'
#' @format A data frame with 10800 rows and 9 columns with class \code{"CVmortalityData"} including the following information
#' * `ccaa` a vector containing all the regions of Spain. Indeed, the column takes the following information: Spain, Andalucia, Aragon, Asturias, Baleares, Canarias, Cantabria, Castillayla Mancha, CastillayLeon, Cataluna, ComunidadValenciana, Extremadura, Galicia, Madrid, Murcia, Navarra, PaisVasco, and LaRioja.
#' * `years` a vector containing the periods of the dataset from 1991 to 2020.
#' * `ages` a vector containing the abridged ages considered in the dataset, 0, <1, 1-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, 80-84, 85-89, and 90-94.
#' * `qx_male` mortality rates for the males in every region of Spain including Nation data.
#' * `qx_female` mortality rates for the females in every region of Spain including Nation data.
#' * `lx_male` survivor function considered for the males in every region of Spain including Nation data.
#' * `lx_female` survivor function considered for the females in every region of Spain including Nation data.
#' * `series` information for the series of data provided.
#' * `label` the assigned tag to the data frame.
#'
#' @docType data
#'
#' @usage SpainRegions
#'
#' @references
#' Spanish National Institute of Statistics (INE) (2023). Tablas de mortalidad, metodologia.
#' Technical report, Instituto Nacional de Estadistica
#'
#' @examples
#' #The example takes more than 5 seconds because it includes
#' #several fitting and forecasting process and hence all
#' #the process is included in donttest
#' \donttest{
#' #In this case, we show the region dataset applying it to a multipopulation model.
#' #First, we present the dataset
#' SpainRegions
#'
#' #An example to how the additive multi-population model fits the data
#' ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)
#'
#' library(gnm)
#' multiplicative_Spainmales <- fitLCmulti(model = "multiplicative",
#'                                         qxt = SpainRegions$qx_male,
#'                                         periods = c(1991:2020),
#'                                         ages = c(ages),
#'                                         nPop = 18,
#'                                         lxt = SpainRegions$lx_male)
#'
#' multiplicative_Spainmales
#' }
"SpainRegions"
print.SpainRegionsData <- function(x) {
  cat("Spain Regions data\n")
  cat("geometry to contruct maps")
}
