#' Spain National Mortality data
#'
#' Data from the Spanish national of Spain from the Spanish National Institute of Statistics (INE) for both genders years 1991-2020 and abridged ages from 0 to 90.
#' This dataset contains mortality rates for the total national population of Spain.
#' Additionally, the dataset includes the number of people alive (lxt) for each age and period.
#'
#' @name SpainNat
#'
#' @format A data frame with 600 rows and 9 columns with class \code{"CVmortalityData"} including the following information
#' * `ccaa` a vector containing all the regions of Spain. Indeed, the column takes the following information: Spain.
#' * `years` a vector containing the periods of the dataset from 1991 to 2020.
#' * `ages` a vector containing the abridged ages considered in the dataset, 0, <1, 1-4, 5-9, 10-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, 45-49, 50-54, 55-59, 60-64, 65-69, 70-74, 75-79, 80-84, 85-89, and 90-94.
#' * `qx_male` mortality rates for the males in the Spain Nation.
#' * `qx_female` mortality rates for the females in the Spain Nation.
#' * `lx_male` survivor function considered for the males of Spain Nation.
#' * `lx_female` survivor function considered for the females of Spain Nation.
#' * `series` information for the series of data provided.
#' * `label` the assigned tag to the data frame.
#'
#' @docType data
#' @usage SpainNat
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
#' SpainNat
#' #An example to how the additive multi-population model fits the data
#' ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)
#' library(gnm)
#' LC_Spainmales <- fitLCmulti(qxt = SpainNat$qx_male,
#'                             periods = c(1991:2020),
#'                             ages = ages,
#'                             nPop = 1)
#'
#' LC_Spainmales
#' }
"SpainNat"
#' @export
print.CVmortalityData <- function(x, ...) {
  cat("Mortality Data\n")
  cat(x$label, "including", x$series ,"\n")
  cat("Periods", c(min(x$periods),":", max(x$periods)),"\n")
  cat("Abridged Ages", c(min(x$ages),":", max(x$ages)), "\n")
}
