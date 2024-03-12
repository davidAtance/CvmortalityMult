#' Spain National Mortality data
#'
#' Data from the spanish national of Spain from the Spanish National Institute of Statistics (INE) for both genders years 1991-2020 and abridged ages from 0 to 90.
#' This dataset contains mortality rates for the total national population of Spain.
#' Additionally, the dataset includes the number of people alive (lxt) for each age and period.
#'
#' @name SpainNat
#'
#' @format A data frame with 600 rows and 9 columns with the following information
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
#' SpainNat
#' LC_Spainmales <- fit.additive.LC.multi(qxt = SpainNat$qx_male,
#'                                        periods = c(1991:2020),
#'                                        ages = ages,
#'                                        nPop = 1)
#' plot.fit.LC.multi(LC_Spainmales)
#'
"SpainNat"
