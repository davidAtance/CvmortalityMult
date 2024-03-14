#' Spain National map information
#' @description
#' This data contains information to plot the percentiles plot in Spanish regions. Therefore, the users only have to provide a specific variable to show in regions of Spain.
#' @param regionvalue vector with the values that you want to plot in percentiles in the Spain map.
#' @param main the specific title of the map plot
#' @param name the assigned name for the legend in map plot.
#'
#' @return a map from the regions of Spain colored with the variable provided by the user.
#'
#' @references
#' Spanish National Institute of Statistics (INE) (2023). Tablas de mortalidad, metodologia.
#' Technical report, Instituto Nacional de Estadistica
#'
#' @importFrom tmap tm_shape tm_polygons tm_layout tm_borders
#' @importFrom stats quantile
#' @importFrom utils install.packages
#'
#' @examples
#' ages <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40,
#' 45, 50, 55, 60, 65, 70, 75, 80, 85, 90)
#' #We fit the additive multi-population mortality model for spain males
#' additive_Spainmales <- fit_additive.LC.multi(qxt = SpainRegions$qx_male,
#'                               periods = c(1991:2020),
#'                               ages = c(ages),
#'                               nPop = 18,
#'                               lxt = SpainRegions$lx_male)
#' name <- c("Ii")
#' main <- c("Multiplicative for males")
#' regionvalue <- additive_Spainmales$Ii[2:18]
#'
#' SpainMap(regionvalue, main, name)
#'
#' @export
SpainMap <- function(regionvalue, main, name){

  if(length(regionvalue) != 17){
    stop("The regionvalue is not a vector of length 17.")
  }
  regions <- NULL
  autonomias <- regions
  autonomias$Ii <- regionvalue
  names(autonomias)[4] <- name

  tm_shape(autonomias) +
    tm_polygons(col=name, palette = c("#FF0000", "#FF9D00", "#FFD800", "#FFFFAF"),
                breaks=quantile(autonomias[[4]]), border.col = "black",
                title = name) +
    tm_layout(main,
              legend.title.fontfamily = "serif",
              legend.title.size = 1.3,
              legend.text.size = 0.9) +
    tm_borders(lwd = 2)

}
