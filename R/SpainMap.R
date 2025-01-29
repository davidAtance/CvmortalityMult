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
#' @importFrom sf st_sf st_crs st_transform
#'
#' @examples
#' name <- c("Ii")
#' main <- c("Multiplicative for males")
#' regionvalue <- c(0.131867619, -0.063994652,  0.088094096,
#'                  0.019685552,  0.064671498,   0.012212161,
#'                 -0.088864474, -0.146079884, -0.017703536,
#'                  0.050376721,  0.052476852, -0.022871202,
#'                 -0.093952332,  0.049266816, -0.101224890,
#'                  0.001481333, -0.078523511)
#' library(sf)
#'
#' SpainMap(regionvalue, main, name)
#'
#' @export
SpainMap <- function(regionvalue, main, name){

  if(length(regionvalue) != 17){
    stop("The regionvalue is not a vector of length 17.")
  }

  autonomias <- st_sf(CvmortalityMult::regions)
  autonomias$Ii <- regionvalue
  names(autonomias)[5] <- name
  sf::st_crs(autonomias) <- 32630
  autonomias <- st_transform(autonomias, crs = 32630)

  #if (is.na(st_crs(autonomias))) {
  #  st_crs(autonomias) <- 32630
  #} else {
  #  # Si ya tiene un CRS asignado, y no es el correcto, cambiar a 32630
  #  if (st_crs(autonomias)$epsg != 32630) {
  #    autonomias <- st_transform(autonomias, crs = 32630)
  #  }
  #}

  tm_shape(autonomias) +
    tm_polygons(col=name, palette = c("#FFFFAF", "#FFD800", "#FF9D00", "#FF0000"),
                  breaks=quantile(autonomias[[5]]), border.col = "black",
                  title = name) +
    tm_layout(main,
                legend.title.fontfamily = "serif",
                legend.title.size = 1.3,
                legend.text.size = 0.9,
                legend.width = 0.4) +
    tm_shape(autonomias) +
    tm_borders(lwd = 1)



}
