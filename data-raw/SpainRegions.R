## code to prepare `DATASET` dataset goes here
load("datos.brutos.ccaa5.RData")
datos.final.ccaa5 <- datos.brutos.ccaa5[(datos.brutos.ccaa5$edad<95),]

Spain_Regions <- datos.final.ccaa5
nombre<-c("Spain", "Andalucia", "Aragon", "Asturias", "Baleares", "Canarias", "Cantabria", "Castillayla Mancha",
          "CastillayLeon", "Cataluna", "ComunidadValenciana", "Extremadura", "Galicia", "Madrid", "Murcia",
          "Navarra", "PaisVasco", "LaRioja")
i <- 1
for(i in 1:18){
  Spain_Regions <- within(Spain_Regions, ccaa[ccaa == (i-1)] <- nombre[i])
}

Spain_Regions$tiempo <- NULL
Spain_Nat <- Spain_Regions[Spain_Regions$ccaa == "Spain",]

Spain_Regions <- as.data.frame(Spain_Regions)
SpainRegions <- Spain_Regions

SpainNat <- Spain_Nat

SpainRegions <- structure(list(ccaa = Spain_Regions$ccaa,
                                years = Spain_Regions$periodo,
                                ages = Spain_Regions$edad,
                                qx_male = Spain_Regions$homqx,
                                qx_female = Spain_Regions$mujqx,
                                lx_male = Spain_Regions$homriesgo,
                                lx_female = Spain_Regions$mujriesgo,
                                series = c('males and females'),
                                label = c('Spain regions')))
print.CVmortalityData <- function(x) {
  cat("Mortality Data\n")
  cat(x$label, "including", x$series ,"\n")
  cat("Years", c(min(x$years),":", max(x$years)),"\n")
  cat("Abridged Ages", c(min(x$ages),":", max(x$ages)), "\n")
}

class(SpainRegions) <- "CVmortalityData"
SpainRegions
save(SpainRegions, file = "data/SpainRegions.rda")

SpainNat <- structure(list(ccaa = Spain_Nat$ccaa,
                               years = Spain_Nat$periodo,
                               ages = Spain_Nat$edad,
                               qx_male = Spain_Nat$homqx,
                               qx_female = Spain_Nat$mujqx,
                               lx_male = Spain_Nat$homriesgo,
                               lx_female = Spain_Nat$homriesgo,
                               lx_female = Spain_Nat$mujriesgo,
                               series = c('males and females'),
                               label = c('Spain Total population')))
class(SpainNat) <- "CVmortalityData"
SpainNat
SpainRegions
save(SpainNat, file = "data/SpainNat.rda")

regions <- load(file = "data/regions.rda")
regions <- structure(regions,
                     class = "SpainRegionsData")

print.SpainRegionsData <- function(x) {
  cat("Spain Regions data\n")
  cat("geometry to contruct maps")
}
class(regions) <- "SpainRegionsData"
regions
SpainNat
SpainRegions
save(regions, file = "data/regions.rda")

#include in the file NAMESPACE:
#S3method(print,CVmortalityData)
#S3method(print,SpainRegionsData)

roxygenise()
devtools::build()
devtools::document()


