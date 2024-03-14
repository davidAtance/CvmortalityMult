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
Spain_Regions <- structure(list(ccaa = Spain_Regions$ccaa,
                                years = Spain_Regions$periodo,
                                ages = Spain_Regions$edad,
                                qx_male = Spain_Regions$homqx,
                                qx_female = Spain_Regions$mujqx,
                                lx_male = Spain_Regions$homriesgo,
                                lx_female = Spain_Regions$mujriesgo,
                                series = c('male and female'),
                                label = c('Spain regions')),
                           class = "myclass")
print.myclass <- function(x) {
  cat("Mortality Data\n")
  cat("Spain Regions for males and females\n")
  cat("Years", c(min(x$years),":", max(x$years)),"\n")
  cat("Abridged Ages", c(min(x$ages),":", max(x$ages)), "\n")
}
SpainRegions <- Spain_Regions
SpainRegions

Spain_Nat <- as.data.frame(Spain_Nat)
Spain_Nat <- structure(list(ccaa = Spain_Nat$ccaa,
                            years = Spain_Nat$periodo,
                            ages = Spain_Nat$edad,
                            qx_male = Spain_Nat$homqx,
                            qx_female = Spain_Nat$mujqx,
                            lx_male = Spain_Nat$homriesgo,
                            lx_female = Spain_Nat$mujriesgo,
                            series = c('male and female'),
                            label = c('Spain National population')),
                       class = "myclass1")
print.myclass1 <- function(x) {
  cat("Mortality Data\n")
  cat("Spain National population for males and females\n")
  cat("Years", c(min(x$years),":", max(x$years)),"\n")
  cat("Abridged Ages", c(min(x$ages),":", max(x$ages)), "\n")
}


SpainNat <- Spain_Nat
SpainNat
SpainRegions

save(SpainNat, file = "data/SpainNat.RData")
save(SpainRegions, file = "data/SpainRegions.RData")

usethis::use_data_raw()
#regions <- autonomias

save(regions, file = "regions.RData")
regions <- structure(regions,
                     class = "myclass2")

print.myclass2 <- function(x) {
  cat("Spain Regions\n")
  cat("geometry to contruct maps")
}
regions
SpainNat
SpainRegions

save(regions, file = "data/regions.RData")


print.myclass(autonomias)
autonomias
usethis::use_data(regions, overwrite = TRUE)

usethis::use_data(SpainRegions, SpainNat, overwrite = TRUE)
library(devtools)
use_data(SpainRegions, SpainNat, overwrite = TRUE)
