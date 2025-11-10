# Relative Humidity Approximation from Dewpoint Temp and Surface Temp
# Using Magnus Approximation
# Temperature (T) and Dewpoint Temperature (Td) must be in degrees celcius
library(terra)
library(climetrics)
library(maps)
library(vars)

setwd("C:/Users/Charlotte/Desktop/RStudio")

T <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_tas_C.nc")
Td <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_dewpoint_C.nc")

numer <- exp((17.625*Td)/(243.04+Td))
denom <- exp((17.625*T)/(243.04+T))
RH <- 100 * (numer/denom)
print(RH)

writeCDF(RH,
         "ERA5_rh.nc",
         varname = "RH",
         varunit = "%",
         longname = "Relative Humidity",
         overwrite = TRUE)