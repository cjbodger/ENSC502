# Relative Humidity Approximation from Dewpoint Temp and Surface Temp
# Using Magnus Approximation
# Temperature (T) and Dewpoint Temperature (Td) must be in degrees celcius
library(terra)
library(climetrics)
library(maps)
library(vars)

setwd("C:/Users/Charlotte/Documents/ENSC502")

T <- terra::rast("C:/Users/Charlotte/Documents/ENSC502/F2MT_2001_2022_K.nc")
Td <- terra::rast("C:/Users/Charlotte/Documents/ENSC502/FDPT_2001_2022_K.nc")

Tc <- T - 273.15
Tdc <- Td - 273.15

numer <- exp((17.625*Tdc)/(243.04+Tdc))
denom <- exp((17.625*Tc)/(243.04+Tc))
RH <- 100 * (numer/denom)
print(RH)

writeCDF(RH,
         "FRH_2001_2022.nc",
         varname = "RH",
         varunit = "%",
         longname = "Relative Humidity",
         overwrite = TRUE)
