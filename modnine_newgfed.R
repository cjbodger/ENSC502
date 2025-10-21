# Module 9 Replay

library(terra)
library(climetrics)
library(maps)
library(vars)


# GFED
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/burntFrac_fix.nc")
crs(data) <- "EPSG:4326"
gfed <- data
mask <- mean(gfed) # create mask
mask <- mask - mask + 1
rm(data)

# Precipitation
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_pr_fix.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
pr <- data
pr <- data * 1000 # convert unit from m to mm
rm(data)

# Near Surface Temp
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_tas_fix.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
tas <- data
rm(data)

# Surface Downwelling
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_rsds_fix.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
rsds <- data / 86400 # Convert unit from J m-2 day-1 to W m-2 
rm(data)

# Wind Speed @ 10m
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_wind_fix.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
wind <- data
rm(data)

# Volumetric Soil Mositure Layer 1 (0-7cm)
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_VSM_fix.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
vsm <- data
rm(data)

# Leaf Area Index
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_lai_fix.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
lai <- data
rm(data)

data <- pr
# Create a sequence of dates
start_date <- as.Date("2001-01-01")
end_date <- as.Date("2015-12-01")
dates <- seq(from = start_date, to = end_date, by = "month")

data <- pr
# Raster Time Series
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts,'mean')
my.col <- rev(map.pal("magma", n = 100))
plot(data.12, col = my.col)

# Anomalies
n <- length(dates)/12
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
my.col <- rev(map.pal("magma", n = 100))
plot(subset(data.anom, 1:1), col = my.col)
map("world2", add = TRUE)


# Detrending Anomalies

# Removing Linearity
detrend.fun <- function(x) {
  time <- 1:length(x)
  # If a grid cell contains NA, then set the result to NA
  if (is.na(mean(x))) {
    return(rep(NA, length(time)))
  } else {
    time <- 1:length(x)
    linear.model <- lm(x ~ time)
    detrended.series <- stats::residuals(linear.model)
    return(detrended.series)
  }
}


data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
pr.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(pr.anom.detrend, 1:1), col = my.col, main ="pr")


# Repetition for TAS
data <- tas
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
tas.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(tas.anom.detrend, 1:1), col = my.col, main = "tas")

# Repetition for DSR
data <- rsds
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
rsds.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(rsds.anom.detrend, 1:1), col = my.col, main = "rsds")


# Repetition for GFED
data <- gfed
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
gfed.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("ryg", n = 100))
plot(subset(gfed.anom.detrend, 1), col = my.col, main = "gfed")


# Rep for wind
data <- wind
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
wind.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(wind.anom.detrend, 1), col = my.col, main = "wind speed")


# Rep for VSM
data <- vsm
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
vsm.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(vsm.anom.detrend, 1), col = my.col, main = "volumetric soil moisture layer 1")

# Rep for LAI
data <- lai
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
lai.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(vsm.anom.detrend, 1), col = my.col, main = "Leaf Area Index - Low Vegetation")

granger.fun <- function(x) {
  
  # If a grid cell contains NA, then set the result to NA
  if (is.na(mean(x))) {
    return(NA)
  } else {
    
    n <- length(x)
    # Get first (a) and second (b) variable
    a <- x[1:(n/2)]
    b <- x[(n/2+1):n]
    # Convert to time series
    a <- ts(a)
    b <- ts(b)
    tsDat <- ts.union(a, b)
    tsVAR <- vars::VAR(tsDat, p = 2)
    # Apply Granger causality test 
    p.value <- c(vars::causality(tsVAR, cause = "a")$Granger[3]$p.value)
    return(p.value)
  }
}



granger.fun.2 <- function(x) {
  
  # If the cell has missing or constant data, skip
  if (any(is.na(x))) return(NA)
  
  n <- length(x)
  a <- x[1:(n/2)]
  b <- x[(n/2 + 1):n]
  
  # Convert to time series
  a <- ts(a)
  b <- ts(b)
  tsDat <- ts.union(a, b)
  
  # Attempt VAR and Granger causality safely
  p.value <- tryCatch({
    tsVAR <- vars::VAR(tsDat, p = 2)
    causality_res <- vars::causality(tsVAR, cause = "a")
    
    # Return the p-value safely
    as.numeric(causality_res$Granger$p.value)
  },
  error = function(e) {
    # Catch singular matrix, NA, or model fitting errors
    return(NA)
  })
  
  return(p.value)
}


# PR GC GFED
data <- c(pr.anom.detrend, gfed.anom.detrend)

# Messing with two different functions
p.value <- app(x = data, fun = granger.fun.2)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.pr <- p.value
plot(p.value.pr, col = c("white", "orange"), main = "Precipitation Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)

# TAS GC GFED
data <- c(tas.anom.detrend, gfed.anom.detrend)
p.value <- app(x = data, fun = granger.fun.2)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.tas <- p.value
plot(p.value, col = c("white", "red"), main = "Temperature Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)

# RSDS GC GFED
data <- c(rsds.anom.detrend, gfed.anom.detrend)
p.value <- app(x = data, fun = granger.fun.2)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.rsds <- p.value
plot(p.value, col = c("white", "cyan"), main = "Downwelling Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)

# Wind GC GFED
data <- c(wind.anom.detrend, gfed.anom.detrend)
p.value <- app(x = data, fun = granger.fun.2)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.wind <- p.value
plot(p.value, col = c("white", "purple"), main = "Wind Speed Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)

# VSM GC GFED
data <- c(vsm.anom.detrend, gfed.anom.detrend)
p.value <- app(x = data, fun = granger.fun.2)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.vsm <- p.value
plot(p.value, col = c("white", "pink"), main = "VSM Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)

# LAI GC GFED
data <- c(lai.anom.detrend, gfed.anom.detrend)
p.value <- app(x = data, fun = granger.fun.2)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.lai <- p.value
plot(p.value, col = c("white", "green"), main = "LAI Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)

# Summarizing
p.value <- p.value.pr + p.value.tas + p.value.rsds + p.value.wind + p.value.vsm + p.value.lai
plot(p.value, main = "PR, TAS, RSDS, WS, VSM and LAI Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)

# Cross Correlation
lon <- 131
lat <- -17
coords <- matrix(c(lon, lat), ncol = 2, byrow = TRUE)
location <- vect(coords, type = "points")

plot(p.value)
map("world2", add = TRUE, interior = FALSE)
points(location, col = "red", pch = 16, cex = 1.0)

gfed.anom.detrend.gc <- extract(gfed.anom.detrend, location)
pr.anom.detrend.gc <- extract(pr.anom.detrend, location)
tas.anom.detrend.gc <- extract(tas.anom.detrend, location)
rsds.anom.detrend.gc <- extract(rsds.anom.detrend, location)
wind.anom.detrend.gc <- extract(wind.anom.detrend, location)
vsm.anom.detrend.gc <- extract(vsm.anom.detrend, location)
lai.anom.detrend.gc <- extract(lai.anom.detrend, location)

gfed.anom.detrend.gc <- unlist(unname(as.vector(gfed.anom.detrend.gc)))
pr.anom.detrend.gc <- unlist(unname(as.vector(pr.anom.detrend.gc)))
tas.anom.detrend.gc <- unlist(unname(as.vector(tas.anom.detrend.gc)))
rsds.anom.detrend.gc <- unlist(unname(as.vector(rsds.anom.detrend.gc)))
wind.anom.detrend.gc <- unlist(unname(as.vector(wind.anom.detrend.gc)))
vsm.anom.detrend.gc <- unlist(unname(as.vector(vsm.anom.detrend.gc)))
lai.anom.detrend.gc <- unlist(unname(as.vector(lai.anom.detrend.gc)))

n <- length(gfed.anom.detrend.gc)
gfed.anom.detrend.gc <- gfed.anom.detrend.gc[2:n]
pr.anom.detrend.gc <- pr.anom.detrend.gc[2:n]
tas.anom.detrend.gc <- tas.anom.detrend.gc[2:n]
rsds.anom.detrend.gc <- rsds.anom.detrend.gc[2:n]
wind.anom.detrend.gc <- wind.anom.detrend.gc[2:n]
vsm.anom.detrend.gc <- vsm.anom.detrend.gc[2:n]
lai.anom.detrend.gc <- lai.anom.detrend.gc[2:n]


gfed.ts <- ts(gfed.anom.detrend.gc)
pr.ts <- ts(pr.anom.detrend.gc)
tas.ts <- ts(tas.anom.detrend.gc)
rsds.ts <- ts(rsds.anom.detrend.gc)
wind.ts <- ts(wind.anom.detrend.gc)
vsm.ts <- ts(vsm.anom.detrend.gc)
lai.ts <- ts(lai.anom.detrend.gc)

tsDat <- ts.union(gfed.ts, pr.ts, tas.ts, rsds.ts, wind.ts, vsm.ts, lai.ts)

# Time Series Graph
plot(tsDat)

# Cross Correlation Comuptation
ccf_result <- ccf(gfed.ts, pr.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
ccf_result <- ccf(gfed.ts, tas.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
ccf_result <- ccf(gfed.ts, rsds.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
ccf_result <- ccf(gfed.ts, wind.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
ccf_result <- ccf(gfed.ts, vsm.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
ccf_result <- ccf(gfed.ts, lai.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)

# Index of maximum correlation
max_index <- which.max(ccf_result$acf)
max_correlation <- ccf_result$acf[max_index]
corresponding_lag <- ccf_result$lag[max_index]
result <- data.frame(Max_Correlation = max_correlation, Lag = corresponding_lag)
print(result)


# Finding time lag with largest correlation coefficient
findLag.fun <- function(x) {
  # Handle NA / constant series
  if (any(is.na(x))) return(NA_real_)
  
  n <- length(x)
  a <- x[1:(n/2)]
  b <- x[(n/2 + 1):n]
  
  if (sd(a) == 0 || sd(b) == 0) return(NA_real_)
  
  a <- ts(as.numeric(a))
  b <- ts(as.numeric(b))
  
  # Compute cross-correlation safely
  ccf_result <- stats::ccf(a, b, lag.max = 12, plot = FALSE)
  
  # Extract numeric vectors explicitly
  acf_vals <- as.numeric(ccf_result$acf)
  lag_vals <- as.numeric(ccf_result$lag)
  
  # Find index of max absolute correlation
  max_index <- which.max(abs(acf_vals))
  
  # Extract lag as a pure number
  corresponding_lag <- lag_vals[max_index]
  
  # Ensure a single double is returned
  return(as.numeric(corresponding_lag))
}

findMaxCorr.fun <- function(x) {
  # If a grid cell contains NA, then set the result to NA
  if (is.na(mean(x))) {
    return(NA)
  } else {
    n <- length(x)
    
    # Get first (a) and second (b) variable
    a <- x[1:(n / 2)]
    b <- x[(n / 2 + 1):n]
    
    if (sd(a) == 0 || sd(b) == 0) return(NA_real_)
    
    # Convert to time series
    a <- ts(as.numeric(a))
    b <- ts(as.numeric(b))
    
    # Conduct cross-correlation
    ccf_result <- ccf(a, b, lag.max = 6, plot = FALSE)
    
    # Extract numeric vectors explicitly
    acf_vals <- as.numeric(ccf_result$acf)
    lag_vals <- as.numeric(ccf_result$lag)
    
    # Find the index of the maximum correlation
    max_index <- which.max(abs(ccf_result$acf))
    
    max_correlation <- ccf_result$acf[max_index]
    
    return(as.numeric(max_correlation))

  }
}

data <- c(gfed.anom.detrend, pr.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.pr <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.pr <- maxCorr

par(mfrow = c(1, 2))

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.pr,
  main = "Monthly time lags with largest correlation coefficient \n (Precipitation and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)

# Plot Correlation Coefficients
breaks <- seq(-1, 1, 0.1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  maxCorr.pr,
  main = "Corresponding Correlation coefficient \n (Precipitation and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)

# Repeat for TAS
data <- c(gfed.anom.detrend, tas.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.tas <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.tas <- maxCorr

par(mfrow = c(1, 2))

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.tas,
  main = "Monthly time lags with largest correlation coefficient \n (Temperature and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)

# Plot Correlation Coefficients
breaks <- seq(-1, 1, 0.1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  maxCorr.tas,
  main = "Corresponding Correlation coefficient \n (Temperature and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)

# Repeat for RSDS
data <- c(gfed.anom.detrend, rsds.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.rsds <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.rsds <- maxCorr

par(mfrow = c(1, 2))

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.rsds,
  main = "Monthly time lags with largest correlation coefficient \n (RSDS and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)

# Plot Correlation Coefficients
breaks <- seq(-1, 1, 0.1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  maxCorr.rsds,
  main = "Corresponding Correlation coefficient \n (RSDS and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)

# Repeat for Wind Speed
data <- c(gfed.anom.detrend, wind.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.wind <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.wind <- maxCorr

par(mfrow = c(1, 2))

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.wind,
  main = "Monthly time lags with largest correlation coefficient \n (Wind Speed and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)

# Plot Correlation Coefficients
breaks <- seq(-1, 1, 0.1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  maxCorr.wind,
  main = "Corresponding Correlation coefficient \n (Wind Speed and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)

# Repeat for VSM
data <- c(gfed.anom.detrend, vsm.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.vsm <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.vsm <- maxCorr

par(mfrow = c(1, 2))

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.vsm,
  main = "Monthly time lags with largest correlation coefficient \n (VSM and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)

# Plot Correlation Coefficients
breaks <- seq(-1, 1, 0.1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  maxCorr.vsm,
  main = "Corresponding Correlation coefficient \n (VSM and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)

# Repeat for LAI
data <- c(gfed.anom.detrend, lai.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.lai <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.lai <- maxCorr

par(mfrow = c(1, 2))

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.lai,
  main = "Monthly time lags with largest correlation coefficient \n (LAI and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)

# Plot Correlation Coefficients
breaks <- seq(-1, 1, 0.1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  maxCorr.lai,
  main = "Corresponding Correlation coefficient \n (LAI and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 1.0)


# Combining GC test and cross correlation
maxR2.pr <- maxCorr.pr^2 * p.value.pr
maxR2.tas <- maxCorr.tas^2 * p.value.tas
maxR2.rsds <- maxCorr.rsds^2 * p.value.rsds
maxR2.wind <- maxCorr.wind^2 * p.value.wind
maxR2.vsm <- maxCorr.vsm^2 * p.value.vsm
maxR2.lai <- maxCorr.lai^2 * p.value.lai


stacked_rasters <- c(maxR2.pr, maxR2.tas, maxR2.rsds, maxR2.wind, maxR2.vsm, maxR2.lai)
names(stacked_rasters) <- c("PR", "TAS", "RSDS", "Wind Speed", "VSM", "LAI")

# Apply a function to find the index of the max value for each cell
highest_index <- app(stacked_rasters, which.max)

# Set the levels to replace numeric values with the names of the rasters
levels(highest_index) <- data.frame(id = 1:6, name = names(stacked_rasters))

# Plot the result to visualize the highest value index for each cell
plot(highest_index, main = "Met variables that dominate GFED anomaly variability")
map("world2", add = TRUE, interior = FALSE)
