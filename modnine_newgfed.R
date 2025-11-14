# Module 9 Replay
library(rmarkdown)
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

# Relative Humidity
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_rh.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
rh <- data

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


data <- pr
# Create a sequence of dates
start_date <- as.Date("2001-01-01")
end_date <- as.Date("2015-12-01")
dates <- seq(from = start_date, to = end_date, by = "month")

data <- pr
# Raster Time Series
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts,'mean')
pr.12 <- data.12
my.col <- rev(map.pal("magma", n = 100))
my.title <- c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
plot(pr.12, col = my.col, main = my.title)
mtext("Mean Monthly PR (mm) from 2001 - 2015", side = 3, line = 1, outer = FALSE)

# Anomalies
n <- length(dates)/12
data.clim <- rep(pr.12, n)
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
tas.12 <- data.12
my.col <- rev(map.pal("magma", n = 100))
my.title <- c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
plot(tas.12, col = my.col, main = my.title)
mtext("Mean Monthly Surface Temp. (K) from 2001 - 2015", side = 3, line = 1, outer = FALSE)
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
tas.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(tas.anom.detrend, 1:1), col = my.col, main = "tas")

# Repetition for RH
data <- rh
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
rh.12 <- data.12
my.col <- rev(map.pal("magma", n = 100))
my.title <- c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
plot(rh.12, col = my.col, main = my.title)
mtext("Mean Monthly Relative Humidity (%) from 2001 - 2015", side = 3, line = 1, outer = FALSE)
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
rh.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(rh.anom.detrend, 1:1), col = my.col, main = "RH")


# Repetition for GFED
data <- gfed
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
gfed.12 <- data.12
my.col <- rev(map.pal("magma", n = 100))
my.title <- c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
plot(gfed.12, col = my.col, main = my.title)
mtext("Mean Monthly Frational Area Burned from 2001 - 2015", side = 3, line = 1, outer = FALSE)

data <- gfed
gfed.sum <- sum(gfed)
my.col <- rev(map.pal("inferno", n = 100))
plot(gfed.sum, col = my.col, main = "Total Frational Area Burned from 2001 - 2015")
map("world2", add = TRUE)

data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
gfed.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("viridis", n = 100))
plot(subset(gfed.anom.detrend, 1), col = my.col, main = "gfed")


# Rep for wind
data <- wind
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
wind.12 <- data.12
my.col <- rev(map.pal("magma", n = 100))
my.title <- c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
plot(wind.12, col = my.col, main = my.title)
mtext("Mean Monthly Wind Speed (m/s) from 2001 - 2015", side = 3, line = 1, outer = FALSE)
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
vsm.12 <- data.12
my.col <- rev(map.pal("magma", n = 100))
my.title <- c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
plot(vsm.12, col = my.col, main = my.title)
mtext("Mean Monthly Volumetric Soil Moisture (m^3/m^3) from 2001 - 2015", side = 3, line = 1, outer = FALSE)
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
vsm.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(vsm.anom.detrend, 1), col = my.col, main = "volumetric soil moisture layer 1")

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

# RH GC GFED
data <- c(rh.anom.detrend, gfed.anom.detrend)
p.value <- app(x = data, fun = granger.fun.2)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.rh <- p.value
plot(p.value, col = c("white", "cyan"), main = "Relative Humidity Granger-causes Burned Area Anomalies")
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
plot(p.value, col = c("white", "magenta4"), main = "VSM Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)


# Summarizing
p.value <- p.value.pr + p.value.tas + p.value.rh + p.value.wind + p.value.vsm
my.col <- c("#fff","#ff595e", "#ffca3a", "#8ac926", "#1982c4", "#6a4c93")
plot(p.value, main = "PR, TAS, RH, WS, VSM Granger-causes Burned Area Anomalies", col = my.col)
map("world2", add = TRUE, interior = FALSE)

# Cross Correlation (Switching Location to Central Africa)
lon <- 35
lat <- 0
coords <- matrix(c(lon, lat), ncol = 2, byrow = TRUE)
location <- vect(coords, type = "points")

plot(p.value, col = my.col)
map("world2", add = TRUE, interior = FALSE)
points(location, col = "black", pch = 4, cex = 2, lwd = 2)

gfed.anom.detrend.gc <- extract(gfed.anom.detrend, location)
pr.anom.detrend.gc <- extract(pr.anom.detrend, location)
tas.anom.detrend.gc <- extract(tas.anom.detrend, location)
rh.anom.detrend.gc <- extract(rh.anom.detrend, location)
wind.anom.detrend.gc <- extract(wind.anom.detrend, location)
vsm.anom.detrend.gc <- extract(vsm.anom.detrend, location)

gfed.anom.detrend.gc <- unlist(unname(as.vector(gfed.anom.detrend.gc)))
pr.anom.detrend.gc <- unlist(unname(as.vector(pr.anom.detrend.gc)))
tas.anom.detrend.gc <- unlist(unname(as.vector(tas.anom.detrend.gc)))
rh.anom.detrend.gc <- unlist(unname(as.vector(rh.anom.detrend.gc)))
wind.anom.detrend.gc <- unlist(unname(as.vector(wind.anom.detrend.gc)))
vsm.anom.detrend.gc <- unlist(unname(as.vector(vsm.anom.detrend.gc)))


n <- length(gfed.anom.detrend.gc)
gfed.anom.detrend.gc <- gfed.anom.detrend.gc[2:n]
pr.anom.detrend.gc <- pr.anom.detrend.gc[2:n]
tas.anom.detrend.gc <- tas.anom.detrend.gc[2:n]
rh.anom.detrend.gc <- rh.anom.detrend.gc[2:n]
wind.anom.detrend.gc <- wind.anom.detrend.gc[2:n]
vsm.anom.detrend.gc <- vsm.anom.detrend.gc[2:n]


gfed.ts <- ts(gfed.anom.detrend.gc)
pr.ts <- ts(pr.anom.detrend.gc)
tas.ts <- ts(tas.anom.detrend.gc)
rh.ts <- ts(rh.anom.detrend.gc)
wind.ts <- ts(wind.anom.detrend.gc)
vsm.ts <- ts(vsm.anom.detrend.gc)

tsDat <- ts.union(gfed.ts, pr.ts, tas.ts, rh.ts, wind.ts, vsm.ts)

# Time Series Graph
plot(tsDat)

# Cross Correlation Computation
ccf_result <- ccf(gfed.ts, pr.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
max_index <- which.max(ccf_result$acf)
max_correlation <- ccf_result$acf[max_index]
corresponding_lag <- ccf_result$lag[max_index]
result <- data.frame(Max_Correlation = max_correlation, Lag = corresponding_lag)
print ("Max Correlation Lag PR")
print(result)

ccf_result <- ccf(gfed.ts, tas.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
max_index <- which.max(ccf_result$acf)
max_correlation <- ccf_result$acf[max_index]
corresponding_lag <- ccf_result$lag[max_index]
result <- data.frame(Max_Correlation = max_correlation, Lag = corresponding_lag)
print ("Max Correlation Lag TAS")
print(result)

ccf_result <- ccf(gfed.ts, rh.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
max_index <- which.max(ccf_result$acf)
max_correlation <- ccf_result$acf[max_index]
corresponding_lag <- ccf_result$lag[max_index]
result <- data.frame(Max_Correlation = max_correlation, Lag = corresponding_lag)
print ("Max Correlation Lag RH")
print(result)

ccf_result <- ccf(gfed.ts, wind.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
max_index <- which.max(ccf_result$acf)
max_correlation <- ccf_result$acf[max_index]
corresponding_lag <- ccf_result$lag[max_index]
result <- data.frame(Max_Correlation = max_correlation, Lag = corresponding_lag)
print ("Max Correlation Lag Wind")
print(result)

ccf_result <- ccf(gfed.ts, vsm.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
max_index <- which.max(ccf_result$acf)
max_correlation <- ccf_result$acf[max_index]
corresponding_lag <- ccf_result$lag[max_index]
result <- data.frame(Max_Correlation = max_correlation, Lag = corresponding_lag)
print ("Max Correlation Lag VSM")
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

# Plot time lags
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.pr,
  main = "Monthly time lags with largest correlation coefficient \n (Precipitation and GFED)",
  cex.main = 0.7,
  col = my.col,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Plot Correlation Coefficients
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  maxCorr.pr,
  main = "Corresponding Correlation coefficient \n (Precipitation and GFED)",
  cex.main = 0.7,
  col = my.col,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Repeat for TAS
data <- c(gfed.anom.detrend, tas.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.tas <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.tas <- maxCorr

# Plot time lags
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.tas,
  main = "Monthly time lags with largest correlation coefficient \n (Temperature and GFED)",
  cex.main = 0.7,
  col = my.col,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Plot Correlation Coefficients
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  maxCorr.tas,
  main = "Corresponding Correlation coefficient \n (Temperature and GFED)",
  cex.main = 0.7,
  col = my.col,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Repeat for RH
data <- c(gfed.anom.detrend, rh.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.rh <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.rh <- maxCorr

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.rh,
  main = "Monthly time lags with largest correlation coefficient \n (RH and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Plot Correlation Coefficients
breaks <- seq(-1, 1, 0.1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  maxCorr.rh,
  main = "Corresponding Correlation coefficient \n (RH and GFED)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Repeat for Wind Speed
data <- c(gfed.anom.detrend, wind.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.wind <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.wind <- maxCorr

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
points(location, pch = 1, cex = 2, lwd = 2)

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
points(location, pch = 1, cex = 2, lwd = 2)

# Repeat for VSM
data <- c(gfed.anom.detrend, vsm.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.vsm <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.vsm <- maxCorr

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
points(location, pch = 1, cex = 2, lwd = 2)

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
points(location, pch = 1, cex = 2, lwd = 2)

# Combining GC test and cross correlation
maxR2.pr <- maxCorr.pr^2 * p.value.pr
maxR2.tas <- maxCorr.tas^2 * p.value.tas
maxR2.rh <- maxCorr.rh^2 * p.value.rh
maxR2.wind <- maxCorr.wind^2 * p.value.wind
maxR2.vsm <- maxCorr.vsm^2 * p.value.vsm

# Set all zeros to NA
maxR2.pr[maxR2.pr==0] <- NA 
maxR2.tas[maxR2.tas==0] <- NA 
maxR2.rh[maxR2.rh==0] <- NA 
maxR2.wind[maxR2.wind==0] <- NA 
maxR2.vsm[maxR2.vsm==0] <- NA 

stacked_rasters <- c(maxR2.pr, maxR2.tas, maxR2.rh, maxR2.wind, maxR2.vsm)
names(stacked_rasters) <- c("PR", "TAS", "RH", "Wind Speed", "VSM")

# Apply a function to find the index of the max value for each cell
highest_index <- app(stacked_rasters, which.max, na.rm = TRUE)

# Set the levels to replace numeric values with the names of the rasters
levels(highest_index) <- data.frame(id = 1:5, name = names(stacked_rasters))

# Plot the result to visualize the highest value index for each cell
my.col <- c("#de324c","#f4895f", "#f8e16f", "#95cf92", "#369acc", "#9656a2")
plot(highest_index, main = "Met variables that dominate GFED anomaly variability", col = my.col)
map("world2", add = TRUE, interior = FALSE)

