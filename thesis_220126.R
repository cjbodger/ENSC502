# Module 9 Replay
library(rmarkdown)
library(terra)
library(climetrics)
library(maps)
library(vars)
library(rnaturalearth)


# Area Burned (ESA FIRE CCI)
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/FB_2001_2019.nc")
crs(data) <- "EPSG:4326"
BA <- data
mask <- mean(BA) # create mask
mask <- mask - mask + 1
rm(data)
plot(mean(BA))


# Precipitation
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/FPR_2001_2022.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
pr <- data
pr <- data * 1000 # convert unit from m to mm
pr <- subset(pr, 1:228)
rm(data)

# Near Surface Temp
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/F2MT_2001_2022_K.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
tas <- data
tas <- subset(tas, 1:228)
rm(data)

# Relative Humidity
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_rh_2001_2022.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
rh <- data
rh <- subset(rh, 1:228)

# Wind Speed @ 10m
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/FWS_2001_2022.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
wind <- data
wind <- subset(wind, 1:228)
rm(data)

# Volumetric Soil Mositure Layer 1 (0-7cm)
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/FVSM_2001_2022.nc")
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
vsm <- data
vsm <- subset(vsm, 1:228)
rm(data)


data <- pr
# Create a sequence of dates
start_date <- as.Date("2001-01-01")
end_date <- as.Date("2019-12-01")
dates <- seq(from = start_date, to = end_date, by = "month")

data <- pr
# Raster Time Series
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts,'mean')
pr.12 <- data.12
my.col <- rev(map.pal("magma", n = 100))
my.title <- c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
plot(pr.12, col = my.col, main = my.title)
mtext("Mean Monthly PR (mm) from 2001 - 2019", side = 3, line = 3, outer = FALSE)

# Anomalies
n <- length(dates)/12
data.clim <- rep(pr.12, n)
data.anom <- data - data.clim
my.col <- rev(map.pal("magma", n = 100))
plot(subset(data.anom, 1:1), col = my.col)

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
mtext("Mean Monthly Surface Temp. (K) from 2001 - 2019", side = 3, line = 3, outer = FALSE)
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
tas.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(tas.anom.detrend, 1:1), col = my.col, main = "tas")

tas.sum <- sum(tas)
my.col <- rev(map.pal("inferno", n = 100))
plot(tas.sum, col = my.col, main = "Total Temperature from 2001 - 2019")

# Repetition for RH
data <- rh
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
rh.12 <- data.12
my.col <- rev(map.pal("magma", n = 100))
my.title <- c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
plot(rh.12, col = my.col, main = my.title)
mtext("Mean Monthly Relative Humidity (%) from 2001 - 2019", side = 3, line = 3, outer = FALSE)
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
rh.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(rh.anom.detrend, 1:1), col = my.col, main = "RH")


# Repetition for BA
data <- BA
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
BA.12 <- data.12
my.col <- rev(map.pal("magma", n = 100))
my.title <- c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
plot(BA.12, col = my.col, main = my.title)
mtext("Mean Monthly Total Area Burned from 2001 - 2019", side = 3, line = 3, outer = FALSE)

data <- BA
BA.sum <- sum(BA)
my.col <- rev(map.pal("inferno", n = 100))
plot(BA.sum, col = my.col, main = "Total Area Burned from 2001 - 2019")
map("world2", add = TRUE)

data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
BA.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("inferno", n = 40))
plot(subset(BA.anom.detrend, 1), col = my.col, main = "BA")


# Rep for wind
data <- wind
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
wind.12 <- data.12
my.col <- rev(map.pal("magma", n = 100))
my.title <- c("Jan", "Feb", "March", "April", "May", "June", "July", "Aug", "Sept", "Oct", "Nov", "Dec")
plot(wind.12, col = my.col, main = my.title)
mtext("Mean Monthly Wind Speed (m/s) from 2001 - 2019", side = 3, line = 3, outer = FALSE)
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
mtext("Mean Monthly Volumetric Soil Moisture (m^3/m^3) from 2001 - 2019", side = 3, line = 3, outer = FALSE)
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
vsm.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(vsm.anom.detrend, 1), col = my.col, main = "volumetric soil moisture layer 1")


granger.fun.3 <- function(x) {
  
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
    tsDat <- ts.union(a, b)
    tsVAR <- vars::VAR(tsDat, p = 12)
    # Apply Granger causality test 
    p.value <- c(vars::causality(tsVAR, cause = "a")$Granger[3]$p.value)
    p.value <- as.numeric(p.value)
  },
  error = function(e) {
    # Catch singular matrix, NA, or model fitting errors
    return(NA)
  })
  
  return(p.value)
}

# PR GC BA
data <- c(pr.anom.detrend, BA.anom.detrend)

# Messing with two different functions
p.value <- app(x = data, fun = granger.fun.3)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.pr <- p.value
plot(p.value.pr, col = c("white", "orange"), main = "Precipitation Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)

# TAS GC BA
data <- c(tas.anom.detrend, BA.anom.detrend)
p.value <- app(x = data, fun = granger.fun.3)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.tas <- p.value
plot(p.value, col = c("white", "red"), main = "Temperature Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)

# RH GC BA
data <- c(rh.anom.detrend, BA.anom.detrend)
p.value <- app(x = data, fun = granger.fun.3)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.rh <- p.value
plot(p.value, col = c("white", "cyan"), main = "Relative Humidity Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)

# Wind GC BA
data <- c(wind.anom.detrend, BA.anom.detrend)
p.value <- app(x = data, fun = granger.fun.3)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.wind <- p.value
plot(p.value, col = c("white", "purple"), main = "Wind Speed Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)

# VSM GC BA
data <- c(vsm.anom.detrend, BA.anom.detrend)
p.value <- app(x = data, fun = granger.fun.3)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.vsm <- p.value
plot(p.value, col = c("white", "magenta4"), main = "VSM Granger-causes Burned Area Anomalies")
map("world2", add = TRUE, interior = FALSE)


# Summarizing
p.value <- p.value.pr + p.value.tas + p.value.rh + p.value.wind + p.value.vsm
my.col <- c("#fff","#D8DE27", "#f89540", "#cc4778", "#7e03a8", "#0d0887")
plot(p.value, main = "PR, TAS, RH, WS, VSM Granger-causes Burned Area Anomalies", col = my.col)
map("world2", add = TRUE, interior = FALSE)

# Cross Correlation
lon <- 35
lat <- 0
coords <- matrix(c(lon, lat), ncol = 2, byrow = TRUE)
location <- vect(coords, type = "points")

plot(p.value, col = my.col)
map("world2", add = TRUE, interior = TRUE)
points(location, col = "black", pch = 1, cex = 2, lwd = 3)

my.col <- rev(map.pal("inferno", n = 100))
plot(BA.sum, col = my.col, main = "Total Area Burned from 2001 - 2019")
map("world2", add = TRUE)
points(location, col = "black", pch = 1, cex = 2, lwd = 3)

BA.anom.detrend.gc <- extract(BA.anom.detrend, location)
pr.anom.detrend.gc <- extract(pr.anom.detrend, location)
tas.anom.detrend.gc <- extract(tas.anom.detrend, location)
rh.anom.detrend.gc <- extract(rh.anom.detrend, location)
wind.anom.detrend.gc <- extract(wind.anom.detrend, location)
vsm.anom.detrend.gc <- extract(vsm.anom.detrend, location)

BA.anom.detrend.gc <- unlist(unname(as.vector(BA.anom.detrend.gc)))
pr.anom.detrend.gc <- unlist(unname(as.vector(pr.anom.detrend.gc)))
tas.anom.detrend.gc <- unlist(unname(as.vector(tas.anom.detrend.gc)))
rh.anom.detrend.gc <- unlist(unname(as.vector(rh.anom.detrend.gc)))
wind.anom.detrend.gc <- unlist(unname(as.vector(wind.anom.detrend.gc)))
vsm.anom.detrend.gc <- unlist(unname(as.vector(vsm.anom.detrend.gc)))

n <- length(BA.anom.detrend.gc)
BA.anom.detrend.gc <- BA.anom.detrend.gc[2:n]
pr.anom.detrend.gc <- pr.anom.detrend.gc[2:n]
tas.anom.detrend.gc <- tas.anom.detrend.gc[2:n]
rh.anom.detrend.gc <- rh.anom.detrend.gc[2:n]
wind.anom.detrend.gc <- wind.anom.detrend.gc[2:n]
vsm.anom.detrend.gc <- vsm.anom.detrend.gc[2:n]


BA.ts <- ts(BA.anom.detrend.gc)
pr.ts <- ts(pr.anom.detrend.gc)
tas.ts <- ts(tas.anom.detrend.gc)
rh.ts <- ts(rh.anom.detrend.gc)
wind.ts <- ts(wind.anom.detrend.gc)
vsm.ts <- ts(vsm.anom.detrend.gc)

tsDat <- ts.union(BA.ts, pr.ts, tas.ts, rh.ts, wind.ts, vsm.ts)

# Time Series Graph
plot(tsDat)


# Cross Correlation Computation
ccf_result <- ccf(BA.ts, pr.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
max_index <- which.max(ccf_result$acf)
max_correlation <- ccf_result$acf[max_index]
corresponding_lag <- ccf_result$lag[max_index]
result <- data.frame(Max_Correlation = max_correlation, Lag = corresponding_lag)
print ("Max Correlation Lag PR")
print(result)

ccf_result <- ccf(BA.ts, tas.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
max_index <- which.max(ccf_result$acf)
max_correlation <- ccf_result$acf[max_index]
corresponding_lag <- ccf_result$lag[max_index]
result <- data.frame(Max_Correlation = max_correlation, Lag = corresponding_lag)
print ("Max Correlation Lag TAS")
print(result)

ccf_result <- ccf(BA.ts, rh.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
max_index <- which.max(ccf_result$acf)
max_correlation <- ccf_result$acf[max_index]
corresponding_lag <- ccf_result$lag[max_index]
result <- data.frame(Max_Correlation = max_correlation, Lag = corresponding_lag)
print ("Max Correlation Lag RH")
print(result)

ccf_result <- ccf(BA.ts, wind.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
max_index <- which.max(ccf_result$acf)
max_correlation <- ccf_result$acf[max_index]
corresponding_lag <- ccf_result$lag[max_index]
result <- data.frame(Max_Correlation = max_correlation, Lag = corresponding_lag)
print ("Max Correlation Lag Wind")
print(result)

ccf_result <- ccf(BA.ts, vsm.ts, lag.max = 12, plot = FALSE)
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
    ccf_result <- ccf(a, b, lag.max = 12, plot = FALSE)
    
    # Extract numeric vectors explicitly
    acf_vals <- as.numeric(ccf_result$acf)
    lag_vals <- as.numeric(ccf_result$lag)
    
    # Find the index of the maximum correlation
    max_index <- which.max(abs(ccf_result$acf))
    
    max_correlation <- ccf_result$acf[max_index]
    
    return(as.numeric(max_correlation))

  }
}

data <- c(BA.anom.detrend, pr.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.pr <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.pr <- maxCorr

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.pr,
  main = "Monthly time lags with largest correlation coefficient \n (Precipitation and BA)",
  cex.main = 0.7,
  col = my.col,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Plot Correlation Coefficients
breaks <- seq(-1, 1, 0.1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  maxCorr.pr,
  main = "Corresponding Correlation coefficient \n (Precipitation and BA)",
  cex.main = 0.7,
  col = my.col,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Repeat for TAS
data <- c(tas.anom.detrend, BA.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.tas <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.tas <- maxCorr

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.tas,
  main = "Monthly time lags with largest correlation coefficient \n (Temperature and BA)",
  cex.main = 0.7,
  col = my.col,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Plot Correlation Coefficients
breaks <- seq(-1, 1, 0.1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  maxCorr.tas,
  main = "Corresponding Correlation coefficient \n (Temperature and BA)",
  cex.main = 0.7,
  col = my.col,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Repeat for RH
data <- c(BA.anom.detrend, rh.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.rh <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.rh <- maxCorr

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.rh,
  main = "Monthly time lags with largest correlation coefficient \n (RH and BA)",
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
  main = "Corresponding Correlation coefficient \n (RH and BA)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Repeat for Wind Speed
data <- c(BA.anom.detrend, wind.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.wind <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.wind <- maxCorr

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.wind,
  main = "Monthly time lags with largest correlation coefficient \n (Wind Speed and BA)",
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
  main = "Corresponding Correlation coefficient \n (Wind Speed and BA)",
  cex.main = 0.7,
  col = my.col,
  breaks = breaks,
  type = "continuous"
)
map("world2", add = TRUE, interior = FALSE)
points(location, pch = 1, cex = 2, lwd = 2)

# Repeat for VSM
data <- c(BA.anom.detrend, vsm.anom.detrend)
timeLags <- app(x = data, fun = findLag.fun)
timeLags.vsm <- timeLags
maxCorr <- app(x = data, fun = findMaxCorr.fun)
maxCorr.vsm <- maxCorr

# Plot time lags
breaks <- seq(-12, 12, 1)
my.col <- rev(map.pal("differences", n = length(breaks) - 1))

plot(
  timeLags.vsm,
  main = "Monthly time lags with largest correlation coefficient \n (VSM and BA)",
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
  main = "Corresponding Correlation coefficient \n (VSM and BA)",
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
plot(highest_index, main = "Met variables that dominate BA anomaly variability", col = my.col)
map("world2", add = TRUE, interior = FALSE)



# Summarizing by Frequency

freq_table <- freq(highest_index)

# Calculate percentages
freq_table$percentage <- (freq_table$count / sum(freq_table$count)) * 100

# Create pie chart
pie(freq_table$count, 
    labels = paste0(freq_table$value, "\n", round(freq_table$percentage, 1), "%"),
    col = my.col[1:5],
    main = "Global Dominant Variable for Burned Area Anomalies",
    cex = 0.8)


# Highest Index by Region Summary
tropics <- ext(0, 360, -23.5, 23.5)
highest_tropics <- crop(highest_index, tropics)
levels(highest_tropics) <- data.frame(id = 1:5, name = names(stacked_rasters))

freq_table <- freq(highest_tropics)

# Calculate percentages
freq_table$percentage <- (freq_table$count / sum(freq_table$count)) * 100

# Create pie chart
pie(freq_table$count, 
    labels = paste0(freq_table$value, "\n", round(freq_table$percentage, 1), "%"),
    col = my.col[1:5],
    main = "Global Dominant Variable for Burned Area Anomalies In Tropics",
    cex = 0.8)


temperate_north <- ext(0, 360, 23.5, 66.5)
highest_temp_north <- crop(highest_index, temperate_north)
levels(highest_temp_north) <- data.frame(id = 1:5, name = names(stacked_rasters))

freq_table <- freq(highest_temp_north)

# Calculate percentages
freq_table$percentage <- (freq_table$count / sum(freq_table$count)) * 100

# Create pie chart
pie(freq_table$count, 
    labels = paste0(freq_table$value, "\n", round(freq_table$percentage, 1), "%"),
    col = my.col[1:5],
    main = "Global Dominant Variable for Burned Area Anomalies in Temperate North",
    cex = 0.8)

polar_north <- ext(0, 360, 66.5, 90)
highest_polar_north <- crop(highest_index, polar_north)
levels(highest_polar_north) <- data.frame(id = 1:5, name = names(stacked_rasters))
plot(highest_polar_north)


freq_table <- freq(highest_polar_north)

# Calculate percentages
freq_table$percentage <- (freq_table$count / sum(freq_table$count)) * 100

# Create pie chart
pie(freq_table$count, 
    labels = paste0(freq_table$value, "\n", round(freq_table$percentage, 1), "%"),
    col = my.col[1:5],
    main = "Global Dominant Variable for Burned Area Anomalies in Polar North",
    cex = 0.8)

temperate_south <- ext(0, 360, -66.5, -23.5)
highest_temp_south <- crop(highest_index, temperate_south)
levels(highest_temp_south) <- data.frame(id = 1:5, name = names(stacked_rasters))


freq_table <- freq(highest_temp_south)

# Calculate percentages
freq_table$percentage <- (freq_table$count / sum(freq_table$count)) * 100

# Create pie chart
pie(freq_table$count, 
    labels = paste0(freq_table$value, "\n", round(freq_table$percentage, 1), "%"),
    col = my.col[1:5],
    main = "Global Dominant Variable for Burned Area Anomalies in Temperate South",
    cex = 0.8)

# Australia
aust <- ext(100, 170, -50, -10)
highest_aus <- crop(highest_index, aust)
levels(highest_aus) <- data.frame(id = 1:5, name = names(stacked_rasters))
plot(highest_aus)
freq_table <- freq(highest_aus)

# Calculate percentages
freq_table$percentage <- (freq_table$count / sum(freq_table$count)) * 100

# Create pie chart
pie(freq_table$count, 
    labels = paste0(freq_table$value, "\n", round(freq_table$percentage, 1), "%"),
    col = my.col[1:5],
    main = "Global Dominant Variable for Burned Area Anomalies in Australia",
    cex = 0.8)

# South America
souam <- ext(270, 340, -60, 15)
highest_souam <- crop(highest_index, souam)
levels(highest_souam) <- data.frame(id = 1:5, name = names(stacked_rasters))
plot(highest_souam)
freq_table <- freq(highest_souam)

# Calculate percentages
freq_table$percentage <- (freq_table$count / sum(freq_table$count)) * 100

# Create pie chart
pie(freq_table$count, 
    labels = paste0(freq_table$value, "\n", round(freq_table$percentage, 1), "%"),
    col = my.col[1:5],
    main = "Global Dominant Variable for Burned Area Anomalies in South America",
    cex = 0.8)

# North America
noram <- ext(190, 320, 15, 80)
highest_noram <- crop(highest_index, noram)
levels(highest_noram) <- data.frame(id = 1:5, name = names(stacked_rasters))
plot(highest_noram)
freq_table <- freq(highest_noram)

# Calculate percentages
freq_table$percentage <- (freq_table$count / sum(freq_table$count)) * 100

# Create pie chart
pie(freq_table$count, 
    labels = paste0(freq_table$value, "\n", round(freq_table$percentage, 1), "%"),
    col = my.col[1:5],
    main = "Global Dominant Variable for Burned Area Anomalies in North America",
    cex = 0.8)


# Time Series Analysis
BA_crop <- crop(BA.anom.detrend, tropics)
mean_BA <- global(BA_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_BA), main = "BA")
BA.ts <- (ts(mean_BA))

pr_crop <- crop(pr.anom.detrend, tropics)
mean_pr <- global(pr_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_pr), main = "PR")
pr.ts <- (ts(mean_pr))

rh_crop <- crop(rh.anom.detrend, tropics)
mean_rh <- global(rh_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_rh), main = "RH")
rh.ts <- (ts(mean_rh))

tas_crop <- crop(tas.anom.detrend, tropics)
mean_tas <- global(tas_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_tas), main = "TAS")
tas.ts <- (ts(mean_tas))

wind_crop <- crop(wind.anom.detrend, tropics)
mean_wind <- global(wind_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_wind), main = "WS")
wind.ts <- (ts(mean_wind))

vsm_crop <- crop(vsm.anom.detrend, tropics)
mean_vsm <- global(vsm_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_vsm), main = "VSM")
vsm.ts <- (ts(mean_vsm))


tsDat <- ts.union(BA.ts, pr.ts, tas.ts, rh.ts, wind.ts, vsm.ts)

plot(x = BA.ts, y = tas.ts)


# Time Series Graph
plot(tsDat, main = "Time Series (Mean) for Tropics")

# For Temp N
BA_crop <- crop(BA, temperate_north)
mean_BA <- global(BA_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_BA), main = "BA")
BA.ts <- (ts(mean_BA))

pr_crop <- crop(pr, temperate_north)
mean_pr <- global(pr_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_pr), main = "PR")
pr.ts <- (ts(mean_pr))

rh_crop <- crop(rh, temperate_north)
mean_rh <- global(rh_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_rh), main = "RH")
rh.ts <- (ts(mean_rh))

tas_crop <- crop(tas, temperate_north)
mean_tas <- global(tas_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_tas), main = "TAS")
tas.ts <- (ts(mean_tas))

wind_crop <- crop(wind, temperate_north)
mean_wind <- global(wind_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_wind), main = "WS")
wind.ts <- (ts(mean_wind))

vsm_crop <- crop(vsm, temperate_north)
mean_vsm <- global(vsm_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_vsm), main = "VSM")
vsm.ts <- (ts(mean_vsm))


tsDat <- ts.union(BA.ts, pr.ts, tas.ts, rh.ts, wind.ts, vsm.ts)

# Time Series Graph
plot(tsDat, main = "Time Series (Mean) for Temperate North")


# Polar North
BA_crop <- crop(BA, polar_north)
mean_BA <- global(BA_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_BA), main = "BA")
BA.ts <- (ts(mean_BA))

pr_crop <- crop(pr, polar_north)
mean_pr <- global(pr_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_pr), main = "PR")
pr.ts <- (ts(mean_pr))

rh_crop <- crop(rh, polar_north)
mean_rh <- global(rh_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_rh), main = "RH")
rh.ts <- (ts(mean_rh))

tas_crop <- crop(tas, polar_north)
mean_tas <- global(tas_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_tas), main = "TAS")
tas.ts <- (ts(mean_tas))

wind_crop <- crop(wind, polar_north)
mean_wind <- global(wind_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_wind), main = "WS")
wind.ts <- (ts(mean_wind))

vsm_crop <- crop(vsm, polar_north)
mean_vsm <- global(vsm_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_vsm), main = "VSM")
vsm.ts <- (ts(mean_vsm))


tsDat <- ts.union(BA.ts, pr.ts, tas.ts, rh.ts, wind.ts, vsm.ts)

# Time Series Graph
plot(tsDat, main = "Time Series (Mean) for Polar North")

# Temp South
BA_crop <- crop(BA, temperate_south)
mean_BA <- global(BA_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_BA), main = "BA")
BA.ts <- (ts(mean_BA))

pr_crop <- crop(pr, temperate_south)
mean_pr <- global(pr_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_pr), main = "PR")
pr.ts <- (ts(mean_pr))

rh_crop <- crop(rh, temperate_south)
mean_rh <- global(rh_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_rh), main = "RH")
rh.ts <- (ts(mean_rh))

tas_crop <- crop(tas, temperate_south)
mean_tas <- global(tas_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_tas), main = "TAS")
tas.ts <- (ts(mean_tas))

wind_crop <- crop(wind, temperate_south)
mean_wind <- global(wind_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_wind), main = "WS")
wind.ts <- (ts(mean_wind))

vsm_crop <- crop(vsm, temperate_south)
mean_vsm <- global(vsm_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_vsm), main = "VSM")
vsm.ts <- (ts(mean_vsm))


tsDat <- ts.union(BA.ts, pr.ts, tas.ts, rh.ts, wind.ts, vsm.ts)

# Time Series Graph
plot(tsDat, main = "Time Series (Mean) for Temperate South")


# Australia
BA_crop <- crop(BA, aust)
mean_BA <- global(BA_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_BA), main = "BA")
BA.ts <- (ts(mean_BA))

pr_crop <- crop(pr, aust)
mean_pr <- global(pr_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_pr), main = "PR")
pr.ts <- (ts(mean_pr))

rh_crop <- crop(rh, aust)
mean_rh <- global(rh_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_rh), main = "RH")
rh.ts <- (ts(mean_rh))

tas_crop <- crop(tas, aust)
mean_tas <- global(tas_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_tas), main = "TAS")
tas.ts <- (ts(mean_tas))

wind_crop <- crop(wind, aust)
mean_wind <- global(wind_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_wind), main = "WS")
wind.ts <- (ts(mean_wind))

vsm_crop <- crop(vsm, aust)
mean_vsm <- global(vsm_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_vsm), main = "VSM")
vsm.ts <- (ts(mean_vsm))


tsDat <- ts.union(BA.ts, pr.ts, tas.ts, rh.ts, wind.ts, vsm.ts)

# Time Series Graph
plot(tsDat, main = "Time Series (Mean) for Australia")


# South America
BA_crop <- crop(BA, souam)
mean_BA <- global(BA_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_BA), main = "BA")
BA.ts <- (ts(mean_BA))

pr_crop <- crop(pr, souam)
mean_pr <- global(pr_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_pr), main = "PR")
pr.ts <- (ts(mean_pr))

rh_crop <- crop(rh, souam)
mean_rh <- global(rh_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_rh), main = "RH")
rh.ts <- (ts(mean_rh))

tas_crop <- crop(tas, souam)
mean_tas <- global(tas_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_tas), main = "TAS")
tas.ts <- (ts(mean_tas))

wind_crop <- crop(wind, souam)
mean_wind <- global(wind_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_wind), main = "WS")
wind.ts <- (ts(mean_wind))

vsm_crop <- crop(vsm, souam)
mean_vsm <- global(vsm_crop, fun = "mean", na.rm = TRUE)
plot(ts(mean_vsm), main = "VSM")
vsm.ts <- (ts(mean_vsm))


tsDat <- ts.union(BA.ts, pr.ts, tas.ts, rh.ts, wind.ts, vsm.ts)

# Time Series Graph
plot(tsDat, main = "Time Series (Mean) for South America")