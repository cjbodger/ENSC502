# Module 9 Replay

library(terra)
library(climetrics)
library(maps)
library(vars)

domain <- c(110, 160, -45, -10) # Set domain to Australia

# LAI
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_lai_2002-2023.nc")
data <- crop(x=data, y=domain)
data <- aggregate(x = data, fact=2, fun="mean") # More coarse spatial resolution
crs(data) <- "EPSG:4326"
lai <- data
mask <- mean(lai) # create mask
mask <- mask - mask + 1
rm(data)

# Precipitation
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_pr_2002-2023.nc")
data <- crop(x = data, y = domain)
data <- aggregate(x = data, fact=2, fun="mean") # More coarse spatial resolution
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
pr <- data
pr <- data * 1000 # convert unit from m to mm
rm(data)


print(lai)
print(data)
print(mask)
print(pr)
print(tas)
print(rsds)

# Near Surface Temp
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_tm_2002-2023.nc")
data <- crop(x = data, y = domain)
data <- aggregate(x = data, fact=2, fun="mean") # More coarse spatial resolution
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
tas <- data
rm(data)

# Surface Downwelling
data <- terra::rast("C:/Users/Charlotte/Desktop/RStudio/ERA5_ssrd_2002-2023.nc")
data <- crop(x = data, y = domain)
data <- aggregate(x = data, fact=2, fun="mean") # More coarse spatial resolution
crs(data) <- "EPSG:4326" # assign a map projection (WGS84)
data <- data * mask # exclude gridcells that both data sets do not have in common
rsds <- data / 86400 # Convert unit from J m-2 day-1 to W m-2 
rm(data)

data <- pr
# Create a sequence of dates
start_date <- as.Date("2001-01-15")
end_date <- as.Date("2023-12-15")
dates <- seq(from = start_date, to = end_date, by = "month")

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
trend.fun <- function(x) {
  # If a grid cell contains NA, then set the result to NA
  if (is.na(mean(x))) {
    return(NA)
  } else {
    time <- 1:length(x)
    linear.model <- lm(x ~ time)
    trend <- coef(linear.model)["time"]
    return(trend)
  }
}
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
pr.anom.detrend <- data.anom.detrend

my.col <- rev(map.pal("magma", n = 100))
plot(subset(data.anom.detrend, 1:1), col = my.col)
map("world2", add = TRUE)


# Repetition for TAS
data <- tas
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
tas.anom.detrend <- data.anom.detrend
my.col <- rev(map.pal("magma", n = 100))
plot(subset(data.anom.detrend, 1:1), col = my.col)
map("world2", add = TRUE)

# Repetition for DSR
data <- rsds
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
rsds.anom.detrend <- data.anom.detrend

min_val <- 0
max_val <- 2
my.col <- rev(map.pal("magma", n = 100))
plot(subset(data.anom.detrend, 1:1), col = my.col)
map("world2", add = TRUE, col = "white")


# Repetition for LAI
data <- lai
data.ts <- rts(data, dates)
data.12 <- apply.months(data.ts, 'mean')
data.clim <- rep(data.12, n)
data.anom <- data - data.clim
data.anom.detrend <- app(x = data.anom, fun = detrend.fun)
lai.anom.detrend <- data.anom.detrend
lai.det.2 <- subset(lai.anom.detrend, 1:1)
min_val <- -2.2304e-20
max_val <- 2.2304e-20 
my.col <- rev(map.pal("magma", n = 100))
plot(subset(lai.anom.detrend, 1), col = my.col)
map("world2", add = TRUE, col = "white")


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
  # Handle missing or tiny values
  if (is.na(mean(x)) || max(abs(x)) < 1e-10) {
    return(NA_real_)
  }
  
  n <- length(x)
  a <- ts(x[1:(n/2)])
  b <- ts(x[(n/2+1):n])
  tsDat <- ts.union(a, b)
  
  # Try VAR model safely
  tsVAR <- try(vars::VAR(tsDat, p = 2), silent = TRUE)
  if (inherits(tsVAR, "try-error")) return(NA_real_)
  
  # Try causality test safely
  caus <- try(vars::causality(tsVAR, cause = "a"), silent = TRUE)
  if (inherits(caus, "try-error")) return(NA_real_)
  
  # Ensure only one value is returned (first p-value)
  return(caus$Granger$p.value[1])
}


# PR GC LAI
data <- c(pr.anom.detrend, lai.anom.detrend)
p.value <- app(x = data, fun = granger.fun.2)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.pr <- p.value
plot(p.value, col = c("white", "orange"), main = "Precipitation Granger-causes LAI Anomalies")
map("world2", add = TRUE, interior = FALSE)

# TAS GC LAI
data <- c(tas.anom.detrend, lai.anom.detrend)
p.value <- app(x = data, fun = granger.fun.2)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.tas <- p.value
plot(p.value, col = c("white", "red"), main = "Temperature Granger-causes LAI Anomalies")
map("world2", add = TRUE, interior = FALSE)

# RSDS GC LAI
data <- c(rsds.anom.detrend, lai.anom.detrend)
p.value <- app(x = data, fun = granger.fun.2)

p.value[p.value >= 0.05] <- 0
p.value[p.value > 0] <- 1
p.value.rsds <- p.value
plot(p.value, col = c("white", "blue"), main = "Downwelling Granger-causes LAI Anomalies")
map("world2", add = TRUE, interior = FALSE)


# Summarizing
p.value <- p.value.pr + p.value.tas + p.value.rsds
plot(p.value, main = "PR, TAS, and RSDS Granger-causes LAI Anomalies")
map("world2", add = TRUE, interior = FALSE)

# Cross Correlation
lon <- 131
lat <- -17
coords <- matrix(c(lon, lat), ncol = 2, byrow = TRUE)
location <- vect(coords, type = "points")

plot(p.value)
map("world2", add = TRUE, interior = FALSE)
points(location, col = "red", pch = 16, cex = 1.0)

lai.anom.detrend.gc <- extract(lai.anom.detrend, location)
pr.anom.detrend.gc <- extract(pr.anom.detrend, location)
tas.anom.detrend.gc <- extract(tas.anom.detrend, location)
rsds.anom.detrend.gc <- extract(rsds.anom.detrend, location)

lai.anom.detrend.gc <- unlist(unname(as.vector(lai.anom.detrend.gc)))
pr.anom.detrend.gc <- unlist(unname(as.vector(pr.anom.detrend.gc)))
tas.anom.detrend.gc <- unlist(unname(as.vector(tas.anom.detrend.gc)))
rsds.anom.detrend.gc <- unlist(unname(as.vector(rsds.anom.detrend.gc)))

n <- length(lai.anom.detrend.gc)
lai.anom.detrend.gc <- lai.anom.detrend.gc[2:n]
pr.anom.detrend.gc <- pr.anom.detrend.gc[2:n]
tas.anom.detrend.gc <- tas.anom.detrend.gc[2:n]
rsds.anom.detrend.gc <- rsds.anom.detrend.gc[2:n]

lai.ts <- ts(lai.anom.detrend.gc)
pr.ts <- ts(pr.anom.detrend.gc)
tas.ts <- ts(tas.anom.detrend.gc)
rsds.ts <- ts(rsds.anom.detrend.gc)

tsDat <- ts.union(lai.ts, pr.ts, tas.ts, rsds.ts)

# Time Series Graph
plot(tsDat)

# Cross Correlation Comuptation
ccf_result <- ccf(lai.ts, pr.ts, lag.max = 12, plot = FALSE)
plot(ccf_result)
ccf_result <- ccf(lai.ts, tas.ts, lag.max = 12, plot = FALSE)
ccf_result <- ccf(lai.ts, rsds.ts, lag.max = 12, plot = FALSE)

# Index of maximum correlation
max_index <- which.max(ccf_result$acf)
max_correlation <- ccf_result$acf[max_index]
corresponding_lag <- ccf_result$lag[max_index]
result <- data.frame(Max_Correlation = max_correlation, Lag = corresponding_lag)
print(result)


# Finding time lag with largest correlation coefficient
findLag.fun <- function(x) {
  # If a grid cell contains NA, then set the result to NA
  if (is.na(mean(x))) {
    return(NA)
  } else {
    n <- length(x)
    
    # Get first (a) and second (b) variable
    a <- x[1:(n / 2)]
    b <- x[(n / 2 + 1):n]
    
    # Convert to time series
    a <- ts(a)
    b <- ts(b)
    
    # Conduct cross-correlation
    ccf_result <- ccf(a, b, lag.max = 12, plot = FALSE)
    
    # Find the index of the maximum correlation
    max_index <- which.max(abs(ccf_result$acf))
    
    # Extract the maximum correlation and corresponding lag
    max_correlation <- ccf_result$acf[max_index]
    corresponding_lag <- ccf_result$lag[max_index]
    
    # Combine results for display
    result <- corresponding_lag
    return(result)
  }
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
    
    # Convert to time series
    a <- ts(a)
    b <- ts(b)
    
    # Conduct cross-correlation
    ccf_result <- ccf(a, b, lag.max = 6, plot = FALSE)
    
    # Find the index of the maximum correlation
    max_index <- which.max(abs(ccf_result$acf))
    
    # Extract the maximum correlation and corresponding lag
    max_correlation <- ccf_result$acf[max_index]
    corresponding_lag <- ccf_result$lag[max_index]
    
    # Combine results for display
    result <- max_correlation
    return(result)
  }
}

data <- c(lai.anom.detrend, pr.anom.detrend)
