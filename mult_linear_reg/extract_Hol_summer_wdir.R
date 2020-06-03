### Extract summer wind dir anomalies ######

setwd("C:/Users/sarah/Documents/PhD/analyses/")

# load necessary packages

#library(dplyr)
library(ncdf4)

# Load surface wind direction
# dim = 72 (lon) x 46 (lat) x 8 (time slice) x 120 (monthly mean for 10 decades)
ncin <- nc_open("C:/Users/sarah/Documents/PhD/analyses/data/GISS_clim_vars.nc")
usurf <- ncvar_get(ncin, "usurf")
vsurf <- ncvar_get(ncin, "vsurf")

# Replace fill values with NA
fillvalue <- ncatt_get(ncin,"usurf","_FillValue")
usurf[usurf==fillvalue$value] <- NA
vsurf[vsurf==fillvalue$value] <- NA
nc_close(ncin); rm(ncin)

# Calculate monthly mean

monmean_usurf <- array(NA, dim = c(72,46,8,12))
monmean_vsurf <- array(NA, dim = c(72,46,8,12))
for (i in 1:12){ # each month
  dimvals <- numeric()
  for (j in 1:(dim(usurf)[4]/12)){
    val <- ((j*12)-(12-i))
    dimvals[j] <- val # select i month in j decade
  }
  mon_u <- usurf[,,,dimvals] # filter to month
  mon_v <- vsurf[,,,dimvals]
  
  monmean_u <- apply(mon_u, c(1,2,3), mean, na.rm = T) # calc monmean
  monmean_v <- apply(mon_v, c(1,2,3), mean, na.rm = T)
  
  monmean_usurf[,,,i] <- monmean_u # save
  monmean_vsurf[,,,i] <- monmean_v
}

# Get lon, lat and mask data
ncin <- nc_open("C:/Users/sarah/Documents/PhD/analyses/data/0ka_clim.nc")
lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
mask <- ncvar_get(ncin, "frac_land")
mask <- ifelse(mask < 50, 1, NA) # 1  = ocean grids, NA = land grids
nc_close(ncin)


# Calculate summer means

# MJJAS mean (NH summer)
MJJAS_u <- monmean_usurf[,,,5:9]; MJJAS_v <- monmean_vsurf[,,,5:9]
MJJAS_u <- apply(MJJAS_u, c(1:3), mean); MJJAS_v <- apply(MJJAS_v, c(1:3), mean)

# Calculate wind direction (as degrees from N)
MJJAS_winddir <- array(NA, dim = c(72,46,8))
for (i in 1:8){
  wind_abs = sqrt(MJJAS_u[,,i]^2 + MJJAS_v[,,i]^2) # calc abs windspeed: a^2 + b^2 = c^2 reorganised for c
  wind_dir_trig_to = atan2(MJJAS_u[,,i]/wind_abs, MJJAS_v[,,i]/wind_abs) # normalise u and v components using abs windspeed, conv m/s to radians
  wind_dir_trig_to_degrees = wind_dir_trig_to * 180/pi # radians to degrees
  wind_dir_PI <- ifelse(wind_dir_trig_to_degrees < 0, wind_dir_trig_to_degrees + 360, wind_dir_trig_to_degrees) # 0 to 360
  MJJAS_winddir[,,i] = wind_dir_PI
}

# Calculate anomalies
MJJAS_winddir_0ka <- MJJAS_winddir[,,1]; 
for (i in 1:8){
  MJJAS_winddir[,,i] <- MJJAS_winddir[,,i] - MJJAS_winddir_0ka
}

MJJAS_winddir <- MJJAS_winddir[,,-1]

# Change wind direction anomalies from -360->360 to -180->180, i.e. angle change should always be smaller one
MJJAS_winddir[which(MJJAS_winddir > 180)] <- MJJAS_winddir[which(MJJAS_winddir > 180)] - 360
MJJAS_winddir[which(MJJAS_winddir < -180)] <- 360 + MJJAS_winddir[which(MJJAS_winddir < -180)]

#apply mask to data
for (i in 1:dim(MJJAS_winddir)[3]){
  MJJAS_winddir[,,i] <- MJJAS_winddir[,,i]*mask
}


# NDJFM mean (SH summer)
NDJFM_u <- monmean_usurf[,,,c(11,12,1,2,3)]; NDJFM_v <- monmean_vsurf[,,,c(11,12,1,2,3)]
NDJFM_u <- apply(NDJFM_u, c(1:3), mean); NDJFM_v <- apply(NDJFM_v, c(1:3), mean)

# Calculate wind direction (as degrees from N)
NDJFM_winddir <- array(NA, dim = c(72,46,8))
for (i in 1:8){
  wind_abs = sqrt(NDJFM_u[,,i]^2 + NDJFM_v[,,i]^2)
  wind_dir_trig_to = atan2(NDJFM_u[,,i]/wind_abs, NDJFM_v[,,i]/wind_abs) 
  wind_dir_trig_to_degrees = wind_dir_trig_to * 180/pi 
  wind_dir_PI <- ifelse(wind_dir_trig_to_degrees < 0, wind_dir_trig_to_degrees + 360, wind_dir_trig_to_degrees)
  NDJFM_winddir[,,i] = wind_dir_PI
}

# Calculate anomalies
NDJFM_winddir_0ka <- NDJFM_winddir[,,1]
for (i in 1:8){
  NDJFM_winddir[,,i] <- NDJFM_winddir[,,i] - NDJFM_winddir_0ka
}

NDJFM_winddir <- NDJFM_winddir[,,-1]

# Change wind direction anomalies from -360->360 to -180->180, i.e. angle change should always be smaller one
NDJFM_winddir[which(NDJFM_winddir > 180)] <- NDJFM_winddir[which(NDJFM_winddir > 180)] - 360
NDJFM_winddir[which(NDJFM_winddir < -180)] <- 360 + NDJFM_winddir[which(NDJFM_winddir < -180)]


#apply mask to data
for (i in 1:dim(NDJFM_winddir)[3]){
  NDJFM_winddir[,,i] <- NDJFM_winddir[,,i]*mask
}


# region limits (source areas)
reg_lims <- data.frame(region = c("ISM","EAM","IAM","CAM","SW-SAM","NE-SAM","SAfM"),
                       min_lon = c(50,75,85,-80,-38,-35,20), max_lon = c(75,120,140,-40,-20,-13,55),
                       min_lat = c(0,0,-18,10,-20,-17,-35), max_lat = c(25,22,0,20,-10,-7,-10))

# create dataframe to fill with regional mean wind dir anomalies
source_a <- data.frame(region = rep(c("ISM","EAM","IAM","CAM","SW-SAM","NE-SAM","SAfM"),7),
                       t_slice = rep(c(1:6,9), each = 7),
                       wind_dir = NA) # NH summer

source_a2 <- data.frame(region = rep(c("ISM","EAM","IAM","CAM","SW-SAM","NE-SAM","SAfM"),7),
                        t_slice = rep(c(1:6,9), each = 7),
                        wind_dir = NA) # SH summer


# Subset to each region, calculate mean for each time slice
for (i in 1:dim(MJJAS_winddir)[3]){ # each time slice
  dat_NH <- MJJAS_winddir[,,i]
  dat_SH <- NDJFM_winddir[,,i]
  
  for (k in 1:nrow(reg_lims)){ # each region
    
    # get region limits
    lon_vals <- reg_lims[k,2:3]
    lat_vals <- reg_lims[k,4:5]
    
    # get col and row pos in matrix for region lims
    lon_no <- sapply(lon_vals, function(x) which.min(abs(lon - x)))
    lat_no <- sapply(lat_vals, function(x) which.min(abs(lat - x)))
    
    # subset to region, calculate mean
    dt_NH <- mean(as.numeric(dat_NH[lon_no[1]:lon_no[2], lat_no[2]:lat_no[1]]), na.rm = T)
    dt_SH <- mean(as.numeric(dat_SH[lon_no[1]:lon_no[2], lat_no[2]:lat_no[1]]), na.rm = T)
    
    # save in df
    time <- ifelse(i %in% 1:6, i, 9)
    
    source_a[which(source_a$region == reg_lims$region[k] & source_a$t_slice == time),]$wind_dir <- dt_NH
    source_a2[which(source_a2$region == reg_lims$region[k] & source_a2$t_slice == time),]$wind_dir <- dt_SH
    
  }
}

# filter MJJAS to NH monsoon source regions, and NDJFM to SH regions
source_a2 <- source_a2 %>% filter(region %in% c("SW-SAM","NE-SAM","SAfM","IAM"))
source_a <- source_a %>% filter(region %in% c("CAM","EAM","ISM"))

# Combine NH and SH
GISS_all <- rbind(source_a, source_a2)


# output
write.csv(GISS_all, "data/multiple_regression/source_summer_wdir.csv", row.names = F)
