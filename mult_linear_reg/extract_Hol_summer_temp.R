### Extract source summer temp ####

setwd("C:/Users/sarah/Documents/PhD/analyses/")

# load necessary packages

library(ncdf4)

# Load temp
# dim = 72 (lon) x 46 (lat) x 8 (time slice) x 120 (monthly mean for 10 decades)
ncin <- nc_open("C:/Users/sarah/Documents/PhD/analyses/data/GISS_clim_vars.nc")
temp <- ncvar_get(ncin, "tsurf")

# Replace fill values with NA
fillvalue <- ncatt_get(ncin,"tsurf","_FillValue")
temp[temp==fillvalue$value] <- NA
nc_close(ncin)

# Calculate monthly mean

monmean_vals <- array(NA, dim = c(72,46,8,12))
for (i in 1:12){ # each month
  dimvals <- numeric()
  for (j in 1:(dim(temp)[4]/12)){
    val <- ((j*12)-(12-i))
    dimvals[j] <- val # select i month in j decade
  }
  mon_temp <- temp[,,,dimvals] # filter to month
  monmean <- apply(mon_temp, c(1,2,3), mean, na.rm = T) # calc monmean
  monmean_vals[,,,i] <- monmean # save
}

# Get lon, lat and mask data
ncin <- nc_open("C:/Users/sarah/Documents/PhD/analyses/data/0ka_clim.nc")
lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
mask <- ncvar_get(ncin, "frac_land")
mask <- ifelse(mask < 50, 1, NA) # 1  = ocean grids, NA = land grids
nc_close(ncin)


# Calculate summer means

# MJJAS mean (NH summer)
MJJAS <- monmean_vals[,,,5:9]
MJJAS <- apply(MJJAS, c(1:3), mean) # mean across MJJAS

MJJAS_0ka <- MJJAS[,,1] # select PI

# MJJAS as anomaly to PI
for (i in 1:8){
  MJJAS[,,i] <- MJJAS[,,i] - MJJAS_0ka
}

MJJAS <- MJJAS[,,-1]

#apply mask to data (select ocean grids only)

for (i in 1:dim(MJJAS)[3]){
  MJJAS[,,i] <- MJJAS[,,i]*mask
}


# NDJFM (SH summer)

NDJFM <- monmean_vals[,,,c(11:12,1:3)]
NDJFM <- apply(NDJFM, c(1:3), mean)
NDJFM_0ka <- NDJFM[,,1]

# calculate anomalies to PI
for (i in 1:8){
  NDJFM[,,i] <- NDJFM[,,i] - NDJFM_0ka
}

NDJFM <- NDJFM[,,-1]

#apply mask to data
for (i in 1:dim(NDJFM)[3]){
  NDJFM[,,i] <- NDJFM[,,i]*mask
}


# Set (source) region limits
reg_lims <- data.frame(region = c("ISM","EAM","IAM","CAM","SW-SAM","NE-SAM","SAfM"),
                       min_lon = c(50,75,85,-80,-38,-35,20), max_lon = c(75,120,140,-40,-20,-13,55),
                       min_lat = c(0,0,-18,10,-20,-17,-35), max_lat = c(25,22,0,20,-10,-7,-10))

# create dataframe to put mean summer ppt values into
source_a <- data.frame(region = rep(c("ISM","EAM","IAM","CAM","SW-SAM","NE-SAM","SAfM"),7),
                       t_slice = rep(c(1:6,9), each = 7),
                       source_temp = NA) # NH summer vals

source_a2 <- data.frame(region = rep(c("ISM","EAM","IAM","CAM","SW-SAM","NE-SAM","SAfM"),7),
                       t_slice = rep(c(1:6,9), each = 7),
                       source_temp = NA) # SH summer vals

# Subset to each region, calculate mean for each time slice
for (i in 1:dim(MJJAS)[3]){ # each time slice
  dat_NH <- MJJAS[,,i] # subset to time slice
  dat_SH <- NDJFM[,,i]
  
  for (k in 1:nrow(reg_lims)){ # each region
    
    # Get region limits
    lon_vals <- reg_lims[k,2:3]
    lat_vals <- reg_lims[k,4:5]
    
    # Matrix col and row values of region limits
    lon_no <- sapply(lon_vals, function(x) which.min(abs(lon - x)))
    lat_no <- sapply(lat_vals, function(x) which.min(abs(lat - x)))
    
    # subset to region, calculate mean
    dt_NH <- mean(as.numeric(dat_NH[lon_no[1]:lon_no[2], lat_no[2]:lat_no[1]]), na.rm = T)
    dt_SH <- mean(as.numeric(dat_SH[lon_no[1]:lon_no[2], lat_no[2]:lat_no[1]]), na.rm = T)
    
    # save into df
    time <- ifelse(i %in% 1:6, i, 9)
    
    source_a[which(source_a$region == reg_lims$region[k] & source_a$t_slice == time),]$source_temp <- dt_NH
    source_a2[which(source_a2$region == reg_lims$region[k] & source_a2$t_slice == time),]$source_temp <- dt_SH
    
  }
}

# filter MJJAS to NH monsoon source regions, and NDJFM to SH regions
source_a2 <- source_a2 %>% filter(region %in% c("SW-SAM","NE-SAM","SAfM","IAM")) # SH
source_a <- source_a %>% filter(region %in% c("CAM","EAM","ISM")) # NH

# Combine NH and SH
GISS_all <- rbind(source_a, source_a2)

# Output
write.csv(GISS_all, "data/multiple_regression/source_summer_temp.csv", row.names = F)
