### Extract regional summer ppt ####

setwd("C:/Users/sarah/Documents/PhD/analyses/")

# load necessary packages

library(ncdf4)
library(data.table)
library(dplyr)

# Load ppt
# dim = 72 (lon) x 46 (lat) x 8 (time slice) x 120 (monthly mean for 10 decades)
ncin <- nc_open("data/GISS_clim_vars.nc")
pre <- ncvar_get(ncin, "ppt")

# Replace fill values with NA
fillvalue <- ncatt_get(ncin,"ppt","_FillValue")
pre[pre==fillvalue$value] <- NA
nc_close(ncin)

# Calculate monthly mean

monmean_vals <- array(NA, dim = c(72,46,8,12))
for (i in 1:12){ # each month
  dimvals <- numeric()
  for (j in 1:(dim(pre)[4]/12)){ # each decade
    val <- ((j*12)-(12-i)) 
    dimvals[j] <- val # select i month in j decade
  }
  dat <- pre[,,,dimvals] # filter to month
  monmean <- apply(dat, c(1,2,3), mean, na.rm = T) # calc monmean
  monmean_vals[,,,i] <- monmean # save
}


# Get lon, lat and mask data
ncin <- nc_open("data/0ka_clim.nc")
lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
mask <- ncvar_get(ncin, "frac_land")
mask <- ifelse(mask > 50, 1, NA) # 1  = land grids, NA = ocean grids
nc_close(ncin)


# Calculate summer means

# MJJAS (NH summer)
MJJAS <- monmean_vals[,,,5:9]
MJJAS <- apply(MJJAS, c(1:3), mean) # mean across MJJAS

MJJAS_0ka <- MJJAS[,,1] # select PI

# calculate anomalies to PI
for (i in 1:8){
  MJJAS[,,i] <- MJJAS[,,i] - MJJAS_0ka
}

MJJAS <- MJJAS[,,-1]

#apply mask to data (select land grids only)
for (i in 1:dim(MJJAS)[3]){ # each time slice
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


# Set region limits
region_lims <- data.frame(region = c('ISM','EAM','SW-SAM1','SW-SAM2','IAM','NE-SAM','CAM','SAfM'),
                          min_lat = c(11,20,-10,-30,-24,-10,10,-30),
                          max_lat = c(32,39,0,-10,5,0,33,-17),
                          min_lon = c(50,100,-80,-68,95,-60,-115,10),
                          max_lon = c(95,125,-64,-40,135,-30,-58,40))

# create dataframe to put mean summer ppt values into
region_ppt_GISS <- data.frame(region = rep(c('ISM','EAM','SW-SAM','IAM','NE-SAM','CAM','SAfM'),7),
                              t_slice = rep(c(1:6,9), each = 7),
                              region_ppt = NA) # NH summer values
region_ppt_GISS2 <- data.frame(region = rep(c('ISM','EAM','SW-SAM','IAM','NE-SAM','CAM','SAfM'),7),
                               t_slice = rep(c(1:6,9), each = 7),
                               region_ppt = NA) # SH summer values

# Subset to each region, calculate mean for each time slice
for (i in 1:dim(MJJAS)[3]){ # each time slice
  dat_NH <- MJJAS[,,i]
  dat_SH <- NDJFM[,,i]
  for (k in 1:nrow(region_lims)){ # each region
    if (k %in% 3:4){ next } # don't calculate SW-SAM in this section
    
    # get region limits
    lon_vals <- region_lims[k,4:5]
    lat_vals <- region_lims[k,2:3]
    
    # get matrix col and row vals for region limits
    lon_no <- sapply(lon_vals, function(x) which.min(abs(lon - x)))
    lat_no <- sapply(lat_vals, function(x) which.min(abs(lat - x)))
    
    # subset to region, calculate mean
    dt_NH <- mean(as.numeric(dat_NH[lon_no[1]:lon_no[2], lat_no[2]:lat_no[1]]), na.rm = T)
    dt_SH <- mean(as.numeric(dat_SH[lon_no[1]:lon_no[2], lat_no[2]:lat_no[1]]), na.rm = T)
    
    # save in df
    time <- ifelse(i %in% 1:6, i, 9)
    reg_name <- as.character(region_lims$region[k])
    
    region_ppt_GISS[which(region_ppt_GISS$region == reg_name & region_ppt_GISS$t_slice == time),]$region_ppt <- dt_NH
    region_ppt_GISS2[which(region_ppt_GISS2$region == reg_name & region_ppt_GISS2$t_slice == time),]$region_ppt <- dt_SH
  }
}


# calc mean for SW-SAM

# get region limits
lon_reg_SAM1 <- region_lims[3, 4:5]; lon_reg_SAM2 <- region_lims[4, 4:5]
lat_reg_SAM1 <- region_lims[3, 2:3]; lat_reg_SAM2 <- region_lims[4, 2:3]

# get matrix col and row vals for region limits
SAM1_lon <- sapply(lon_reg_SAM1, function(x) which.min(abs(lon-x))); SAM2_lon <- sapply(lon_reg_SAM2, function(x) which.min(abs(lon-x)))
SAM1_lat <- sapply(lat_reg_SAM1, function(x) which.min(abs(lat-x))); SAM2_lat <- sapply(lat_reg_SAM2, function(x) which.min(abs(lat-x)))

for (i in 1:dim(MJJAS)[3]){ # each time slice
  dat_SH <- NDJFM[,,i]
  df <- data.frame(expand.grid(long = lon, lati = lat), dat = as.numeric(dat_SH))
  for (l in 1:nrow(df)){ # subset to region
    ln <- df$long[l]
    lt <- df$lati[l]
    if (ln %in% lon[SAM1_lon[1]:SAM1_lon[2]] & lt %in% lat[SAM1_lat[1]:SAM1_lat[2]] | 
        ln %in% lon[SAM2_lon[1]:SAM2_lon[2]] & lt %in% lat[SAM2_lat[1]:SAM2_lat[2]]){
      df[l,3] <- df[l,3] 
    } else {
      df[l,3] <- NA
    }
  }
  dt <- mean(df$dat, na.rm = T)
  
  # save in df
  time <- ifelse(i %in% 1:6, i, 9)
  
  region_ppt_GISS2[which(region_ppt_GISS2$region == "SW-SAM" & region_ppt_GISS2$t_slice == time),]$region_ppt <- dt
}

# filter MJJAS to NH monsoon source regions, and NDJFM to SH regions
region_ppt_GISS2 <- region_ppt_GISS2 %>% filter(region %in% c("SW-SAM","NE-SAM","SAfM","IAM"))
region_ppt_GISS <- region_ppt_GISS %>% filter(region %in% c("CAM","EAM","ISM"))

# Combine NH and SH
GISS_all <- rbind(region_ppt_GISS, region_ppt_GISS2)

# output
write.csv(GISS_all, "data/multiple_regression/region_summer_ppt.csv", row.names = F)
