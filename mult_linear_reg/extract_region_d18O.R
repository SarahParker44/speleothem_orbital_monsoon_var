### Extract source area d18O of ppt for each region and time slice ##############

setwd("C:/Users/sarah/Documents/PhD/analyses/")

# load necessary packages

library(dplyr)
library(ncdf4)

# Load GISS simulated d18O
ncin <- nc_open("data/GISS_wiso.nc")
wiso <- ncvar_get(ncin, "H2O18_in_prec")
# Replace fill values with NA
fillvalue <- ncatt_get(ncin,"H2O18_in_prec","_FillValue")
wiso[wiso==fillvalue$value] <- NA
lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
nc_close(ncin); rm(ncin)

# Calculate anomalies to PI
wiso <- apply(wiso, c(1,2), function(x) x-wiso[,,1])
for (i in 1:dim(wiso)[3]){
  wiso[,,i] <- wiso[,,i] - wiso[,,1]
}
wiso <- wiso[,,-1]

# get land mask
ncin <- nc_open("data/GISS_landmask.nc")
mask <- ncvar_get(ncin, "frac_land")
nc_close(ncin)

mask <- ifelse(mask > 50, 1, NA) # 1  = land grids, NA = ocean grids

# Select land grids only
wiso_ls <- lapply(wiso_ls, function(x) x*mask)


# Define region limits
region_lims <- data.frame(region = c('ISM','EAM','SW-SAM1','SW-SAM2','IAM','NE-SAM','CAM','SAfM'),
                          min_lat = c(11,20,-10,-30,-24,-10,10,-30),
                          max_lat = c(32,39,0,-10,5,0,33,-17),
                          min_lon = c(50,100,-80,-68,95,-60,-115,10),
                          max_lon = c(95,125,-64,-40,135,-30,-58,40))

# create dataframe to fill with regional d18O average vals
region_d18O_GISS <- data.frame(region = rep(c('ISM','EAM','SW-SAM','IAM','NE-SAM','CAM','SAfM'),7),
                       t_slice = rep(c(1:6,9), each = 7),
                       region_d18O = NA)

#Subset to each region, and calc mean for each time slice
for (i in 1:length(wiso_ls)){ # each time slice
  dat <- wiso_ls[[i]] # filter to time slice
  for (k in 1:nrow(region_lims)){ # each region
    if (k %in% 3:4){ next } # don't calculate SAM in this section
    
    # Get region limits
    lon_vals <- region_lims[k,4:5]
    lat_vals <- region_lims[k,2:3]
    
    # Get row and col positions in matrix of region limits for subsetting
    lon_no <- sapply(lon_vals, function(x) which.min(abs(lon - x)))
    lat_no <- sapply(lat_vals, function(x) which.min(abs(lat - x)))
    
    # Subset to region, calculate mean
    dt <- mean(as.numeric(dat[lon_no[1]:lon_no[2], lat_no[2]:lat_no[1]]), na.rm = T)
    
    # save data in df
    time <- ifelse(i %in% 1:6, i, 9)
    reg_name <- as.character(region_lims$region[k])
    
    region_d18O_GISS[which(region_d18O_GISS$region == reg_name & region_d18O_GISS$t_slice == time),]$region_d18O <- dt
  }
}

# calc mean for SW-SAM

# region limits
lon_reg_SAM1 <- region_lims[3, 4:5]; lon_reg_SAM2 <- region_lims[4, 4:5]
lat_reg_SAM1 <- region_lims[3, 2:3]; lat_reg_SAM2 <- region_lims[4, 2:3]

# Get matrix col and row numbers for region
SAM1_lon <- sapply(lon_reg_SAM1, function(x) which.min(abs(lon-x))); SAM2_lon <- sapply(lon_reg_SAM2, function(x) which.min(abs(lon-x)))
SAM1_lat <- sapply(lat_reg_SAM1, function(x) which.min(abs(lat-x))); SAM2_lat <- sapply(lat_reg_SAM2, function(x) which.min(abs(lat-x)))

for (i in 1:length(wiso_ls)){ # each time slice
  dat <- wiso_ls[[i]] #filter to time slice
  df <- data.frame(expand.grid(long = lon, lati = lat), dat = as.numeric(dat))
  for (l in 1:nrow(df)){
    ln <- df$long[l]
    lt <- df$lati[l]
    if (ln %in% lon[SAM1_lon[1]:SAM1_lon[2]] & lt %in% lat[SAM1_lat[1]:SAM1_lat[2]] | 
        ln %in% lon[SAM2_lon[1]:SAM2_lon[2]] & lt %in% lat[SAM2_lat[1]:SAM2_lat[2]]){
      df[l,3] <- df[l,3] 
    } else {
      df[l,3] <- NA
    }
  }
  # mean over region
  dt <- mean(df$dat, na.rm = T)
  
  # save into df
  time <- ifelse(i %in% 1:6, i, 9)
  
  region_d18O_GISS[which(region_d18O_GISS$region == "SW-SAM" & region_d18O_GISS$t_slice == time),]$region_d18O <- dt
}


# Output data
write.csv(region_d18O_GISS, "data/multiple_regression/region_d18O.csv", row.names = F)

# additional output to compare with (+/- 4° around sites) composites
write.csv(filter(region_d18O_GISS, region %in% c("ISM","EAM","SW-SAM","IAM")), "C:/Users/sarah/Documents/PhD/analyses/rect_region_d18O.csv", row.names = F)
