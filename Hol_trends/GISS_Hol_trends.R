#### Produce regional composites of simulated d18Oprecip from GISS (with confidence intervals)

setwd("C:/Users/sarah/Documents/PhD/analyses/")

# load necessary packages

library(ncdf4)
library(data.table)
library(tidyr)
library(dplyr)

# Load GISS wiso data
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


# select land grids (>50%land)
ncin <- nc_open("data/GISS_landmask.nc")
mask <- ncvar_get(ncin, "frac_land")
nc_close(ncin)

mask <- ifelse(mask > 50, 1, NA) # 1  = land grids, NA = ocean grids

wiso_ls <- lapply(wiso_ls, function(x) x*mask)


# reshape into df
df <- data.frame(expand.grid(long = lon, lati = lat), 
                 `1` = as.numeric(wiso_ls[[1]]),
                 `2` = as.numeric(wiso_ls[[2]]),
                 `3` = as.numeric(wiso_ls[[3]]),
                 `4` = as.numeric(wiso_ls[[4]]),
                 `5` = as.numeric(wiso_ls[[5]]),
                 `6` = as.numeric(wiso_ls[[6]]),
                 `9` = as.numeric(wiso_ls[[7]]))

# 4) Load speleothem site data (from composites)

sites <- read.csv("data/composite_entities.csv")
sites <- unique(sites[,-c(3,4)]) # removes entity cols, select individual sites


# 5) Select areas around each site, apply distance weighting

# Function that selects grids around a site (grids only partially within selection given less weighting), and calculates linear distance weighted means
# x = site longitude, y = site latitude, mlength = +/- x degrees around site
wdist_means <- function(x, y, mlength){
  
  # select grids around site
  lon_range <- sapply(c(x-mlength,x+mlength), function(x) which.min(abs(lon-x)))
  lat_range <- sapply(c(y-mlength,y+mlength), function(x) which.min(abs(lat-x)))
  a <- lon_range[1]:lon_range[2]
  b <- lat_range[1]:lat_range[2]
  ab <- expand.grid(a,b); ab2 <- (ab$Var2-1)*nlon + ab$Var1 # get all grid row numbers for around site
  
  grids <- df[ab2,]
  
  # calculate distance weightings (linear distance cone: 1.01 * max distance within matrix / difference between grid midpoint and site)
  g_dist <- grids %>% mutate(dist = sqrt((long - x)^2 + (lati - y)^2)) # calculate Euclid distance between grid midpoint and site: a^2 + b^2 = c^2
  g_dist <- g_dist %>% mutate(d_weight = (max(dist)*1.01)-dist) # calculate dist wighting
  
  # I need to cutoff grids that are only partially within the +/- mlength grid, give less weighting
  
  grid_lims <- data.frame(lon = c(x, x+mlength, x, x-mlength, x+mlength, x+mlength, x-mlength, x-mlength),
                          lat = c(y+mlength, y, y-mlength, y, y+mlength, y-mlength, y-mlength, y+mlength)) # select grids at edge and corner of selection area
  
  j <- sapply(grid_lims$lon, function(x) which.min(abs(lon - x)))
  k <- sapply(grid_lims$lat, function(x) which.min(abs(lat - x))) # get lon and lat of grids from grid_lims
  
  lat_s <- diff(c(lat[2], lat[3])) # calc lat and long grid lengths
  lon_s <- diff(c(lon[2], lon[3]))
  
  grid_w <- data.frame(
    A = abs(diff(c(grid_lims$lat[1],lat[k[1]]-(lat_s/2))))/lat_s,
    B = abs(diff(c(grid_lims$lon[2],lon[j[2]]-(lon_s/2))))/lon_s,
    C = abs(diff(c(grid_lims$lat[3],lat[k[3]]+(lat_s/2))))/lat_s,
    D = abs(diff(c(grid_lims$lon[4],lon[j[4]]+(lon_s/2))))/lon_s) # calculate edge weightings (proportion of length of grid that lies within mlength of site)/length of total grid
  
  grid_w2 <- data.frame(
    E = (abs(diff(c(grid_lims$lat[1],lat[k[1]]-(lat_s/2)))) * abs(diff(c(grid_lims$lon[2],lon[j[2]]-(lon_s/2)))))/(lon_s*lat_s),
    G = (abs(diff(c(grid_lims$lon[2],lon[j[2]]-(lon_s/2)))) * abs(diff(c(grid_lims$lat[3],lat[k[3]]+(lat_s/2)))))/(lon_s*lat_s),
    H = (abs(diff(c(grid_lims$lat[3],lat[k[3]]+(lat_s/2)))) * abs(diff(c(grid_lims$lon[4],lon[j[4]]+(lon_s/2)))))/(lon_s*lat_s),
    I = (abs(diff(c(grid_lims$lon[4],lon[j[4]]+(lon_s/2)))) * abs(diff(c(grid_lims$lat[1],lat[k[1]]-(lat_s/2)))))/(lon_s*lat_s))  # calculate corner weightings (proportion of area of grid that lies within mlength of site)/ total area of grid
  
  grid_w2 <- data.frame(long = lon[j[5:8]], lati = lat[k[5:8]], edge_w = as.numeric(grid_w2)) #  df of corner grid weightings with lat and lon
  
  if(nrow(g_dist)-nrow(grid_w2) == 5){
    edge1 <- g_dist %>% filter(lati == max(lati) & long != min(long) & long != max(long)) %>% dplyr::select(c("long","lati")); edge1$edge_w <- grid_w$A # select northern edge grids lon and lat. apply edge weightings
    edge2 <- g_dist %>% filter(long == max(long) & lati != max(lati) & lati != min(lati)) %>% dplyr::select(c("long","lati")); edge2$edge_w <- grid_w$B # select eastern edge grids lon and lat
    edge3 <- g_dist %>% filter(lati == min(lati) & long != min(long) & long != max(long)) %>% dplyr::select(c("long","lati")); edge3$edge_w <- grid_w$C # select southern edge grids lon and lat
    edge4 <- g_dist %>% filter(long == min(long) & lati != max(lati) & lati != min(lati)) %>% dplyr::select(c("long","lati")); edge4$edge_w <- grid_w$D # select western edge grids lon and lat
    
    edge_w <- rbind(edge1, edge2, edge3, edge4, grid_w2) # combine all edge grid lon and lats
    
  } else if(nrow(g_dist)-nrow(grid_w2) == 2) {
    edge2 <- g_dist %>% filter(long == max(long) & lati != max(lati) & lati != min(lati)) %>% dplyr::select(c("long","lati")); edge2$edge_w <- grid_w$B # select eastern edge grids lon and lat
    edge4 <- g_dist %>% filter(long == min(long) & lati != max(lati) & lati != min(lati)) %>% dplyr::select(c("long","lati")); edge4$edge_w <- grid_w$D # select western edge grids lon and lat
    
    edge_w <- rbind(edge2, edge4, grid_w2)
  }
  
  # combine edge cut-offs with distance weightings
  g_dist <- left_join(g_dist, edge_w, by = c("long","lati"))
  g_dist$edge_w[which(is.na(g_dist$edge_w) == T)] <- 1 # non edge grids don't have contribution reduced
  
  g_dist$weight <- g_dist$d_weight * g_dist$edge_w # combined weighting of distance weighted cone, and edge cut-offs
  
  g_dist <- g_dist %>% filter(!is.na(X1)) # remove ocean grids (NA vals)
  
  
  # normalise weightings
  g_dist <- g_dist %>% mutate(d_weight_norm = weight / sum(weight)) 
  
  # Calculated weighted means
  g_dist2 <- g_dist %>% mutate_at(paste("X", c(1:6,9), sep = ""), list(var_w = ~.*d_weight_norm)) # multiple each d18O val at each time slice by weighting
  dist_avg <- g_dist2 %>% summarise_at(paste("X", c(1:6,9), "_var_w", sep = ""), list(var_avg = ~sum(.))) # sum across each time slice
  colnames(dist_avg) <- paste("ka_", c(1:6,9), sep = "")
  
  return(dist_avg)
}


# 6) Calculate weighted means for each site and time period

site_means <- list()
for (i in 1:nrow(sites)){ # each site
  site_means[[i]] <- data.frame(sites[i,], wdist_means(x = sites$longitude[i], y = sites$latitude[i], mlength = 4))
}
site_means <- rbindlist(site_means) # list to df


# 7) Bootstrapping

regions <- c("ISM","EAM","IAM","SW-SAM")
t_slices <- c(1:6,9)
nboot <- 1000
conf <- c(0.05,0.95)


# function that bootstrap resamples between time series extracted from each site
bootstrap <- function(DF){
  DF <- gather(DF, key = time_slice, value = var, 7:13) # reshape so that time slice is a column
  reg_list <- list()
  for (i in regions){ # each region
    reg_dat <- filter(DF, region == i) # select data for individual regions
    reg_list[[i]] <- spread(reg_dat[,-c(2:5)], key = site_id, value = var)[-1] # reshape into cols = sites, rows = time slices
  }
  
  locfit_list <- list()
  for (i in regions){
    result <- reg_list[[i]] # select regional data
    result <- result[,-1]
    
    mboot <- matrix(nrow = length(t_slices), ncol = nboot) # matrix to fill with bootstrapped data
    set.seed(1)
    
    for (j in 1:nboot){
      ne <- sample(seq(1,ncol(result),1), ncol(result), replace = T) # randomly sample sites
      y <- as.numeric(as.matrix(result)[,ne]) # select x and y data of randomly selected sites
      x <- as.numeric(rep(t_slices, length(ne)))
      
      mean_vals <- numeric()
      for (k in t_slices){
        val <- ifelse(k %in% 1:6, k, 7)
        mean_vals[val] <- mean(y[which(x == k)]) #mean trend of randomly sampled sites
      }
      mboot[, j] <- mean_vals
      rm(mean_vals)
    }
    
    bootci <- t(apply(mboot, 1, quantile, probs = conf, na.rm = T)) # get confidence intervals for mboot data
    
    rm(x, y)
    
    y <- c(as.matrix(result)) # get fit of all data
    x <- rep(t_slices, ncol(result))
    dat <- as.data.frame(cbind(x,y))
    
    avg <- dat %>% group_by(x) %>% summarise(fit = mean(y))
    
    
    locfit_list[[i]] <- as.data.frame(rbind(c(0,0,0,0), cbind(t_slices, avg$fit, bootci))) # save fit and bootstrapping confidence intervals 
    colnames(locfit_list[[i]]) <- c("age","locfit","conf5","conf95")
    
  }
  
  return(locfit_list)
}


region_wiso <- bootstrap(site_means) # bootstrap data
s <- rbindlist(region_wiso, idcol = T) # merge list into df
colnames(s)[1] <- "region"

# export data
write.csv(s, "data/Bootstrap_GISS_wiso.csv") # save


