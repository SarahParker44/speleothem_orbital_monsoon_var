library(R.matlab)
library(dplyr)
library(ncdf4)
library(tidyr)
library(data.table)

setwd("C:/Users/sarah/Documents/PhD/analyses/data/")

# Load Speleothem d18O data
spel_d18O <- read.csv("speleothem_boxplot_d18O.csv")

# Extract ECHAM d18O, wiso and temp data around sites used in speleothem box plots.
# STEP 1: load data (from ECHAM_dat.mat, generated in "Step3_calc_anom.m")
# STEP 2: Extract model data around each site and calculate distance weighted means
# STEP 3: Save data for plotting
# Steps carried out for a) wiso, b) ppt and c) temp data

degs_site <- 3

# STEP 1: load data and calculate anomalies

# Load Echam wiso data

isoECH <- readMat("C:/Users/sarah/Documents/PhD/analyses/data/ECHAM_dat.mat")

LIG_anom <- isoECH$LIG.wiso.anom
LGM_anom <- isoECH$LGM.wiso.anom
MH_anom <- isoECH$MH.wiso.anom

# get lon and lat
ncin <- nc_open("C:/Users/sarah/Documents/PhD/analyses/data/Echam5-wiso.T106_EXP010_PL_PI.global.monmean.d18O_prec.nc")
lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
lon1 <- lon[162:length(lon)]; lon2 <- lon[1:161]; lon <- c(lon1 - 360, lon2)
nlon <- length(lon); nlat <- length(lat)
nc_close(ncin); rm(ncin)

# use land grids only, ocean grids -> NA
LIG_mask <- isoECH$land.mask.LIG

slm <- ncvar_get(nc_open("T106.slm.w_miss.nc"))
a <- slm[161:320,]; b <- slm[1:160,]; slm <- rbind(a,b)

LIG_anom <- LIG_anom*LIG_mask
LGM_anom <- LGM_anom*slm
MH_anom <- MH_anom*slm



# STEP 2: Extract model data for speleothem sites

# Load sites
sites <- read.csv("boxplot_sites.csv")
sites <- sites[,-c(5,7)]; sites <- unique(sites)


# reshape ECHAM data into dataframes where cols = lon, lat, wiso anom and rows = individual grids
isoMH_df <- data.frame(expand.grid(lon = lon, lat = lat), wiso = as.numeric(MH_anom))
isoLGM_df <- data.frame(expand.grid(lon = lon, lat = lat), wiso = as.numeric(LGM_anom))
isoLIG_df <- data.frame(expand.grid(lon = lon, lat = lat), wiso = as.numeric(LIG_anom))


# Function for extracting data around an individual site, and calculating distance weighted mean
# Inputs: x = site lon, y = site lat, mlength = distance around each site for extracting data, dt = model data, as dataframe
Euclid_dist <- function(x, y, mlength, dt){
  
  # select grids around site
  lon_range <- sapply(c(x-mlength,x+mlength), function(x) which.min(abs(lon-x)))
  lat_range <- sapply(c(y-mlength,y+mlength), function(x) which.min(abs(lat-x)))
  a <- lon_range[1]:lon_range[2]
  b <- lat_range[1]:lat_range[2]
  ab <- expand.grid(a,b); ab2 <- (ab$Var2-1)*nlon + ab$Var1 # get all grid row numbers for around site
  
  grids <- dt[ab2,]
  
  
  # calculate distance weightings
  g_dist <- grids %>% mutate(dist = sqrt((lon - x)^2 + (lat - y)^2))
  g_dist <- g_dist %>% mutate(d_weight = (max(dist)*1.01)-dist)
  
  # Cutoff grids that are only partially within the +/- slength grid, weight according to proportion of grid within area
  
  grid_lims <- data.frame(lon = c(x, x+mlength, x, x-mlength, x+mlength, x+mlength, x-mlength, x-mlength),
                          lat = c(y+mlength, y, y-mlength, y, y+mlength, y-mlength, y-mlength, y+mlength)) # select the edge and corner grids
  
  j <- sapply(grid_lims$lon, function(x) which.min(abs(lon - x)))
  k <- sapply(grid_lims$lat, function(x) which.min(abs(lat - x)))
  
  lat_s <- diff(c(lat[2], lat[1]))
  lon_s <- diff(c(lon[1], lon[2]))
  
  grid_w <- data.frame(
    A = abs(diff(c(grid_lims$lat[1],lat[k[1]]-(lat_s/2))))/lat_s,
    B = abs(diff(c(grid_lims$lon[2],lon[j[2]]-(lon_s/2))))/lon_s,
    C = abs(diff(c(grid_lims$lat[3],lat[k[3]]+(lat_s/2))))/lat_s,
    D = abs(diff(c(grid_lims$lon[4],lon[j[4]]+(lon_s/2))))/lon_s) # calculate edge grid weightings
  
  grid_w2 <- data.frame(
    E = (abs(diff(c(grid_lims$lat[1],lat[k[1]]-(lat_s/2)))) * abs(diff(c(grid_lims$lon[2],lon[j[2]]-(lon_s/2)))))/(lon_s*lat_s),
    G = (abs(diff(c(grid_lims$lon[2],lon[j[2]]-(lon_s/2)))) * abs(diff(c(grid_lims$lat[3],lat[k[3]]+(lat_s/2)))))/(lon_s*lat_s),
    H = (abs(diff(c(grid_lims$lat[3],lat[k[3]]+(lat_s/2)))) * abs(diff(c(grid_lims$lon[4],lon[j[4]]+(lon_s/2)))))/(lon_s*lat_s),
    I = (abs(diff(c(grid_lims$lon[4],lon[j[4]]+(lon_s/2)))) * abs(diff(c(grid_lims$lat[1],lat[k[1]]-(lat_s/2)))))/(lon_s*lat_s))
  
  grid_w2 <- data.frame(lon = lon[j[5:8]], lat = lat[k[5:8]], edge_w = as.numeric(grid_w2))
  
  edge1 <- g_dist %>% filter(lat == max(lat) & lon != min(lon) & lon != max(lon)) %>% dplyr::select(c("lon","lat")); edge1$edge_w <- grid_w$A
  edge2 <- g_dist %>% filter(lon == max(lon) & lat != max(lat) & lat != max(lat)) %>% dplyr::select(c("lon","lat")); edge2$edge_w <- grid_w$B
  edge3 <- g_dist %>% filter(lat == min(lat) & lon != min(lon) & lon != max(lon)) %>% dplyr::select(c("lon","lat")); edge3$edge_w <- grid_w$C
  edge4 <- g_dist %>% filter(lon == min(lon) & lat != max(lat) & lat != min(lat)) %>% dplyr::select(c("lon","lat")); edge4$edge_w <- grid_w$D
  
  edge_w <- rbind(edge1, edge2, edge3, edge4, grid_w2)
  
  # combine edge cut-offs with distance weightings
  g_dist <- left_join(g_dist, edge_w, by = c("lon","lat"))
  g_dist$edge_w[which(is.na(g_dist$edge_w) == T)] <- 1 # non edge grids don't have contribution reduced
  g_dist$weight <- g_dist$d_weight * g_dist$edge_w
  
  var_name <- colnames(g_dist)[3]
  
  g_dist <- na.omit(g_dist, cols = var_name) # remove ocean grid
  
  # normalise
  g_dist <- g_dist %>% mutate(d_weight_norm = weight / sum(weight)) 
  
  # Calculated weighted means
  g_dist$w_var <- g_dist[,paste(var_name)]*g_dist$d_weight_norm
  dist_avg <- g_dist %>% summarise(w_avg = sum(w_var))
  
  return(dist_avg)
}

# Calculate distance weighted mean for each site in each time period
#ECHAM
site_means_E <- list()
for (i in 1:nrow(sites)){
  site_means_E[[i]] <- data.frame(sites[i,], 
                                  Euclid_dist(x = sites$lon[i], y = sites$lat[i], mlength = degs_site, dt = isoMH_df),
                                  Euclid_dist(x = sites$lon[i], y = sites$lat[i], mlength = degs_site, dt = isoLGM_df),
                                  Euclid_dist(x = sites$lon[i], y = sites$lat[i], mlength = degs_site, dt = isoLIG_df))
}
site_means_E <- rbindlist(site_means_E) # Combine list into dataframe
colnames(site_means_E)[6:8] <- c("MH","LGM","LIG") # set column names

wiso_Echam <- site_means_E %>% gather(key = time_slice, value = wiso_anom, 6:8) # reshape so that time slices are a variable in df



# SAME STEPS AGAIN FOR PRE DATA
rm(list= ls()[!(ls() %in% c('wiso_Echam','spel_d18O','degs_site','Euclid_dist'))])

# STEP 1: load data and calculate anomalies

# Load (resampled) Echam ppt data
ECH <- readMat("ECHAM_dat.mat")

LIG_anom <- ECH$LIG.pre.anom
LGM_anom <- ECH$LGM.pre.anom
MH_anom <- ECH$MH.pre.anom

LIG_mask <- ECH$land.mask.LIG

slm <- ncvar_get(nc_open("T106.slm.w_miss.nc"), "slm")
a <- slm[161:320,]; b <- slm[1:160,]; slm <- rbind(a,b)

LIG_anom <- LIG_anom*LIG_mask
LGM_anom <- LGM_anom*slm
MH_anom <- MH_anom*slm


# Load lon and lat data

ncin <- nc_open("Echam5-wiso.T106_EXP010_PL_PI.global.monmean.d18O_prec.nc")
lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
lon1 <- lon[162:length(lon)]; lon2 <- lon[1:161]; lon <- c(lon1 - 360, lon2)
nlon <- length(lon); nlat <- length(lat)
nc_close(ncin); rm(ncin)


# STEP 2: Extract model data for speleothem sites

# Load sites
sites <- read.csv("boxplot_sites.csv")
sites <- sites[,-c(5,7)]; sites <- unique(sites)

# Reshape data
#ECHAM
MH_anom_df <- data.frame(expand.grid(lon = lon, lat = lat), pre = as.numeric(MH_anom))
LGM_anom_df <- data.frame(expand.grid(lon = lon, lat = lat), pre = as.numeric(LGM_anom))
LIG_anom_df <- data.frame(expand.grid(lon = lon, lat = lat), pre = as.numeric(LIG_anom))


# Calculate distance weighted mean for each site in each time period
# ECHAM
site_means_E <- list()
for (i in 1:nrow(sites)){
  site_means_E[[i]] <- data.frame(sites[i,], 
                                  Euclid_dist(x = sites$lon[i], y = sites$lat[i], mlength = degs_site, dt = MH_anom_df),
                                  Euclid_dist(x = sites$lon[i], y = sites$lat[i], mlength = degs_site, dt = LGM_anom_df),
                                  Euclid_dist(x = sites$lon[i], y = sites$lat[i], mlength = degs_site, dt = LIG_anom_df))
}
site_means_E <- rbindlist(site_means_E) # Combine list into dataframe
colnames(site_means_E)[6:8] <- c("MH","LGM","LIG") # Set column names

pre_Echam <- site_means_E %>% gather(key = time_slice, value = pre_anom, 6:8) # reshape so that time slices are a variable in df


# SAME STEPS AGAIN FOR TEMP DATA
rm(list= ls()[!(ls() %in% c('spel_d18O','wiso_Echam','pre_Echam','degs_site','Euclid_dist'))])

# STEP 1: load data and calculate anomalies

# Load (resampled) Echam ppt data
ECH <- readMat("ECHAM_dat.mat")

LIG_anom <- ECH$LIG.temp.anom
LGM_anom <- ECH$LGM.temp.anom
MH_anom <- ECH$MH.temp.anom

LIG_mask <- ECH$land.mask.LIG

slm <- ncvar_get(nc_open("T106.slm.w_miss.nc"), "slm")
a <- slm[161:320,]; b <- slm[1:160,]; slm <- rbind(a,b)

LIG_anom <- LIG_anom*LIG_mask
LGM_anom <- LGM_anom*slm
MH_anom <- MH_anom*slm


# Load lon and lat data

ncin <- nc_open("C:/Users/sarah/Documents/PhD/analyses/data/Echam5-wiso.T106_EXP010_PL_PI.global.monmean.d18O_prec.nc")
lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
lon1 <- lon[162:length(lon)]; lon2 <- lon[1:161]; lon <- c(lon1 - 360, lon2)
nlon <- length(lon); nlat <- length(lat)
nc_close(ncin); rm(ncin)



# STEP 2: Extract model data for speleothem sites

# Load sites
sites <- read.csv("boxplot_sites.csv")
sites <- sites[,-c(5,7)]; sites <- unique(sites)

# Reshape data
#ECHAM
MH_anom_df <- data.frame(expand.grid(lon = lon, lat = lat), temp = as.numeric(MH_anom))
LGM_anom_df <- data.frame(expand.grid(lon = lon, lat = lat), temp = as.numeric(LGM_anom))
LIG_anom_df <- data.frame(expand.grid(lon = lon, lat = lat), temp = as.numeric(LIG_anom))


# Calculate distance weighted mean for each site in each time period
# ECHAM
site_means_E <- list()
for (i in 1:nrow(sites)){
  site_means_E[[i]] <- data.frame(sites[i,], 
                                  Euclid_dist(x = sites$longitude[i], y = sites$latitude[i], mlength = degs_site, dt = MH_anom_df),
                                  Euclid_dist(x = sites$longitude[i], y = sites$latitude[i], mlength = degs_site, dt = LGM_anom_df),
                                  Euclid_dist(x = sites$longitude[i], y = sites$latitude[i], mlength = degs_site, dt = LIG_anom_df))
}
site_means_E <- rbindlist(site_means_E) # Combine list into dataframe
colnames(site_means_E)[6:8] <- c("MH","LGM","LIG") # Set column names

temp_Echam <- site_means_E %>% gather(key = time_slice, value = temp_anom, 6:8) # reshape so that time slices are a variable in df



# STEP 3: Combine speleothem and ECHAM into 1 dataframe

colnames(spel_d18O)[7] <- "val"
colnames(wiso_Echam)[7] <- "val"
colnames(pre_Echam)[7] <- "val"
colnames(temp_Echam)[7] <- "val"

spel_d18O$dat_type <- "spel_d18O"
wiso_Echam$dat_type <- "Echam_d18O"
pre_Echam$dat_type <- "Echam_pre"
temp_Echam$dat_type <- "Echam_temp"

dat_all <- rbind(wiso_Echam, pre_Echam, temp_Echam)
dat_all <- dat_all %>% filter(site_id != 104) # remove Liang Luar (is not located within land grids, so all vals = 0)
dat_all <- rbind(dat_all, spel_d18O)


# export data
write.csv(dat_all, "3_degs_boxplot_dat.csv", row.names = F)
