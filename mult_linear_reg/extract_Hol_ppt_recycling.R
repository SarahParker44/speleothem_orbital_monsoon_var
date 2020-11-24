### Calculate precipitation recycling index for each region and time slice ####

setwd("C:/Users/sarah/Documents/PhD/analyses/")

# load necessary packages

library(ncdf4)
library(dplyr)
library(tidyr)
library(geosphere)

# Load mask

ncin <- nc_open("data/GISS_landmask.nc")
mask <- ncvar_get(ncin, "frac_land")
mask <- ifelse(mask > 50, 1, NA) # 1  = land grids, NA = ocean grids

# load lon and lat data
lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")

# Load grid surface area data
g_area <- ncvar_get(ncin, "axyp") # surface area of grids in m^2
nc_close(ncin); rm(ncin)

# calculate gridbox length (lon) and width (lat) for each grid using distGeo
# output = global matrix of lon lengths and lat lengths
#lon_mat <- matrix(rep(lon, 46), nrow = 72, ncol = 46)
#lat_mat <- matrix(rep(lat, each = 72), nrow = 72, ncol = 46)

lonlen <- matrix(NA, nrow = 72, ncol = 46)
latlen <- matrix(NA, nrow = 72, ncol = 46)

for (i in 1:72){ # each lon
  for (j in 1:46){ # each lat
    
    lon_val <- lon[i]
    lat_val <- lat[j]
    
    lon_int <- lon[3] - lon[2] # get lon resolution
    lat_int <- lat[3] - lat[2] # get lat resolution
    
    if (!j %in% c(1,46)){
      lon_len <- distGeo(p1 = c(lon_val-(lon_int/2), lat_val), p2 = c(lon_val+(lon_int/2), lat_val)) #get lon length of grids in m
      lat_len <- distGeo(p1 = c(lon_val,lat_val-(lat_int/2)), p2 = c(lon_val, lat_val+(lat_int/2))) #get lat length of grids in m 
    } else if (j %in% c(1,46)){ 
      lon_len <- distGeo(p1 = c(lon_val-(lon_int/2), lat_val), p2 = c(lon_val+(lon_int/2), lat_val))
      lat_len <- distGeo(p1 = c(lon_val,lat_val-1), p2 = c(lon_val,lat_val+1))
    }
    lonlen[i,j] <- lon_len
    latlen[i,j] <- lat_len
  }
}

# Load GISS moisture flux and evap
# dim = 72 (lon) x 46 (lat) x 8 (time slices) x 120 (monthly means for 10 decades)
ncin <- nc_open("data/GISS_clim_vars.nc")
mois_NS <- ncvar_get(ncin, "pvq") # mb*m/s
mois_EW <- ncvar_get(ncin, "puq") # mb*m/s

evap <- ncvar_get(ncin, "evap") # mm/day

nc_close(ncin); rm(ncin)


# define regions limits

region_lims <- data.frame(region = c('ISM','EAM','SW-SAM1','SW-SAM2','IAM','NE-SAM','CAM','SAfM'),
                          min_lat = c(11,20,-10,-30,-24,-10,10,-30),
                          max_lat = c(32,39,0,-10,5,0,33,-17),
                          min_lon = c(50,100,-80,-68,95,-60,-115,10),
                          max_lon = c(95,125,-64,-40,135,-30,-58,40))


# arrays to fill with calculated ppt recycling for each time dimension and region
dat_prec <- array(NA, dim = c(8,120,7))

# Calculate boundary total moisture flux and region total evap for each time dimension and region:
for (i in 1:8){ # each time slice
  for (j in 1:120){ # each monthly mean
    
    # subset to time slice and month, apply mask (ocean grids to NA)
    moisNS <- mois_NS[,,i,j] * mask
    moisEW <- mois_EW[,,i,j] * mask
    evp <- evap[,,i,j] * mask
    
    # mois flux unit conversion (mb*m/s -> mm/day)
    
    moisNS <- moisNS*10.197; moisNS <- moisNS*86400; # mb->mm, s->day
    moisEW <- moisEW*10.197; moisEW <- moisEW*86400;
    
    for (k in c(1:2,5:8)){ # regions (exluding SW-SAM)
      
      # Select regional data limits
      lon_reg <- region_lims[k, 4:5]
      lat_reg <- region_lims[k, 2:3]
      
      # matrix values for subsetting to region
      region_lon <- sapply(lon_reg, function(x) which.min(abs(lon-x)))
      region_lat <- sapply(lat_reg, function(x) which.min(abs(lat-x)))
      
      ## Calculate total evap over region
      # subset evap data to region
      evap_reg <- evp[region_lon[1]:region_lon[2],region_lat[1]:region_lat[2]]
      # subset grid area to region
      g_area_reg <- g_area[region_lon[1]:region_lon[2],region_lat[1]:region_lat[2]]
      # multiply evap in each grid box by surface area of grid box
      evap_reg_x <- evap_reg*g_area_reg 
      # sum vals
      tot_evap <- sum(evap_reg_x, na.rm = T)
      
      
      ## Calculate moisture flux across each boundary
      
      # select N,S,E,W margin grids (won't be straight along coast lines)
      
      moisEW_reg <- moisEW[region_lon[1]:region_lon[2],region_lat[1]:region_lat[2]] # subset to region
      moisEW_reg <- rbind(matrix(NA, ncol = ncol(moisEW_reg)), moisEW_reg, matrix(NA, ncol = ncol(moisEW_reg))) # add NAs around outside of region subset
      moisEW_reg <- cbind(matrix(NA, nrow=nrow(moisEW_reg)), moisEW_reg, matrix(NA, nrow=nrow(moisEW_reg)))
      
      moisNS_reg <- moisNS[region_lon[1]:region_lon[2],region_lat[1]:region_lat[2]]
      moisNS_reg <- rbind(matrix(NA, ncol = ncol(moisNS_reg)), moisNS_reg, matrix(NA, ncol = ncol(moisNS_reg))) 
      moisNS_reg <- cbind(matrix(NA, nrow=nrow(moisNS_reg)), moisNS_reg, matrix(NA, nrow=nrow(moisNS_reg)))
      
      lonlen_reg <- lonlen[region_lon[1]:region_lon[2],region_lat[1]:region_lat[2]]
      lonlen_reg <- rbind(matrix(NA, ncol = ncol(lonlen_reg)), lonlen_reg, matrix(NA, ncol = ncol(lonlen_reg))) # add NAs around outside of land evap vals
      lonlen_reg <- cbind(matrix(NA, nrow=nrow(lonlen_reg)), lonlen_reg, matrix(NA, nrow=nrow(lonlen_reg)))
      
      latlen_reg <- latlen[region_lon[1]:region_lon[2],region_lat[1]:region_lat[2]]
      latlen_reg <- rbind(matrix(NA, ncol = ncol(latlen_reg)), latlen_reg, matrix(NA, ncol = ncol(latlen_reg))) # add NAs around outside of land evap vals
      latlen_reg <- cbind(matrix(NA, nrow=nrow(latlen_reg)), latlen_reg, matrix(NA, nrow=nrow(latlen_reg)))

      
      # get grid positions for each border grid
      
      S_border <- numeric()
      N_border <- numeric()
      E_border <- numeric()
      W_border <- numeric()
      
      vec_pos <- 1
      for (a in 1:nrow(moisNS_reg)){
        for (b in 2:(ncol(moisNS_reg))){
          if (is.na(moisNS_reg[a,b]) == F & is.na(moisNS_reg[a,b-1]) == T){ # S margin = left border of data in matrix, so vals where grid != NA, but grid in column to left = NA
            val <- moisNS_reg[a,b]*lonlen_reg[a,b] # multiply mois flux by lon length
            S_border[vec_pos] <- val # save value
            vec_pos <- vec_pos + 1
          } 
        }
      }
      S_border <- sum(S_border) # Flux over entire border
      vec_pos <- 1
      for (a in 1:nrow(moisNS_reg)){
        for (b in 1:(ncol(moisNS_reg)-1)){
          if (is.na(moisNS_reg[a,b]) == F & is.na(moisNS_reg[a,b+1]) == T){ # N margin = right border of data in matrix, so where grid !-NA, but grid to right = NA
            val <- moisNS_reg[a,b]*lonlen_reg[a,b] # multiply mois flux by lat length
            N_border[vec_pos] <- val # save value
            vec_pos <- vec_pos + 1
          }
        }
      }
      N_border <- sum(N_border)
      vec_pos <- 1
      for (a in 1:(nrow(moisEW_reg)-1)){
        for (b in 1:ncol(moisEW_reg)){
          if (is.na(moisEW_reg[a,b]) == F & is.na(moisEW_reg[a+1,b]) == T){ # E margin = bottom border of data in matrix, so where grid != NA, but grid below = NA
            val <- moisEW_reg[a,b]*latlen_reg[a,b] # multiply mois flux by lat length
            E_border[vec_pos] <- val # save val
            vec_pos <- vec_pos + 1
          } 
        }
      }
      E_border <- sum(E_border)
      vec_pos <- 1
      for (a in 2:(nrow(moisEW_reg))){
        for (b in 1:ncol(moisEW_reg)){
          if (is.na(moisEW_reg[a,b]) == F & is.na(moisEW_reg[a-1,b]) == T){ # W margin = top border of data in matrix, so where grid != NA, but grid box above = NA
            val <- moisEW_reg[a,b]*latlen_reg[a,b] # multiply mois flux by lat length
            W_border[vec_pos] <- val # save val
            vec_pos <- vec_pos + 1
          }
        }
      }
      W_border <- sum(W_border)
      
      # Incoming fluxes
      if (N_border >0){ N_border = 0}
      if (S_border <0){ S_border = 0}
      if (E_border >0){ E_border = 0}
      if (W_border <0){ W_border = 0}
      
      # Calculate flux
      Flux <- W_border - E_border + S_border - N_border
      
      # Calculate precipitation recycling
      prec = tot_evap/(tot_evap + 2*Flux);
      
      # Save val
      reg_val <- ifelse(k %in% 1:2, k, k-1)
      dat_prec[i,j, reg_val] <- prec
      
      rm(lon_reg, lat_reg, region_lon, region_lat, evap_reg, g_area_reg, evap_reg_x, tot_evap, moisEW_reg, moisNS_reg, 
         lonlen_reg, latlen_reg, N_border, S_border, E_border, W_border, Flux, prec)
    
    }

    
    # Same again, but for SW-SAM (2 rectangles together)
    
    # Get region limits
    lon_reg_SAM1 <- region_lims[3, 4:5]; lon_reg_SAM2 <- region_lims[4, 4:5]
    lat_reg_SAM1 <- region_lims[3, 2:3]; lat_reg_SAM2 <- region_lims[4, 2:3]
    
    # Matrix values for subsetting to region
    SAM1_lon <- sapply(lon_reg_SAM1, function(x) which.min(abs(lon-x))); SAM2_lon <- sapply(lon_reg_SAM2, function(x) which.min(abs(lon-x)))
    SAM1_lat <- sapply(lat_reg_SAM1, function(x) which.min(abs(lat-x))); SAM2_lat <- sapply(lat_reg_SAM2, function(x) which.min(abs(lat-x)))
    
    # subset to region
    df <- data.frame(expand.grid(long = lon, lati = lat), evap = as.numeric(evp), mois_EW = as.numeric(moisEW), mois_NS = as.numeric(moisNS), lonlen = as.numeric(lonlen), latlen = as.numeric(latlen), g_area = as.numeric(g_area))
    for (l in 1:nrow(df)){
      ln <- df$long[l]
      lt <- df$lati[l]
      if (ln %in% lon[SAM1_lon[1]:SAM1_lon[2]] & lt %in% lat[SAM1_lat[1]:SAM1_lat[2]] | 
          ln %in% lon[SAM2_lon[1]:SAM2_lon[2]] & lt %in% lat[SAM2_lat[1]:SAM2_lat[2]]){
        df[l,3:8] <- df[l,3:8] 
      } else {
        df[l,3:8] <- NA
      }
    }
    
    # df to matrix
    evap_SAM <- matrix(df$evap, nrow = nrow(evap), ncol = ncol(evap))
    moisNS_SAM <- matrix(df$mois_NS, nrow = nrow(mois_NS), ncol = ncol(mois_NS))
    moisEW_SAM <- matrix(df$mois_EW, nrow = nrow(mois_EW), ncol = ncol(mois_EW))
    lonlen_SAM <- matrix(df$lonlen, nrow = nrow(lonlen), ncol = ncol(lonlen))
    latlen_SAM <- matrix(df$latlen, nrow = nrow(latlen), ncol = ncol(latlen))
    g_area_SAM <- matrix(df$g_area, nrow = nrow(g_area), ncol = ncol(g_area))
    
    # Zoom to SAM
    evap_SAM <- evap_SAM[19:29,15:24]; 
    moisNS_SAM <- moisNS_SAM[19:29,15:24]; 
    moisEW_SAM <- moisEW_SAM[19:29,15:24]; 
    lonlen_SAM <- lonlen_SAM[19:29,15:24]; 
    latlen_SAM <- latlen_SAM[19:29,15:24]; 
    g_area_SAM <- g_area_SAM[19:29,15:24]; 
    
    # Calculate total evap over region
    evap_SAM_x <- evap_SAM*g_area_SAM # multiply evap in each grid box with surface area of grid box
    tot_evap <- sum(evap_SAM_x, na.rm = T) # sum vals
    
    # Calculate moisture flux across each boundary
    
    S_border <- numeric()
    N_border <- numeric()
    E_border <- numeric()
    W_border <- numeric()
    
    vec_pos <- 1
    for (a in 1:nrow(moisNS_SAM)){
      for (b in 2:(ncol(moisNS_SAM))){
        if (is.na(moisNS_SAM[a,b]) == F & is.na(moisNS_SAM[a,b-1]) == T){ # S margin = left border of data in matrix, so vals where grid != NA, but grid in same row but previous column = NA
          val <- moisNS_SAM[a,b]*lonlen_SAM[a,b]
          S_border[vec_pos] <- val
          vec_pos <- vec_pos + 1
        } 
      }
    }
    S_border <- sum(S_border)
    vec_pos <- 1
    for (a in 1:nrow(moisNS_SAM)){
      for (b in 1:(ncol(moisNS_SAM)-1)){
        if (is.na(moisNS_SAM[a,b]) == F & is.na(moisNS_SAM[a,b+1]) == T){ # N margin = right border of data in matrix, so where grid !-NA, but grid to right = NA
          val <- moisNS_SAM[a,b]*lonlen_SAM[a,b]
          N_border[vec_pos] <- val
          vec_pos <- vec_pos + 1
        }
      }
    }
    N_border <- sum(N_border)
    vec_pos <- 1
    for (a in 1:(nrow(moisEW_SAM)-1)){
      for (b in 1:ncol(moisEW_SAM)){
        if (is.na(moisEW_SAM[a,b]) == F & is.na(moisEW_SAM[a+1,b]) == T){ # E margin = bottom border of data in matrix, so where grid != NA, but grid below = NA
          val <- moisEW_SAM[a,b]*latlen_SAM[a,b]
          E_border[vec_pos] <- val
          vec_pos <- vec_pos + 1
        } 
      }
    }
    E_border <- sum(E_border)
    vec_pos <- 1
    for (a in 2:(nrow(moisEW_SAM))){
      for (b in 1:ncol(moisEW_SAM)){
        if (is.na(moisEW_SAM[a,b]) == F & is.na(moisEW_SAM[a-1,b]) == T){ # W margin = top border of data in matrix, so where grid != NA, but grid box above = NA
          val <- moisEW_SAM[a,b]*latlen_SAM[a,b]
          W_border[vec_pos] <- val
          vec_pos <- vec_pos + 1
        }
      }
    }
    W_border <- sum(W_border)
    
    # Influx vals
    if (N_border >0){ N_border = 0}
    if (S_border <0){ S_border = 0}
    if (E_border >0){ E_border = 0}
    if (W_border <0){ W_border = 0}
    
    # Calculate overall flux
    Flux <- W_border - E_border + S_border - N_border
    
    # Calculate precipitation recycling
    prec = tot_evap/(tot_evap + 2*Flux);
    
    # save val
    dat_prec[i,j, 3] <- prec
    
  }
}

# Calculate monthly mean ppt recycling
ppt_rec_df <- data.frame(region = c('ISM','EAM','SW-SAM','IAM','NE-SAM','CAM','SAfM'),
                         year = rep(rep(c(0:6,9), each = 12),7),
                         month = rep(1:12, 56),
                         ppt_recycling = NA) # empty df to fill with monmean ppt rec 

for (i in 1:8){ # time slice
  for (j in 1:7){ # region
    for (k in 1:12){ # month
      vals <- numeric()
      
      for (l in 1:10){ vals[l] <- (12*l)-(12-k)}
      
      monmean <- mean(dat_prec[i,vals,j], na.rm = T)
      
      reg_name <- c('ISM','EAM','SW-SAM','IAM','NE-SAM','CAM','SAfM')[j]
      tslice <- ifelse(i %in% 1:7, i-1, 9)
      
      ppt_rec_df[with(ppt_rec_df, which(region ==  reg_name & year == tslice & month == k)),"ppt_recycling"] <- monmean
    }
  }
}


# Calculate summer averages

# MJJAS (NH summer)
MJJAS_ppt_rec <- ppt_rec_df %>% filter(region %in% c("CAM","ISM","EAM") & month %in% 5:9) %>% group_by(region, year) %>% summarise(ppt_rec = mean(ppt_recycling))
# NDJFM (SH summer)
NDJFM_ppt_rec <- ppt_rec_df %>% filter(region %in% c("SW-SAM","NE-SAM","SAfM","IAM") & month %in% c(11:12,1:3)) %>% group_by(region, year) %>% summarise(ppt_rec = mean(ppt_recycling))

seas_ppt_rec <- rbind(MJJAS_ppt_rec, NDJFM_ppt_rec)


# Calculate anomalies
seas_ppt_rec <- seas_ppt_rec %>% spread(key = year, value = ppt_rec)

for (i in 1:nrow(seas_ppt_rec)){ # each region
  for (j in 3:ncol(seas_ppt_rec)){ # each time slice
    seas_ppt_rec[i,j] <- seas_ppt_rec[i,j] - seas_ppt_rec[i,2] # calc anom
  }
}

seas_ppt_rec <- seas_ppt_rec[,-2] %>% gather(key = "t_slice", value = "ppt_rec", 2:8)

# Save data
write.csv(seas_ppt_rec, "data/multiple_regression/region_ppt_recycling.csv", row.names = F)
