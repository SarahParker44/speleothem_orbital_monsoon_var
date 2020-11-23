### Suplementary figs: modelled d18Oprecip versus SISAL ###

library(R.matlab)
library(dplyr)
library(tidyr)
library(raster)
library(rgdal)
library(rasterVis)
library(RMySQL)
library(ncdf4)

## Connect to SISAL db for speleothem data

db_user <- 'sarah'
db_password <- ''
db_name <- 'sisalv2'
db_host <- '134.225.64.103' 

mydb <-  dbConnect(MySQL(), user = db_user, password = db_password,
                   dbname = db_name, host = db_host)  # Connect to MySQL db

dbExecute(mydb, 'SET NAMES UTF8;')

# Load spel data
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING(site_id) JOIN sample USING(entity_id) JOIN original_chronology USING(sample_id) JOIN d18o USING(sample_id) 
                       WHERE (interp_age BETWEEN 0 AND 127000);")

## ECHAM figs (MH, LGM, LIG)

for (i in c("MH","LGM","LIG")){
  
  # load processed data
  
  if (i %in% c("MH","LGM")){
    x <- readMat(paste("C:/Users/ph805612/OneDrive - University of Reading/Documents/MATLAB/", i, "_files.mat", sep = ""))
    wiso <- x[[paste("wiso",i,"ave", sep = ".")]]
    temp <- x[[paste("temp",i,"ave", sep = ".")]]
    x <- readMat(paste("C:/Users/ph805612/OneDrive - University of Reading/Documents/MATLAB/", i, "_PI_files.mat", sep = ""))
    wiso_PI <- x[[paste("wiso",i,"PI.ave", sep = ".")]]
    temp_PI <- x[[paste("temp",i,"PI.ave", sep = ".")]]
  } else {
    x <- readMat("C:/Users/ph805612/OneDrive - University of Reading/Documents/MATLAB/LIG.mat")
    wiso <- x[["wiso.LIG"]]
    temp <- x[["temp.LIG"]]
    x <- readMat("C:/Users/ph805612/OneDrive - University of Reading/Documents/MATLAB/LIG_PI.mat")
    wiso_PI <- x[["wiso.LIG.PI"]]
    temp_PI <- x[["temp.LIG.PI"]]
  }
  
  # calc as anom
  anom <- wiso - wiso_PI
  
  #reshape
  anom <- anom[,ncol(anom):1] # lat from 90->-90 to -90->90
  temp <- temp[,ncol(temp):1]; temp_PI <- temp_PI[,ncol(temp_PI):1]
  a <- anom[((nrow(anom)/2)+1):nrow(anom),]; b <- anom[1:(nrow(anom)/2),]; anom <- rbind(a,b) # convert lon from 0->360 to -180->180
  a <- temp[((nrow(temp)/2)+1):nrow(temp),]; b <- temp[1:(nrow(temp)/2),]; temp <- rbind(a,b)
  a <- temp_PI[((nrow(temp_PI)/2)+1):nrow(temp_PI),]; b <- temp_PI[1:(nrow(temp_PI)/2),]; temp_PI <- rbind(a,b)
  
  # extract SISAL data for i period and PI
  if (i == "MH"){ period_lims = c(5500, 6500) } else if (i == "LGM") { period_lims = c(20000,22000) } else {
    period_lims = c(124000,126000) } # period limits
    
  dat <- Raw_Data %>% 
    filter(interp_age >= -40 & interp_age <= 100 |
             interp_age >= period_lims[1] & interp_age <= period_lims[2]) %>%
    mutate(t_slice = ifelse(interp_age >= -40 & interp_age <= 100, "PI", "period")) %>%
    dplyr::select(site_id, latitude, longitude, mineralogy, arag_corr, t_slice, d18O_measurement) %>%
    group_by(site_id, latitude, longitude) %>%
    filter(length(unique(t_slice)) > 1) 
  
  dat_PI <- dat %>% filter(t_slice == "PI"); dat_period <- dat %>% filter(t_slice == "period")

  
  if (i %in% c("MH","LGM")){
    ncin <- nc_open("C:/Users/ph805612/OneDrive - University of Reading/Documents/SISAL/R_programming/Data/Echam5-wiso.T106_EXP010_PL_PI.global.monmean.d18O_prec.nc")
  } else {
    ncin <- nc_open("C:/Users/ph805612/OneDrive - University of Reading/Documents/SISAL/R_programming/Data/Eem125-S2_echam5_wiso_wisoaprt_d_ymonmean.nc")
  }
  lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
  
  #reshape
  lat <- lat[length(lat):1] # lat from 90->-90 to -90->90
  a <- (lon[((length(lon)/2)+1):length(lon)])-360; b <- lon[1:(length(lon)/2)]; lon <- c(a,b) # convert lon from 0->360 to -180->180
  
  ## Extract temp vals for each site
  #PI
  lon_no <- sapply(dat_PI$longitude, function(x) which.min(abs(lon - x)))
  lat_no <- sapply(dat_PI$latitude, function(x) which.min(abs(lat - x)))
  ECH_tempPI <- data.frame(expand.grid(lon = lon, lat = lat), temp = as.numeric(temp_PI))
  lon_lat <- (lat_no-1)*length(lon) + lon_no
  ECH_sub <- ECH_tempPI[lon_lat,]
  rownames(ECH_sub) <- 1:nrow(ECH_sub)
  dat_PI <- cbind(data.frame(dat_PI), data.frame(ECH_sub))
  
  # period
  lon_no <- sapply(dat_period$longitude, function(x) which.min(abs(lon - x)))
  lat_no <- sapply(dat_period$latitude, function(x) which.min(abs(lat - x)))
  ECH_temp <- data.frame(expand.grid(lon = lon, lat = lat), temp = as.numeric(temp))
  lon_lat <- (lat_no-1)*length(lon) + lon_no
  ECH_sub <- ECH_temp[lon_lat,]
  rownames(ECH_sub) <- 1:nrow(ECH_sub)
  dat_period <- cbind(data.frame(dat_period), data.frame(ECH_sub))
  
  ## combine and convert d18Ospel PDB to d18Odw SMOW
  dat_all <- rbind(dat_PI, dat_period)
  dat_all$d18O_SMOW <- 1.03092 * dat_all$d18O_measurement + 30.92 #PDB to SMOW
  dat_all$temp_K <- dat_all$temp + 273.15 #temp from celcius to kelvin
  dat_calcite <- dat_all %>% filter(mineralogy == "calcite" |
                                      arag_corr == "yes") %>%
    mutate(d18O_dw_SMOW = d18O_SMOW - (((16.1*1000)/temp_K)-24.6))
  dat_arag <- dat_all %>% filter(mineralogy == "aragonite" & arag_corr %in% c("no","unknown")) %>%
    mutate(d18O_dw_SMOW = d18O_SMOW - (((18.34*1000)/temp_K)-31.954))
  dat_all <- rbind(dat_calcite, dat_arag)
  
  # summarise signal of each spel site
  dat_all <- dat_all %>%
    group_by(site_id, latitude, longitude, t_slice) %>%
    summarise(mean_d18O = mean(d18O_dw_SMOW)) %>%
    spread(key = t_slice, value = mean_d18O, 4:5) %>%
    mutate(diff = period-PI, signal = ifelse(diff >= 0, "pos", "neg"),
           diff_upr = diff + 0.5, diff_lwr = diff -0.5,
           signal_upr = ifelse(diff_upr >= 0, "pos", "neg"),
           signal_lwr = ifelse(diff_lwr >= 0, "pos", "neg"))
  
  ## calc % agreement between SISAL and modelled d18Op
  # extract d18Oprecip for each site
  lon_no <- sapply(dat_all$longitude, function(x) which.min(abs(lon - x)))
  lat_no <- sapply(dat_all$latitude, function(x) which.min(abs(lat - x)))
  ECH_dat <- data.frame(expand.grid(lon = lon, lat = lat), wiso = as.numeric(anom))
  lon_lat <- (lat_no-1)*length(lon) + lon_no
  ECH_sub <- ECH_dat[lon_lat,]
  rownames(ECH_sub) <- 1:nrow(ECH_sub)
  dat_all <- cbind(data.frame(dat_all), data.frame(ECH_sub))

  n_neg_con <- nrow(dat_all %>% filter(wiso <= 0) %>% filter(signal == "neg" | signal_lwr == "neg"))
  n_pos_con <- nrow(dat_all %>% filter(wiso >= 0) %>% filter(signal == "pos" | signal_upr == "pos"))
  
  #n_con <- nrow(dat_all %>% filter(signal == "neg" & wiso <= 0 | signal == "pos" & wiso >= 0))
  n_con <- n_neg_con+n_pos_con
  perc <- n_con/nrow(dat_all)
  
  # plot map
  
  # convert spel point data to coordinates df (for plotting)
  xy <- data.frame(x= dat_all$longitude, y = dat_all$latitude, z = dat_all$signal)
  coordinates(xy) <- ~x+y
  
  wiso <- raster(t(anom)[ncol(anom):1, ])
  projection(wiso) <- CRS("+init=epsg:4326")
  extent(wiso) <- c(min(-180), max(180), min(-90), max(90))
  
  coast_lines <- readOGR("C:/Users/ph805612/OneDrive - University of Reading/Documents/SISAL/R_programming/Data/ne_110m_land/ne_110m_land.shp", verbose = F)
  color.ramp.length <- 100
  neg.length <- round(abs(min(anom))/
                        (max(anom)-min(anom)) *
                        color.ramp.length)
  pos.length <- color.ramp.length - neg.length
  
  cols <- c(colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7"))(neg.length), 
            colorRampPalette(c("#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))(pos.length))
  
  png(paste("C:/Users/ph805612/OneDrive - University of Reading/Documents/monsoon_paper/final_figures/SI_", i, "_ECH_comp_map.png", sep = ""), 
      width = 18, height = 10, units = "cm", res = 300)
  p1 <- levelplot(wiso, margin = F, col.regions = cols, at = seq(min(anom), max(anom), (max(anom)-min(anom))/100),
            xlab = NULL, ylab = NULL) +
    latticeExtra::layer(sp.lines(coast_lines)) +
    latticeExtra::layer({grid.rect(x = 0.08, y = 0.1, height = 0.07, width = 0.1)
                        grid.text(paste(round(perc, 2)*100, "%"), x = 0.08, y = 0.1)})
  p2 <- spplot(xy, col.regions = c("#053061","#67001F"), pch = c(1,3))
  p <- p1 + p2
  l1 <- p1$legend$right
  l2 <- p2$legend$bottom
  ll <- mergedTrellisLegendGrob(l1, l2, vertical = T)
  p$legend$right$fun <- ll
  print(p)
  dev.off()
}



## Same again for GISS

for (i in c(3,6,9)){
  
  # load GISS data
  ncin <- nc_open("C:/Users/ph805612/OneDrive - University of Reading/Documents/allegra_fortran/nc_files/0ka_wiso_100yr.nc")
  
  wiso_PI <- ncvar_get(ncin, "H2O18_in_prec")
  
  lon <- ncvar_get(ncin, "lon")
  lat <- ncvar_get(ncin, "lat")
  
  nc_close(ncin); rm(ncin)
  
  ncin <- nc_open(paste("C:/Users/ph805612/OneDrive - University of Reading/Documents/allegra_fortran/nc_files/", i, "ka_wiso_100yr.nc", sep = ""))
  
  wiso <- ncvar_get(ncin, "H2O18_in_prec")
  
  nc_close(ncin)
  
  # calc as anom
  anom = wiso - wiso_PI
  
  # tsurf data
  ncin <- nc_open("C:/Users/ph805612/OneDrive - University of Reading/Documents/allegra_fortran/GISS_clim_vars.nc")
  temp_PI <- ncvar_get(ncin, "tsurf")[,,1,]; temp_PI <- apply(temp_PI, c(1,2), mean)
  x <- ifelse(i == 3, 4, ifelse(i == 6, 7, 8))
  temp_period <- ncvar_get(ncin, "tsurf")[,,x,]; temp_period <- apply(temp_period, c(1,2), mean, na.rm = T)
  
  # extract SISAL data for i period and PI
  if (i == 3){ period_lims = c(2500, 3500) } else if (i == 6) { period_lims = c(5500,6500) } else {
    period_lims = c(8500,9500) } # period limits
  
  dat <- Raw_Data %>% 
    filter(interp_age >= -40 & interp_age <= 100 |
             interp_age >= period_lims[1] & interp_age <= period_lims[2]) %>%
    mutate(t_slice = ifelse(interp_age >= -40 & interp_age <= 100, "PI", "period")) %>%
    dplyr::select(site_id, latitude, longitude, mineralogy, arag_corr, t_slice, d18O_measurement) %>%
    group_by(site_id, latitude, longitude) %>%
    filter(length(unique(t_slice)) > 1) 
  
  dat_PI <- dat %>% filter(t_slice == "PI"); dat_period <- dat %>% filter(t_slice == "period")
  
  ## Extract temp vals for each site
  #PI
  lon_no <- sapply(dat_PI$longitude, function(x) which.min(abs(lon - x)))
  lat_no <- sapply(dat_PI$latitude, function(x) which.min(abs(lat - x)))
  GISS_tempPI <- data.frame(expand.grid(lon = lon, lat = lat), temp = as.numeric(temp_PI))
  lon_lat <- (lat_no-1)*length(lon) + lon_no
  GISS_sub <- GISS_tempPI[lon_lat,]
  rownames(GISS_sub) <- 1:nrow(GISS_sub)
  dat_PI <- cbind(data.frame(dat_PI), data.frame(GISS_sub))
  
  # period
  lon_no <- sapply(dat_period$longitude, function(x) which.min(abs(lon - x)))
  lat_no <- sapply(dat_period$latitude, function(x) which.min(abs(lat - x)))
  GISS_temp <- data.frame(expand.grid(lon = lon, lat = lat), temp = as.numeric(temp_period))
  lon_lat <- (lat_no-1)*length(lon) + lon_no
  GISS_sub <- GISS_temp[lon_lat,]
  rownames(GISS_sub) <- 1:nrow(GISS_sub)
  dat_period <- cbind(data.frame(dat_period), data.frame(GISS_sub))
  
  ## combine and convert d18Ospel PDB to d18Odw SMOW
  dat_all <- rbind(dat_PI, dat_period)
  dat_all$d18O_SMOW <- 1.03092 * dat_all$d18O_measurement + 30.92 #PDB to SMOW
  dat_all$temp_K <- dat_all$temp + 273.15 #temp from celcius to kelvin
  dat_calcite <- dat_all %>% filter(mineralogy == "calcite" |
                                      arag_corr == "yes") %>%
    mutate(d18O_dw_SMOW = d18O_SMOW - (((16.1*1000)/temp_K)-24.6))
  dat_arag <- dat_all %>% filter(mineralogy == "aragonite" & arag_corr %in% c("no","unknown")) %>%
    mutate(d18O_dw_SMOW = d18O_SMOW - (((18.34*1000)/temp_K)-31.954))
  dat_all <- rbind(dat_calcite, dat_arag)
  
  # summarise signal of each spel site
  dat_all <- dat_all %>%
    group_by(site_id, latitude, longitude, t_slice) %>%
    summarise(mean_d18O = mean(d18O_dw_SMOW)) %>%
    spread(key = t_slice, value = mean_d18O, 4:5) %>%
    mutate(diff = period-PI, signal = ifelse(diff >= 0, "pos", "neg"),
           diff_upr = diff + 0.5, diff_lwr = diff -0.5,
           signal_upr = ifelse(diff_upr >= 0, "pos", "neg"),
           signal_lwr = ifelse(diff_lwr >= 0, "pos", "neg"))

  
  ## calc % agreement between SISAL and modelled d18Op
  lon_no <- sapply(dat_all$longitude, function(x) which.min(abs(lon - x)))
  lat_no <- sapply(dat_all$latitude, function(x) which.min(abs(lat - x)))
  
  GISS_dat <- data.frame(expand.grid(lon = lon, lat = lat), wiso = as.numeric(anom))
  
  lon_lat <- (lat_no-1)*length(lon) + lon_no
  GISS_sub <- GISS_dat[lon_lat,]
  rownames(GISS_sub) <- 1:nrow(GISS_sub)
  
  dat_all <- cbind(data.frame(dat_all), data.frame(GISS_sub))
  
  n_neg_con <- nrow(dat_all %>% filter(wiso <= 0) %>% filter(signal == "neg" | signal_lwr == "neg"))
  n_pos_con <- nrow(dat_all %>% filter(wiso >= 0) %>% filter(signal == "pos" | signal_upr == "pos"))
  
  #n_con <- nrow(dat_all %>% filter(signal == "neg" & wiso <= 0 | signal == "pos" & wiso >= 0))
  n_con <- n_neg_con+n_pos_con
  
  perc <- n_con/nrow(dat_all)
  
  # plot map
  
  # convert spel point data to coordinates df (for plotting)
  xy <- data.frame(x= dat_all$longitude, y = dat_all$latitude, z = dat_all$signal)
  coordinates(xy) <- ~x+y
  
  wiso <- raster(t(anom)[ncol(anom):1, ])
  projection(wiso) <- CRS("+init=epsg:4326")
  extent(wiso) <- c(min(-180), max(180), min(-90), max(90))
  
  coast_lines <- readOGR("C:/Users/ph805612/OneDrive - University of Reading/Documents/SISAL/R_programming/Data/ne_110m_land/ne_110m_land.shp", verbose = F)
  color.ramp.length <- 100
  neg.length <- round(abs(min(anom))/
                        (max(anom)-min(anom)) *
                        color.ramp.length)
  pos.length <- color.ramp.length - neg.length
  
  cols <- c(colorRampPalette(c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7"))(neg.length), 
            colorRampPalette(c("#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))(pos.length))
  
  png(paste("C:/Users/ph805612/OneDrive - University of Reading/Documents/monsoon_paper/final_figures/SI_", i, "ka_GISS_comp_map.png", sep = ""), 
      width = 18, height = 10, units = "cm", res = 300)
  p1 <- levelplot(wiso, margin = F, col.regions = cols, at = seq(min(anom), max(anom), (max(anom)-min(anom))/100),
                  xlab = NULL, ylab = NULL) +
    latticeExtra::layer(sp.lines(coast_lines)) +
    latticeExtra::layer({grid.rect(x = 0.08, y = 0.1, height = 0.07, width = 0.1)
      grid.text(paste(round(perc, 2)*100, "%"), x = 0.08, y = 0.1)})
  #layer(sp.points(xy, pch = ifelse(xy$z == "pos", 3, 1), col = ifelse(xy$z == "pos","#67001F","#053061")))
  p2 <- spplot(xy, col.regions = ifelse(xy$z == "pos","#67001F","#053061"), pch = c(1,3))
  p <- p1 + p2
  l1 <- p1$legend$right
  l2 <- p2$legend$bottom
  ll <- mergedTrellisLegendGrob(l1, l2, vertical = T)
  p$legend$right$fun <- ll
  print(p)
  dev.off()
}
