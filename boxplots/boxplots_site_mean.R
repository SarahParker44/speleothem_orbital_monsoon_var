# Boxplots showing magnitudes of change
# STEP 1: Load SISAL data and filter data
# STEP 2: calculate anomalies (to modern d18O (PDB calcite))
# STEP 3: save data

# Load packages

library(RMySQL)
library(dplyr)
library(ggplot2)
library(ncdf4)
library(tidyr)



# STEP 1: Extract data and apply filters

# Connect to database

db_user <- 'root'
db_password <- ''
db_name <- 'sisalv2'
db_host <- 'localhost'

mydb <-  dbConnect(MySQL(), user = db_user, password = db_password,
                   dbname = db_name, host = db_host)  # Connect to MySQL db

dbExecute(mydb, 'SET NAMES UTF8;')


# Extract data from present to LIG

Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING(site_id) JOIN sample USING(entity_id) JOIN original_chronology USING(sample_id) JOIN d18O USING(sample_id) 
                       WHERE (interp_age <= 128000) AND (latitude BETWEEN -30 AND 45);")  # specify the time range and region that I'm interested in within query


# filter to monsoon regions
Raw_Data$region <- with(Raw_Data, ifelse(latitude >= 10 & latitude <= 35 & longitude >= 50 & longitude <= 100, "ISM",
                                         ifelse(latitude >= 20 & latitude <= 35 & longitude >= 100 & longitude <= 125, "EAM", 
                                                ifelse(latitude >= -10 & latitude <= 0 & longitude >= -80 & longitude <= -70 | latitude >= -30 & latitude <= -10 & longitude >= -80 & longitude <= -40, "SW-SAM",
                                                       ifelse(latitude >= -30 & latitude <= 5 & longitude >= 80 & longitude <= 170, "IAM",
                                                              ifelse(latitude >= -30 & latitude <= 0 & longitude >= 0 & longitude <= 50, "SAfM", 
                                                                     ifelse(latitude >= 10 & latitude <= 25 & longitude >= -100 & longitude <= -58, "CAM", 
                                                                            ifelse(latitude >= -10 & latitude <= 0 & longitude >= -60 & longitude <= -30, "NE-SAM", "other"))))))))


dt <- Raw_Data %>% filter(region != "other")

# filter to current entities
dt <- dt %>% filter(entity_status != "superseded")

# Filter to the time periods
dt <- dt %>% filter(interp_age > 5500 & interp_age < 6500 | # MH = 6ka +/-500yrs
                      interp_age > 20000 & interp_age < 22000 | # LGM = 21ka +/-1000yrs
                      interp_age > 124000 & interp_age < 126000) # LIG = 125ka +/-2000yrs

# Filter to remove high altitude sites (>3500 m asl)
dt <- dt %>% filter(elevation <= 3500)

# REMOVE INCONSISTENT MINERALOGIES
# Identify calcite records
calcite <- dt %>% group_by(site_name, entity_id, mineralogy, arag_corr) %>% 
  filter(mineralogy == "calcite") %>%
  summarise(n()) #n=85

# Identify aragonite (corrected records)
arag_corr <- dt %>% group_by(site_name, entity_id, mineralogy, arag_corr) %>% 
  filter(mineralogy == "aragonite" & arag_corr == "yes") %>%
  summarise(n()) #n=4

x <- rbind(calcite, arag_corr)

# Identify aragonite (uncorrected) records
arag <- dt %>% group_by(site_name, entity_id, mineralogy, arag_corr) %>% 
  filter(mineralogy == "aragonite" & arag_corr == "no") %>%
  summarise(n()) #n=9

# Identify mixed records
mixed <- dt %>% group_by(site_name, entity_id, mineralogy, arag_corr) %>% 
  filter(mineralogy == "mixed") %>%
  summarise(n()) #n=1

# Identify unknown records
unknown <- dt %>% group_by(site_name, entity_id, mineralogy, arag_corr) %>% 
  filter(mineralogy == "unknown") %>%
  summarise(n()) #n=9


# Select entities with variable mineralogy
var_min <- na.omit(full_join(x, arag, by = c("site_name","entity_id"))) #0

# Remove varied mineralogy records from analysis
#dt <- dt %>% filter(!entity_id %in% var_min$entity_id)


# Filter to remove (uncorrected) mixed and unknown mineralogy records
dt <- dt %>% filter(mineralogy != "unknown")
dt[which(dt$mineralogy == "mixed" & dt$arag_corr != "yes"),] #0


# Add time slice column
dt$time_slice <- with(dt, ifelse(interp_age > 5500 & interp_age < 6500, "MH", 
                                        ifelse(interp_age > 20000 & interp_age < 22000, "LGM",
                                               ifelse(interp_age > 124000 & interp_age < 126000, "LIG", "other"))))


# Summarise number of sites in each region's time slice
n_sites <- dt %>% group_by(site_id, site_name, region, time_slice) %>% summarise(n()) %>%
  group_by(region, time_slice) %>% summarise(n())

# filter to region that has 3 time slices
regions <- n_sites %>% group_by(region) %>% summarise(n()) %>%
  filter(`n()` == 3)
dt <- dt %>% filter(region %in% regions$region)

# Remove SW-SAM (models don't get d18O right here)
dt <- filter(dt, region != "SW-SAM")

# Summarise (and save) sites
dt_sites <- dt %>% group_by(site_id, site_name, entity_id, latitude, longitude, mineralogy, region) %>% summarise(n())
write.csv(dt_sites, "C:/Users/sarah/Documents/PhD/analyses/data/boxplot_sites.csv", row.names = F)


# STEP 2: Calculate d18O as anomalies

# Extract OIPC d18O-ppt for site locations
# Extract temp data from CRU-ts for site locations
# Use temp to convert OIPC d18O from drip water to calcite, convert from SMOW to PDB
# Calculate anomalies


# Extract OIPC d18O for each site
ncin <- nc_open("C:/Users/sarah/Documents/PhD/analyses/data/oipc_global10_v2.nc")

lon <- ncvar_get(ncin, "lon"); nlon <- length(lon)
lat <- ncvar_get(ncin, "lat"); nlat <- length(lat)
d18O <- ncvar_get(ncin, "oipc_global10_v2")

nc_close(ncin); rm(ncin)

grid <- expand.grid(lon=lon, lat=lat)
j <- sapply(dt_sites$longitude, function(x) which.min(abs(lon-x))) 
k <- sapply(dt_sites$latitude, function(x) which.min(abs(lat-x)))
d18O_vec <- as.vector(d18O[,,26]) # select annual weighted d18O
jk <- (k-1)*nlon + j
oipc_d18O <- d18O_vec[jk]

pts <- data.frame(dt_sites$site_id, dt_sites$site_name, dt_sites$longitude,
                  dt_sites$latitude, dt_sites$mineralogy, oipc_d18O)
colnames(pts) <- c("site_id","site_name","longitude","latitude","mineralogy","OIPC_d18O")

# NA sites (site location just within ocean grids in oipc, manually find closest value)
pts[31,6] <- -2.997487 # Cave C126
pts[38,6] <- -0.8775846 # Hoq


# Extract temp data for sites (from CRU)
ncin <- nc_open("C:/Users/sarah/Documents/PhD/analyses/data/cru_ts4.01.1901.2016.tmp.dat.nc")
tmp_array <- ncvar_get(ncin, "tmp", start = c(1,1,709), count = c(720,360,683)) # reading only data from 1960 onwards (to match oipc)
lon <- ncvar_get(ncin, "lon"); lat <- ncvar_get(ncin, "lat")
nlon <- length(lon); nlat <- length(lat)
time <- ncvar_get(ncin, "time", start = 709, count = 683); time <- as.Date(time, origin = "1900-01-01")
nt <- length(time)
nc_close(ncin); rm(ncin)
# Reshape temp data
lonlat <- expand.grid(lon=lon, lat=lat)
tmp_vec <- as.vector(tmp_array)
tmp_mat <- matrix(tmp_vec, nrow = nlon*nlat, ncol = nt)

a <- sapply(dt_sites$longitude, function(x) which.min(abs(lon-x)))
b <- sapply(dt_sites$latitude, function(x) which.min(abs(lat-x)))
ab <- (b-1)*nlon + a
sisal_tmp <- tmp_mat[ab,]
mean_tmp <- rowMeans(sisal_tmp)

mean_vals <- cbind(dt_sites[,1:6], temp = mean_tmp)

# Convert weighted d18O precip (SMOW) to carbonate (PDB)
d18O_SMOW <- left_join(pts, mean_vals, by = c("site_id","site_name", "latitude","longitude","mineralogy"))
d18O_SMOW <- d18O_SMOW %>% mutate(tmp_K = mean_tmp + 273.15) # convert temp from celcius to K
d18O_SMOW_calc <- d18O_SMOW %>% mutate(carbonate_SMOW = OIPC_d18O + (((16.1*10^3)/tmp_K) - 24.6)) #convert water to calcite (equation from Tremaine et al,2011)
d18O_SMOW_calc <- d18O_SMOW_calc %>% mutate(carbonate_PDB = 0.97001 * carbonate_SMOW - 29.29) # Convert from SMOW to PDB
d18O_SMOW_arag <- d18O_SMOW %>% mutate(carbonate_SMOW = OIPC_d18O + ((18.34*10^3)/tmp_K) - 31.954)
d18O_SMOW_arag <- d18O_SMOW_arag %>% mutate(carbonate_PDB = 0.97001 * carbonate_SMOW - 29.29) # Convert from SMOW to PDB


# Calculate d18O site anomalies

# calcite
dt_calc <- dt %>% filter(mineralogy == "calcite" | arag_corr == "yes")
dt_anom_calc <- left_join(dt_calc, d18O_SMOW_calc[-5])

# aragonite
dt_arag <- dt %>% filter(mineralogy == "aragonite" & arag_corr != "yes")
dt_anom_arag <- left_join(dt_arag, d18O_SMOW_arag[,-5])

# calc anomalies
dt_anom <- rbind(dt_anom_calc, dt_anom_arag)
dt_anom <- dt_anom %>% mutate(d18O_anom = d18O_measurement - carbonate_PDB) %>%
  group_by(site_id, site_name, longitude, latitude, region, time_slice) %>% summarise(d18O_anom = mean(d18O_anom))


# STEP 3: save data
write.csv(dt_anom, "C:/Users/sarah/Documents/PhD/analyses/data/speleothem_boxplot_d18O.csv", row.names = F)


