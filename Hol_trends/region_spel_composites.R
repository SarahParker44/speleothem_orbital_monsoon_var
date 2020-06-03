
### Generate Holocene regional composites of speleothem d18O #####################

setwd("C:/Users/sarah/Documents/PhD/analyses/")

# Load necessary packages

library(RMySQL)
library(dplyr)
library(tidyr)
library(locfit)


# Connect to database

db_user <- 'root'
db_password <- ''
db_name <- 'sisalv2'
db_host <- 'localhost'

mydb <-  dbConnect(MySQL(), user = db_user, password = db_password,
                   dbname = db_name, host = db_host)  # Connect to MySQL db

dbExecute(mydb, 'SET NAMES UTF8;')


# Load Hol data between -30 and 45 degress lat
Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING(site_id) JOIN sample USING(entity_id) JOIN original_chronology USING(sample_id) JOIN d18O USING(sample_id) 
                       WHERE (interp_age BETWEEN 0 AND 12000) AND (latitude BETWEEN -30 AND 45);")  


# Filter to monsoon regions

region_data <- Raw_Data %>% mutate(region = ifelse(longitude >= 50 & longitude <= 100 & latitude >= 15 & latitude <= 35, "ISM",
                                                  ifelse(longitude >= 100 & longitude <= 125 & latitude >= 25 & latitude <= 45, "EAM",
                                                         ifelse(longitude >= 70 & longitude <= 170 & latitude >= -30 & latitude <= 5, "IAM",
                                                                ifelse(longitude >= -80 & longitude <= -70 & latitude >= -10 & latitude <= 0 | longitude >= -70 & longitude <= -30 & latitude >= -30 & latitude <= -10, "SW-SAM", 
                                                                       "other")))))

region_data <- region_data %>% filter(region != "other")

#remove superseded entities
region_data <- region_data %>% filter(entity_status != "superseded")

# remove unknown mineralogy entities
region_data <- region_data %>% filter(mineralogy != "unknown")

# remove uncorrected mixed mineralogy entities
if (length(which(region_data$mineralogy == "mixed" & region_data$arag_corr == "no", arr.ind = T) != 0)){
  region_data <- region_data[,-which(region_data$mineralogy == "mixed" & region_data$arag_corr == "no", arr.ind = T)]
}

region_data <- region_data %>% filter(entity_id != 410) #remove El Condor composite (contains some input errors: 03/06/20)

# Choose base period for Z-Score calculations 

basebeg <- 3000
baseend <- 7000

# Exclude records that aren't suitable with chosen base period (removes sites that don't have at least 20 data points within base period)

My_Data <- region_data %>%
  group_by(entity_id) %>%
  filter(length(d18O_measurement[interp_age >= basebeg & interp_age <= baseend])> 20) # Selects sites with more than 20 points within base period


# Exclude entities that have variables mineralogies

# Identify calcite or aragonite corrected records
calcite <- My_Data %>% group_by(site_name, entity_name, mineralogy, arag_corr) %>% 
  summarise(n()) %>%
  filter(mineralogy == "calcite" | mineralogy == "aragonite" & arag_corr == "yes")

# Identify mixed and aragonite uncorrected records
arag_mixed <- My_Data %>% group_by(site_name, entity_name, mineralogy, arag_corr) %>% 
  summarise(n()) %>%
  filter(mineralogy == "aragonite" & arag_corr != "yes" | mineralogy == "mixed" & arag_corr != "yes")

# Select entities with variable mineralogy
var_min <- na.omit(full_join(calcite, arag_mixed, by = c("site_name","entity_name"))) #n=0

# Remove varied mineralogy records from analysis
#My_Data <- My_Data %>% filter(!entity_name %in% var_min$entity_name)

# remove (uncorrected) mixed mineralogies
My_Data %>% filter(mineralogy == "mixed" & arag_corr != "yes") #n=0

# summarise sites and entities in this analysis
dt_sites <- My_Data %>% group_by(site_name, site_id, longitude, latitude) %>% summarise(n())
dt_entities <- My_Data %>% group_by(site_id, site_name, entity_id, entity_name, latitude, longitude,elevation, region) %>% summarise(n())

# Extract entities for fig 1
write.csv(dt_entities[,-9], "data/composite_entities.csv", row.names = F)


# Transform speleothem d18O to Z-Scores

# Z-Score = (d18O - mean)/sd
# Where mean and sd are calculated for d18O values within pre-defined base period

ZScore <- function(d18O, interp_age, basebeg, baseend){
  mean <- mean(d18O[interp_age >= basebeg & interp_age <= baseend])
  sd <- sd(d18O[interp_age >= basebeg & interp_age <= baseend])
  Z_Score <- (d18O-mean)/sd
  
  return(Z_Score)
}


# calculate Z-Scores for each site
Trans_Data <- My_Data %>% 
  group_by(site_name) %>%
  mutate(Z_Score = ZScore(d18O_measurement, interp_age, basebeg, baseend))


# Binning, loess and bootstrapping using paleofire pacakage #########

bins <- seq(0,12000,100) # start, end bin size
nboot <- 1000 # number of iterations for bootstrapping
centres <- seq(50,11950,100) # bin centres
hw <- 1500 # half window length for smoothing

# Function for reshaping data, binning, loess fitting and bootstrap resampling by site
# Reshape = Mutating a sequence of numbers for each site, then spreading data out so that cols = site_id, rows = number sequences and fill = Z-Score
# [,-1] at the end removes the number sequences (no longer required))

Comp_get <- function(dt){
  
  #Reshape data
  Trans <- dt[,c("site_id", "Z_Score")] %>% group_by(site_id) %>% mutate(number = row_number()) %>% spread(site_id, Z_Score)
  Trans <- Trans[,-1]
  Trans <- as.matrix(Trans)
  
  Age <- dt[,c("site_id", "interp_age")] %>% group_by(site_id) %>% mutate(number = row_number()) %>% spread(site_id, interp_age)
  Age <- Age[,-1]
  Age <- as.matrix(Age)
  
  #binning
  result <- matrix(ncol = length(Age[1,]), nrow = length(bins)-1)
  for (k in 1:length(Trans[1,])){ # for each site
    c1 <- cut(Age[,k], breaks = bins) 
    tmean <- tapply(Trans[,k], c1, na.omit(mean)) # calc mean for each bin
    result[,k] <- c(as.numeric(tmean))
  }
  
  #bootstrapping 
  mboot <- matrix(nrow = length(bins)-1, ncol = nboot)
  set.seed(1)
  
  for (j in 1:nboot){
    ne <- sample(seq(1,ncol(result),1), ncol(result), replace = T) # randomly sample columns (cols = sites)
    y <- as.vector(as.matrix(result)[,ne])
    x <- as.vector(rep(centres, length(ne)))
    
    locboot <- locfit(y ~ lp(x, deg = 1, h = hw), maxk = 2000, family = "grgauss")
    predboot <- predict(locboot, newdata = centres, se.fit = T) # loess smoothed fit through randomly sampled data
    mboot[, j] <- predboot$fit # save output
  }
  
  # Confidence intervals from bootstrapping
  bootci <- t(apply(mboot, 1, quantile, probs = c(0.05,0.95), na.rm = T))
  
  rm(x, y)
  
  # Loess regression of all sites
  y <- c(result)
  x <- rep(centres, ncol(result))
  dat <- as.data.frame(cbind(x,y))
  dat <- na.omit(dat[order(x),])
  x <- as.vector(dat$x)
  y <- as.vector(dat$y)
  
  locbootA <- locfit(y ~ lp(x, deg = 1, h = hw), maxk = 2000, family = "grgauss")
  predbootA <- predict(locbootA, newdata = centres, se.fit = T)
  LocfitA <- predbootA$fit
  
  # Create output data frame
  locfit_list <- as.data.frame(cbind(centres, LocfitA, bootci))
  colnames(locfit_list) <- c("age","locfit","conf5","conf95")
  
  return(locfit_list)
}

# Composite for each region
COMP_dt <- Trans_Data %>%
  group_by(region) %>%
  nest() %>%
  mutate(data = purrr::map(data, Comp_get)) %>%
  unnest()


# Export data
write.csv(COMP_dt, "data/Sisal_Comp_dt.csv", row.names = F)

