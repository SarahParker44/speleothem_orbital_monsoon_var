## PCoA/RDA on speleothem d18O correlation values ####

library(RMySQL)
library(dplyr)
library(tidyr)
library(ggplot2)
library(vegan)
library(cowplot)
library(Hmisc)
library(Cairo)


## STEP 1: extract data and filter sites ##

## Connect to SISAL db (in MySQL)

db_user <- 'root'
db_password <- ''
db_name <- 'sisalv2'
db_host <- 'localhost'

mydb <-  dbConnect(MySQL(), user = db_user, password = db_password,
                   dbname = db_name, host = db_host)  # Connect to MySQL db

dbExecute(mydb, 'SET NAMES UTF8;')


## Extract Holocene lower latitude data

Raw_Data <- dbGetQuery(mydb, "SELECT * FROM site JOIN entity USING(site_id) JOIN sample USING(entity_id) JOIN original_chronology USING(sample_id) JOIN d18O USING(sample_id) 
                       WHERE (interp_age <= 12000) AND (latitude BETWEEN -30 AND 45);")  # specify the time range and region that I'm interested in within query


## Filter data to current entity status
Hol_Data <- Raw_Data %>% group_by(site_id, site_name, entity_id, entity_name, longitude, latitude, elevation) %>% filter(entity_status != "superseded") 

## Filter to monsoon regions
Hol_Data$region <- with(Hol_Data, ifelse(latitude >= 15 & latitude <= 35 & longitude >= 50 & longitude <= 98, "ISM",
                                         ifelse(latitude >= 20 & latitude <= 45 & longitude >= 100 & longitude <= 125, "EAM", 
                                                ifelse(latitude >= -10 & latitude <= 0 & longitude >= -80 & longitude <= -70 | latitude >= -30 & latitude <= -10 & longitude >= -60 & longitude <= -30, "SW-SAM",
                                                              ifelse(latitude >= -30 & latitude <= 5 & longitude >= 80 & longitude <= 170, "IAM",
                                                                     ifelse(latitude >= -30 & latitude <= 0 & longitude >= 0 & longitude <= 50, "SAfM", 
                                                                            ifelse(latitude >= 0 & latitude <= 35 & longitude >= -110 & longitude <= -50, "CAM", 
                                                                                   ifelse(latitude >= -10 & latitude <= 0 & longitude >= -70 & longitude <= -30, "NE-SAM", "other"))))))))

Hol_Data <- Hol_Data %>% filter(region != "other")


## Filter to sites covering at least 4000 years

source("C:/Users/sarah/Documents/PhD/monsoon_paper/analyses_figs/PCoA_RDA/d18O_record_length.R")
record_lengths <- get_site_coverage(threshold = 4000)
record_lengths <- record_lengths %>% filter(site_id %in% Hol_Data$site_id)

Hol_Data <- Hol_Data %>% filter(site_id %in% record_lengths$site_id)


## Exclude entities that have variables mineralogies

# Identify calcite or aragonite corrected records
calcite <- Hol_Data %>% group_by(site_name, entity_id, mineralogy, arag_corr) %>% 
  filter(mineralogy == "calcite" | mineralogy == "aragonite" & arag_corr == "yes") %>%
  summarise(n()) #n=97

# Identify aragonite uncorrected records
arag <- Hol_Data %>% group_by(site_name, entity_id, mineralogy, arag_corr) %>% 
  filter(mineralogy == "aragonite" & arag_corr != "yes") %>%
  summarise(n()) #n=37

# Identify mixed uncorrected records
mixed <- Hol_Data %>% group_by(site_name, entity_id, mineralogy, arag_corr) %>% 
  filter(mineralogy == "mixed" & arag_corr != "yes") %>%
  summarise(n()) #n=0

# Unknown mineralogy records
unknown <- Hol_Data %>% group_by(site_name, entity_id, mineralogy, arag_corr) %>% 
  filter(mineralogy == "unknown") %>%
  summarise(n()) #1

# Select entities with variable mineralogy
var_min <- na.omit(full_join(calcite, arag, by = c("site_name","entity_id"))) #0

# Remove varied mineralogy records from analysis (and with unknown mineralogy)
#Hol_Data <- Hol_Data %>% filter(!entity_name %in% var_min$entity_name)
Hol_Data <- Hol_Data %>% filter(!entity_id %in% unknown$entity_id)

## Remove some entities to ensure that sites have a constant mineralogy
x <- Hol_Data %>% filter(site_name == "Dos Anas cave" & mineralogy == "calcite"|
                           site_name == "KNI-51 cave" & mineralogy == "aragonite"|
                           site_name == "Mawmluh cave" & mineralogy == "aragonite" | site_name == "Mawmluh cave" & mineralogy == "unknown"|
                           site_name == "Tamboril cave")
Hol_Data <- Hol_Data %>% filter(!sample_id %in% x$sample_id)

Hol_Data <- Hol_Data %>% filter(!entity_id %in% filter(Hol_Data, site_name == "Lianhua cave" & entity_name == "A1") &
                                  !entity_id %in% filter(Hol_Data, site_name == "Botuver¡" & entity_name  == "BTV21a")) # remove 2 entities to ensure all site records have overlapping entities

# Remove El Condor composite record (contains input errors - 01/06/20)
Hol_Data <- Hol_Data %>% filter(entity_id != 410)

# Extract site data
Hol_sites <- Hol_Data %>% group_by(site_id, longitude, latitude, elevation, region) %>% summarise(min_age = min(interp_age), max_age = max(interp_age), 
                                                                                                  length = max(interp_age) - min(interp_age))
Hol_sites <- Hol_sites %>% mutate(mid = min_age + length/2)

# export records for fig 1 (map of records)
hol_ent <- Hol_Data %>% group_by(site_id, site_name, entity_id, entity_name, latitude,longitude, elevation, region) %>% summarise(n())
write.csv(hol_ent[,1:8], "C:/Users/sarah/Documents/PhD/analyses/data/PCoA_entities.csv", row.names = F)



## STEP 2: BINNING ##

# Function that bins data:
binning <- function(dt, bins, binhw){
  #Reshape data
  dat_res <- dt[,c("site_id", "d18O_measurement")] %>% group_by(site_id) %>% mutate(number = row_number()) %>% spread(site_id, d18O_measurement)
  dat_res <- dat_res[,-1]
  dat_res <- as.matrix(dat_res)
  
  Age_res <- dt[,c("site_id", "interp_age")] %>% group_by(site_id) %>% mutate(number = row_number()) %>% spread(site_id, interp_age)
  Age_res <- Age_res[,-1]
  Age_res <- as.matrix(Age_res)
  
  #binning
  result <- matrix(ncol = length(Age_res[1,]), nrow = length(bins))
  for (k in 1:length(dat_res[1,])){
    if(length(dat_res[is.na(dat_res[,k]) == F, k]) != 0){
      for (i in 1:length(bins)){
        t <- na.omit(cbind(as.numeric(Age_res[,k]), as.numeric(dat_res[,k])))
        result[i,k] <- mean(t[t[,1] > bins[i] - binhw & t[,1] < bins[i] + binhw, 2])
      }
    }
  }
  BinnedData <- structure(result, row.names = as.character(bins), col.names = colnames(dat_res), class = "matrix")
  return(BinnedData)
}

bin_dat <- binning(dt = Hol_Data, bins = seq(0,12000,500), binhw = 250)



## STEP 3: PAIRWISE CORRELATION OF ALL SITES ##

# Produce correlation matrix
corr_mat <- cor(bin_dat, use = "pairwise.complete.obs")
corr_mat <- corr_mat[ , apply(corr_mat, 2, function(x) !any(is.na(x)))] # remove columns with NA values (so variables/columns = correlation coefficient of sites with reasonably complete Holocene records)


## STEP 4: PCoA ON CORRELATION DATA

corr_dist <- vegdist(corr_mat, "euclidean") # Convert correlation coefficient matrix to distances

mod <- betadisper(corr_dist, Hol_sites$region) # PCoA using vegan package (alternative function below, same package, same results)
#mod <- cmdscale(d = corr_dist, eig = T)

# Get eigenvalues/variance explained
eig <- data.frame(Eigenvalue = mod$eig)
eig$sum <- sum(eig$Eigenvalue)
eig$var_explained <- (eig$Eigenvalue/eig$sum)*100
eig$cumulative_var <- cumsum(eig$var_explained)

# Export as table
x <- data.frame(t(eig[1:5,-2]))
row.names(x)[2:3] <- c("Explained (%)","Cumulative (%)")

# export results as csv file
write.csv(x, "C:/Users/sarah/Documents/PhD/monsoon_paper/svg_figs/pcoa_results.csv")


## STEP 5: DETERMINE PCoA SIGNIFICANCE (BROKEN STICK APPROACH)

# generate broken stick model (bsm) data
n <- length(eig$Eigenvalue)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))

eig$b_stick <- rev(100*bsm$p/n) # bsm % var explained

b_stick_df <- data.frame(PCoA_axis = 1:33, var_explained = eig$var_explained, b_stick = eig$b_stick) # df of PCoA % explained versus bsm
b_stick_df <- gather(b_stick_df, "var","var_explained", 2:3)

# plot
ggplot(data = b_stick_df, aes(x = PCoA_axis, y = var_explained, col = var)) +
  geom_line() + geom_point() +
  ylab("variance explained (%)") + xlab("PCoA axes") + 
  scale_x_continuous(breaks = 1:31) + scale_y_continuous(breaks = seq(0,60,10)) +
  scale_color_manual(values = c("red","black"), labels = c("broken stick model","observed")) +
  theme(legend.position = c(0.7,0.8), legend.title = element_blank(), text = element_text(size = 12), axis.text = element_text(size = 10))
# First 2 axes sig


## STEP 5: PLOT DATA

scrs <- scores(mod, choices = c(1,2)) # PCoA scores for each site
COA_dt <- data.frame(Hol_sites, scrs$sites) # combine site data and scores into a df

COA_dt$PCoA1 <- COA_dt$PCoA1 * -1; # Flip PCoA1 scores so that latitude plots left to right (step is specific to my data)

COA_dt$time <- with(COA_dt, ifelse(length >= 8000, "Hol",
                                   ifelse(length <= 8000 & mid <= 5000, "MH to LH",
                                          ifelse(length <= 8000 & mid >= 8000, "EH to MH", "MH")))) # assign a time period to my PCoA scores (will be plotting as shapes)
                                   

find_hull <- function(df) df[chull(df$PCoA1, df$PCoA2), ]
hulls <- filter(COA_dt) %>%
  group_by(region) %>%
  nest() %>%
  mutate(data = purrr::map(data, find_hull)) %>%
  unnest() # Generate a line that joins up all the outer points in a region, for plotting

COA_dt$time <- factor(COA_dt$time, c("Hol","EH to MH", "MH", "MH to LH")) # logical order for my factor variables
COA_dt$region <- factor(COA_dt$region, c("ISM","EAM","IAM","SW-SAM","NE-SAM","CAM","SAfM"))

palette <- c("#FF5722","#42D4F4","#FFA000","#F032E6","#8CBD00","#4363d8","#a9a9a9") # custom palette for regions

# Plot PCoA results
p1 <- ggplot() + geom_point(data = COA_dt, aes(x = PCoA1, y = PCoA2, col = region, shape = time), size = 1.5) +
  scale_color_manual(values = palette, guide = F) +
  theme_bw() + xlab("PCoA 1") +
  ylab("PCoA 2") +
  geom_polygon(data = hulls, aes(x = PCoA1, y = PCoA2, col = region), fill = NA, show.legend = F, size = 0.2) +
  theme(legend.title = element_blank(), legend.position = c(0.8,0.25),
        axis.text = element_text(size = 8), axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.8, "lines"))

# Create legend as individual plot
px <- ggplot() + geom_point(data = COA_dt, aes(x = PCoA1, y = PCoA2, col = region), size = 1.5) +
  scale_color_manual(values = palette, guide = guide_legend(title.position = "top", nrow = 1)) +
  theme(legend.box = "horizontal", legend.direction = "horizontal", legend.title = element_blank(), legend.text = element_text(size = 8), legend.key.width = unit(1, "cm"),
        legend.margin = margin(t = 0, r = 5, b = 0, l = 5, unit = "cm"), legend.key = element_rect(fill = NA)) +
  scale_shape(guide = guide_legend(title.position = "top")) 

grobs <- ggplotGrob(px)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]



## STEP 7: RDA ON PCoA SCORES
# Using PCoA 1 and 2 as response variables, as these are the sig axes

# Explanatory variables (latitude, longitude, elevation)
y <- COA_dt[,c(2:3)]

# centre variables on their mean and standardise
y_scale <- scale(y, center = TRUE, scale = T)

# Response variables (PCoA 1 and 2)
x <- COA_dt[,10:11]

# centre variables on their mean and standardise
x_scale <- scale(x, center = T, scale = T)

# RDA (triplot in ggplot amended from https://oliviarata.wordpress.com/2014/09/06/rda-in-ggplot2/)
pcoa.rda <- rda(x_scale ~ . , data = as.data.frame(y_scale), scale = T) # vegan function
rda_sum <- summary(pcoa.rda) # get results from RDA function

scor <- scores(pcoa.rda, 
               display = c("sp","wa","bp"),
               scaling = 2)

# site scores
site_scor <- data.frame(scor$sites)

p2 <- ggplot() +
  geom_point(site_scor, mapping = aes(x = RDA1, y = RDA2, col = COA_dt$region), size = 1.5) +
  scale_colour_manual(values = palette) + theme_bw() +
  geom_vline(mapping = aes(xintercept = 0), lty = 2, size = 0.2) + geom_hline(mapping = aes(yintercept = 0), lty = 2, size = 0.2)

# species
species_scor <- data.frame(scor$species)

p2 <- p2 +
  geom_text(data = species_scor, aes(x = c(2.1,1.6), y= c(0.31,-0.35), label = rownames(species_scor)), col = "red", size = 8/(14/5)) +
  geom_segment(data = species_scor, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.25, "cm")), col = "red", size = 0.2)

# biplot arrows
y_arrows <- data.frame(scor$biplot)
basplot <- plot(pcoa.rda, scaling = 2)
mult <- attributes(basplot$biplot)$arrow.mul

p2 <- p2 +
  geom_segment(data = y_arrows, aes(x = 0, y = 0, xend = mult*RDA1, yend = mult*RDA2),
               arrow = arrow(length = unit(0.25, "cm")), size = 0.2) +
  geom_text(data = y_arrows, aes(x = c(5,4.6), y = c(-3.5,4.5),
                                 label = rownames(y_arrows)), size = 8/(15/5)) +
  theme(legend.position = "none", axis.text = element_text(size = 8), axis.title = element_text(size = 8))


# Table for eig, var
a <- data.frame(rbind(rda_sum$cont$importance))
#write.csv(a, "../R_programming/Results/May_2019/PCoA_RDA/eig_var_RDA.csv")

# table of response var loadings
b <- data.frame(rda_sum$biplot[,1:2])

# Get p values
z <- data.frame(rda_sum$sites[,1:2], y) #df of RDA scores, lon, lat, elev for each site
RDA1_pval <- rcorr(as.matrix(z))$P[3:4,1] # P vals of RDA1 vs lon, lat, elev
RDA2_pval <- rcorr(as.matrix(z))$P[3:4,2] # P vals of RDA2 vs lon, lat, elev
c <- data.frame(b$RDA1, RDA1_pval, b$RDA2, RDA2_pval)
colnames(c) <- c("RDA1", "RDA1 p-value", "RDA2", "RDA2 p-value")
#write.csv(c, "../R_programming/Results/May_2019/PCoA_RDA/constr_var_RDA.csv")


### Multiplot of PCoA, RDA and legend
p <- plot_grid(p1, p2, labels = c("(a)","(b)"), label_size = 10, hjust = -0.2, label_fontface = "plain")
p <- plot_grid(p, legend, nrow = 2, rel_heights = c(1,0.1), align = "v", axis = "b", rel_widths = c(1,10))
print(p)

# export
#Cairo(file = "C:/Users/sarah/Documents/PhD/monsoon_paper/svg_figs/test_pcoa.pdf", type = "pdf", units = "cm", width = 16, height = 8)
pdf("C:/Users/sarah/Documents/PhD/monsoon_paper/svg_figs/test_pcoa.pdf", width = 16/2.54, height = 8/2.54)
print(p)
dev.off()
