## Plot map of sites used in my paper ##################

library(rgdal)
library(ggplot2)
library(dplyr)
library(ggnewscale)

setwd("C:/Users/sarah/Documents/PhD/analyses/data/")


#  load world map data
wmap <- readOGR(dsn = "ne_110m_land", layer = "ne_110m_land")
wmap@data$id <- rownames(wmap@data)
worldMap <- fortify(wmap)
wmap_DF <- merge(worldMap, wmap@data, by = "id")


# load WOKAM (carbonate bedrock) data
wokam <- readOGR(dsn = "WHYMAP_WOKAM/shp", layer = "whymap_karst__v1_poly")
wokam@data$id <- rownames(wokam@data)
wokamMap <- fortify(wokam)
wokam_DF <- merge(wokamMap, wokam@data, by = "id")


# Import site data
PCoA_sites <- read.csv("PCoA_entities.csv")
boxplot_sites <- read.csv("boxplot_sites.csv")
composite_sites <- read.csv("composite_entities.csv")
LIG_sites <- read.csv("LIG_entities.csv")

length(unique(c(PCoA_sites$entity_id, boxplot_sites$entity_id, composite_sites$entity_id, LIG_sites$entity_id)))


# edit data and merge into 1 df
boxplots <- unique(boxplot_sites[,-c(5:7)]); boxplots$boxplot <- T
comps <- unique(composite_sites[,-c(3:4,7:8)]); comps$composite <- T
PCoA <- unique(PCoA_sites[,-c(3,4,7,8)]); PCoA$PCoA <- T

dat <- full_join(boxplots, comps, by = c("site_id","site_name","latitude","longitude")) # merge
dat <- full_join(dat, PCoA, by = c("site_id","site_name","latitude","longitude"))


# create new column of which analyses a site has been used in
dat[which(dat$PCoA == T),"analyses"] <- "used in PCoA/RDA"
dat[is.na(dat$analyses),"analyses"] <- "used in other analyses"


# region limits for monsoon regions (used in multiple regression)
reg_mult <- data.frame(x1 = c(50,100,95,-60,-115,10), x2 = c(95,125,135,-30,-58,40),
                       y1 = c(11,20,-24,-10,10,-30), y2 = c(32,39,5,0,33,-17),
                       type = "monsoon regions") # monsoon regions

# regiona limits for source areas
reg_source <- data.frame(x1 = c(50,75,85,-80,-38,-35,20), x2 = c(75,120,140,-40,-20,-13,55),
                         y1 = c(0,0,-18,10,-20,-17,-35), y2 = c(25,22,0,20,-10,-7,-10),
                         type = "source regions") # source regions

reg <- rbind(reg_mult, reg_source)

# region names to add as annotations to map
reg_names <- data.frame(name = c("CAM","SW-SAM","NE-SAM","SAfM","ISM","EAM","IAM"),
                        x = c(-60,-22,-32,40,60,132,90),
                        y = c(36,-23,4,-33,7,20,-22))

# reorder variables
dat$analyses <- factor(dat$analyses, c("used in PCoA/RDA","used in other analyses"))

# plot map
map <- ggplot() +
  
  # background of world map and carbonate bedrock
  geom_polygon(data = wmap_DF, aes(x = long, y = lat, group = group), fill = "#F1F1F1") + # fill of continents
  geom_polygon(data = wokam_DF, aes(x = long, y = lat, group = group, fill = "Carbonate bedrock")) + # fill of wokam data
  scale_fill_manual(values = "#FDAE6B")  + # set fill of wokam data
  geom_path(data = wmap_DF, aes(x = long, y = lat, group = group), col = "#585757", size = 0.1) + # outline of continents
  
  # plot sites
  new_scale_fill() + # reset fill (2 fill settings here = wokam and sites)
  geom_point(data = dat, aes(x = longitude, y = latitude, fill = analyses), shape = 21, col = "black", size =1.52, stroke = 0.5) + #site data as points
  scale_fill_manual(values = c("#B3589A", "#9BBF85")) +
  
  # plot region limits, with labels
  #geom_rect(data = d, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = NA, col = "black", size = 0.2) + # speleothem region limits
  geom_rect(data = reg, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, col = type), fill = NA, size = 0.2) +
  scale_color_manual(values = c("black","red")) + 
  geom_polygon(mapping = aes(x = c(-80,-64,-64,-40,-40,-68,-68,-80), y = c(0,0,-10,-10,-30,-30,-10,-10)), fill = NA, col = "black", size = 0.2) +
  geom_text(data = reg_names, aes(x = x, y = y, label = name), size = 8/(14/5)) + # region labels
  
  # edit overall appearance of graph
  theme_bw() + # plain theme
  scale_x_continuous(sec.axis = dup_axis(), expand = c(0,0)) +
  scale_y_continuous(sec.axis = dup_axis(), expand = c(0,0)) + # axis ticks on all sides of map
  coord_fixed(ylim = c(-48,45), xlim = c(-170,160)) + # Cut out Pacific and high latitudes
  theme(legend.title = element_blank(), # no legend titles
        axis.title = element_blank(), #no axis titles
        panel.background = element_rect(fill = "#F6FFFF"), panel.grid = element_blank(), # set background colour, remove grid lines
        axis.text = element_blank(), # no axis labels
        axis.ticks.length = unit(-0.1,"cm"), # ticks inside the plot
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA), # remove legend background colours
        legend.text = element_text(size = 8), legend.key.size = unit(0.2, "cm"), # change size of legend
        legend.spacing.y = unit(-0.1, "cm"), # reduce space between 2 keys
        legend.position = c(0.13,0.25) # inset legend
  ) 

print(map)

# save fig  
pdf("sites_map.pdf", width = 15.9/2.54, height = 5/2.54, useDingbats = F)
print(map)
dev.off()
