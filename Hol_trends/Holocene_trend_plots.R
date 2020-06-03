#### Make plots of Holocene d18Ospel, d18Oprecip and summer insolation ####

setwd("C:/Users/sarah/Documents/PhD/analyses/")

# load necessary packages

library(palinsol)
library(ggplot2)
library(ggthemes)
library(cowplot)


# Load d18Ospel composite data, calculated in "region_spel_composites.R"

Sisal_comp <- read.csv("data/Sisal_Comp_dt.csv")

# Calculate summer insolation

# Function that calculates mean insolation for a specified lat and sequence of months
get_insol <- function(Lat, Months){
  mid_month_days_0ka <- seq(15.5, 345.5, by=30)
  tt_present = 0.0
  orbit_present <- astro(tt_present, ber78, degree = FALSE) # present-day atronomical parameters (using Berger, 1978)
  mid_month_tsl_0ka <- day2l(orbit_present, mid_month_days_0ka) #convert days to true solar longitude
  time_BP <- seq(-12000,0,by=100) # Holocene time intervals
  orbital_params <- data.frame(time_BP, t(sapply(time_BP, function(tt) astro(tt, ber78, degree=FALSE)))) # Holocene astro params
  insol_month <- matrix(0, nrow=length(time_BP), ncol=12)
  for (month in seq(1:12)) {
    tsl <- mid_month_tsl_0ka[month] # month solar longitude
    insol_month[,month] <- as.numeric(Insol(orbital_params, long=tsl, lat=Lat*pi/180, S0=1365)) # calc insolation (convert lat to radians) 
  }
  Months_insol <- data.frame(insol_month[,Months]) # filter to specified months
  Months_insol2 <- data.frame(AGE = time_BP*-1, insol = apply(Months_insol, 1, mean)) # mean insolation across months
}

# Calculate Holocene summer insolation for NH and SH
NH_insol <- get_insol(Lat = 30, Month = c(5,6,7,8,9)) # mean MJJAS insolation (NH summer)
SH_insol <- get_insol(Lat = -20, Month = c(11,12,1,2,3)) # mean NDJFM insolation (SH summer)


# Load GISS d18Oprecip trends
GISS <- read.csv("data/Bootstrap_GISS_wiso.csv")


# NH plots

# NH insolation
p1 <- ggplot(data = NH_insol, aes(x = AGE, y = insol)) + 
  annotate("rect", xmin = c(0,4000,8000), xmax = c(2000,6000,10000), ymin = min(NH_insol$insol), ymax = max(NH_insol$insol), fill = "#F1F1F1") +
  geom_line() +
  scale_x_continuous(position = "top", breaks = seq(0,12000,2000)) +
  scale_y_continuous(position = "right") +
  geom_rangeframe(sides = "r") + theme_tufte() +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size = 8), axis.title = element_text(size = 8)) +
  ylab("30°N insolation")

# ISM speleothem composite
p2 <- ggplot(data = filter(Sisal_comp, region == "ISM"), aes(x = age, ymin = conf5, ymax = conf95, y = seq(-2.21,6.3,8.51/119))) +
  annotate("rect", xmin = c(0,4000,8000), xmax = c(2000,6000,10000), ymin = 6, ymax = -2.25, fill = "#F1F1F1") +
  geom_ribbon(fill = "#ffa58a") + 
  geom_line(data = filter(Sisal_comp, region == "ISM"), aes(x = age, y = locfit), col = "#FF5722") +
  geom_rangeframe() + theme_tufte() + scale_y_reverse(breaks = seq(-2,6,2)) +
  ylab(expression(paste(delta^18,"O Z-Score", sep = ""))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8), axis.line.x = element_line(colour = "white"), axis.line.y = element_line(colour = "#FF5722"), axis.ticks.y = element_line(colour = "#FF5722"))

# ISM GISS d18O
p3 <- ggplot() +
  annotate("rect", xmin = c(0,4,8), xmax = c(2,6,10), ymin = -2.6, ymax = 0.1, fill = "#F1F1F1") +
  geom_ribbon(data = filter(GISS, region == "ISM"), aes(x = age, ymin = conf5, ymax = conf95), fill = "#ffa58a", alpha = 0.5) +
  geom_line(data = filter(GISS, region == "ISM"), aes(x = age, y = locfit), col = "#FF5722", linetype = "longdash") +
  geom_rangeframe(sides = "rb", aes(y = c(0.1,-2.6), x = c(0,12))) + theme_tufte() + expand_limits(x = c(0,12)) +
  scale_x_continuous(breaks = seq(0,12,2)) + scale_y_reverse(position = "right", breaks = seq(-2.4,0,0.8)) + 
  ylab(expression(paste(delta^18,"O of ppt"))) +
  theme(axis.title.x = element_blank(), axis.text = element_text(size = 8),
        axis.title = element_text(size = 8))#, axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())

# EAM speleothem composite
p4 <- ggplot(data = filter(Sisal_comp, region == "EAM"), aes(x = age, ymin = conf5, ymax = conf95, y = seq(-2.21,6.3,8.51/119))) +
  annotate("rect", xmin = c(0,4000,8000), xmax = c(2000,6000,10000), ymin = 7, ymax = -2, fill = "#F1F1F1") +
  geom_ribbon(fill = "#a3eafa") + 
  geom_line(data = filter(Sisal_comp, region == "EAM"), aes(x = age, y = locfit), col = "#42D4F4") +
  geom_rangeframe() + theme_tufte() + scale_y_reverse(breaks = seq(-2,6,2)) +
  ylab(expression(paste(delta^18,"O Z-Score", sep = ""))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8), axis.line.x = element_line(colour = "white"), axis.line.y = element_line(colour = "#42D4F4"), axis.ticks.y = element_line(colour = "#42D4F4"))

# EAM GISS d18O
p5 <- ggplot() +
  annotate("rect", xmin = c(0,4,8), xmax = c(2,6,10), ymin = -1, ymax = 0, fill = "#F1F1F1") +
  geom_ribbon(data = filter(GISS, region == "EAM"), aes(x = age, ymin = conf5, ymax = conf95), fill = "#a3eafa", alpha = 0.5) +
  geom_line(data = filter(GISS, region == "EAM"), aes(x = age, y = locfit), col = "#42D4F4", linetype = "longdash") +
  geom_rangeframe(sides = "rb", aes(y = c(-1,0))) + theme_tufte() +
  scale_x_continuous(breaks = seq(0,12,2)) + scale_y_reverse(position = "right", breaks = seq(-1.2,0,0.3)) +
  ylab(expression(paste(delta^18,"O of ppt"))) + expand_limits(x = c(0,12)) +
  theme(axis.title.x = element_blank(), axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())


# SH plots

# SH insolation
p6 <- ggplot(data = SH_insol, aes(x = AGE, y = insol)) + 
  annotate("rect", xmin = c(0,4000,8000), xmax = c(2000,6000,10000), ymin = min(SH_insol$insol), ymax = max(SH_insol$insol), fill = "#F1F1F1") +
  geom_line() +
  scale_x_continuous(position = "top", breaks = seq(0,12000,2000)) +
  scale_y_continuous(position = "right", breaks = c(455,455,465)) +
  geom_rangeframe(sides = "r") + theme_tufte() +
  theme(panel.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size = 8), axis.title = element_text(size = 8)) +
  ylab("20°S insolation")

# IAM speleothem composite
p7 <- ggplot(data = filter(Sisal_comp, region == "IAM"), aes(x = age, ymin = conf5, ymax = conf95, y = seq(-4,11.3,15.3/119))) +
  annotate("rect", xmin = c(0,4000,8000), xmax = c(2000,6000,10000), ymin = 11.3, ymax = -4, fill = "#F1F1F1") +
  geom_ribbon(fill = "#ffdda3") + 
  geom_line(data = filter(Sisal_comp, region == "IAM"), aes(x = age, y = locfit), col = "#FFA000") +
  geom_rangeframe() + theme_tufte() + scale_y_reverse(breaks = seq(-4,8,4)) + scale_x_continuous(breaks = seq(0,12000,2000), labels = seq(0,12,2)) +
  ylab(expression(paste(delta^18,"O Z-Score", sep = ""))) +
  theme(axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8), axis.title.x = element_blank(), 
        axis.line.x = element_line(colour = "white"), axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.line.y = element_line(colour = "#FFA000"), axis.ticks.y = element_line(colour = "#FFA000"), plot.background = element_blank())

# IAM GISS d18O
p8 <- ggplot() +
  annotate("rect", xmin = c(0,4,8), xmax = c(2,6,10), ymin = min(filter(GISS, region == "IAM")$conf5), ymax = max(filter(GISS, region == "IAM")$conf95), fill = "#F1F1F1") +
  geom_line(data = filter(GISS, region == "IAM"), aes(x = age, y = locfit), col = "#FFA000", linetype = "longdash") + 
  geom_ribbon(data = filter(GISS, region == "IAM"), aes(x = age, ymin = conf5, ymax = conf95), alpha = 0.5, fill = "#ffdda3") +
  geom_rangeframe(sides = "rb", mapping = aes(y = c(0.85, -0.36))) + theme_tufte() +
  scale_x_continuous(breaks = seq(0,12,2)) + scale_y_reverse(position = "right", breaks = c(-0.3,0,0.3,0.6)) +
  ylab(expression(paste(delta^18,"O of ppt"))) + expand_limits(x = c(0,12)) +
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 8), axis.title.x = element_blank())


# SW-SAM speleothem composite
p9 <- ggplot() +
  annotate("rect", xmin = c(0,4000,8000), xmax = c(2000,6000,10000), ymin = 3.2, ymax = -2.4, fill = "#F1F1F1") +
  geom_ribbon(data = filter(Sisal_comp, region == "SW-SAM"), aes(x = age, ymin = conf5, ymax = conf95), fill = "#F5B9F1") +
  geom_line(data = filter(Sisal_comp, region == "SW-SAM"), aes(x = age, y = locfit), col = "#F032E6") +
  geom_rangeframe(mapping = aes(y = c(-2.4,3.2))) + theme_tufte() + scale_y_reverse(breaks = seq(-2,4,2)) +
  ylab(expression(paste(delta^18,"O Z-Score", sep = ""))) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size = 8),
        axis.title.y = element_text(size = 8), axis.line.x = element_line(colour = "white"), axis.line.y = element_line(colour = "#F032E6"), axis.ticks.y = element_line(colour = "#F032E6"))

# SW-SAM GISS d18O
p10 <- ggplot() +
  annotate("rect", xmin = c(0,4,8), xmax = c(2,6,10), ymin = min(filter(GISS, region == "SW-SAM")$conf5), ymax = max(filter(GISS, region == "SW-SAM")$conf95), fill = "#F1F1F1") +
  geom_line(data = filter(GISS, region == "SW-SAM"), aes(x = age, y = locfit), col = "#F032E6", linetype = "longdash") + 
  geom_ribbon(data = filter(GISS, region == "SW-SAM"), aes(x = age, ymin = conf5, ymax = conf95), alpha = 0.5, fill = "#F5B9F1") +
  geom_rangeframe(sides = "rb", mapping = aes(y = c(-0.2,0.7), x = c(0,12))) + theme_tufte() +
  scale_x_continuous(breaks = seq(0,12,2)) + scale_y_reverse(position = "right", breaks = seq(-0.2,0.6,0.2)) +
  ylab(expression(paste(delta^18,"O of ppt"))) + expand_limits(x = c(0,12)) +
  theme(axis.text = element_text(size = 8),
        axis.title.y = element_text(size = 8), axis.title.x = element_blank())


# legend
px <- ggplot(data = GISS, aes(x = age)) + 
  geom_line(mapping = aes(y = locfit, lty = "speleothem")) +
  geom_line(mapping = aes(y = locfit, lty = "model")) + 
  labs(lty = NULL) + theme_bw() + theme(legend.position = "bottom") +
  scale_linetype_manual(values = c("longdash","solid"), breaks = c("speleothem","model")) 
grobs <- ggplotGrob(px)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


# export as multi plot

pdf("C:/Users/sarah/Documents/PhD/monsoon_paper/svg_figs/Holocene_comps.pdf", width = 7, height = 10)
plot_grid(p1,p6,p4,p7,p5,p8,p2,p9,p3,p10, ncol = 2, align = "hv", rel_heights = c(2,3,3,3,3))
dev.off()
