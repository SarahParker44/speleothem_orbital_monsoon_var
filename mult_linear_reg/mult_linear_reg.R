
## Multiple regression analysis

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

setwd("C:/Users/sarah/Documents/PhD/analyses/data/multiple_regression/")

# Load dat
region_d18O <- read.csv("region_d18O.csv")
region_ppt <- read.csv("region_summer_ppt.csv")
source_temp <- read.csv("source_summer_temp.csv")
wind_dir <- read.csv("source_summer_wdir.csv")
recycling <- read.csv("region_ppt_recycling.csv")

# Join data into 1 df
dat <- full_join(region_d18O, region_ppt, by = c("region", "t_slice"))
dat <- full_join(dat, source_temp, by = c("region", "t_slice"))
dat <- full_join(dat, wind_dir, by = c("region","t_slice"))
dat <- full_join(dat, recycling, by = c("region","t_slice"))

# multiple regression
fit <- lm(region_d18O ~ region_ppt + source_temp + wind_dir + ppt_rec, data = dat) # can use lm() or glm()
summary(fit)


# extract data for partial residual plots
partial.res <- residuals(fit, "partial") # partial residuals for each explanatory variable
x <- model.matrix(fit) # Partial predicted Y values for each variable
colnames(partial.res) <- paste(colnames(partial.res), "_res", sep = "")

resid_dat <- cbind(dat[,1:2], data.frame(cbind(x, partial.res)))


# partial residual plots
palette <- c("#FF5722","#42D4F4","#FFA000","#F032E6","#8CBD00","#4363d8","#a9a9a9") # custom palette for regions
resid_dat$region <- factor(resid_dat$region, c("ISM","EAM","IAM","SW-SAM","NE-SAM","CAM","SAfM"))

p <- ggplot(data = resid_dat, aes(x = region_ppt, y = region_ppt_res, col = region, #shape = as.factor(t_slice)
                                  )) + 
  geom_point(size = 0.5) + ylab(NULL) +
  geom_smooth(method = "lm", fill = NA, size = 0.5) +
  xlab("precipitation (mm/d)") + theme_bw() + 
  scale_shape_manual(values = 1:8) + scale_colour_manual(values = palette) +
  geom_abline(slope = lm(region_ppt_res ~ region_ppt, data = resid_dat)$coefficients[2], intercept = lm(region_ppt_res ~ region_ppt, data = resid_dat)$coefficients[1]) +
  theme(legend.position = "none", text = element_text(size = 8)) +
  annotate(geom = "text", x = -0.65, y = 1.35, label = "(a)", size = 10/(14/5))
p2 <- ggplot(data = resid_dat, aes(x = source_temp, y = source_temp_res, col = region#, shape = as.factor(t_slice)
                                   )) + 
  geom_smooth(method = "lm", fill = NA, size = 0.5) +
  geom_point(size = 0.5) + ylab(NULL) + xlab("source temperature (°C)") + theme_bw() + 
  scale_shape_manual(values = 1:8) + scale_colour_manual(values = palette) +
  geom_abline(slope = lm(source_temp_res ~ source_temp, data = resid_dat)$coefficients[2], intercept = lm(source_temp_res ~ source_temp, data = resid_dat)$coefficients[1])+
  theme(legend.position = "none", text = element_text(size = 8)) +
  annotate(geom = "text", x = -0.75, y = 0.9, label = "(b)", size = 10/(14/5))
p3 <- ggplot(data = resid_dat, aes(x = wind_dir, y = wind_dir_res, col = region#, shape = as.factor(t_slice)
                                   )) + 
  geom_smooth(method = "lm", fill = NA, size = 0.5) +
  geom_point(size = 0.5) + ylab(NULL) + xlab("wind direction (°)") + theme_bw() + 
  scale_shape_manual(values = 1:8) + scale_colour_manual(values = palette) +
  geom_abline(slope = lm(wind_dir_res ~ wind_dir, data = resid_dat)$coefficients[2], intercept = lm(wind_dir_res ~ wind_dir, data = resid_dat)$coefficients[1]) +
  theme(legend.position = "none", text = element_text(size = 8)) +
  annotate(geom = "text", x = -17, y = 1.25, label = "(c)", size = 10/(14/5))
p4 <- ggplot(data = resid_dat, aes(x = ppt_rec, y = ppt_rec_res, col = region#, shape = as.factor(t_slice)
                                   )) + 
  geom_smooth(method = "lm", fill = NA, size = 0.5) +
  geom_point(size = 0.5) + ylab(NULL) + xlab("RI") + theme_bw() + 
  scale_shape_manual(values = 1:7) + scale_colour_manual(values = palette) +
  geom_abline(slope = lm(ppt_rec_res ~ ppt_rec, data = resid_dat)$coefficients[2], intercept = lm(ppt_rec_res ~ ppt_rec, data = resid_dat)$coefficients[1]) +
  theme(legend.position = "none", text = element_text(size = 8)) +
  annotate(geom = "text", x = -0.033, y = 1.1, label = "(d)", size = 10/(14/5))


# legend
px <- ggplot(data = resid_dat, aes(x = region_ppt, y = region_ppt_res, col = region#, shape = as.factor(t_slice)
                                   )) + 
  geom_point() + ylab("Component + resiudal(wiso)") + xlab("pre") + theme_bw() + 
  scale_shape_manual(values = 1:8) + scale_colour_manual(values = palette) +
  geom_abline(slope = lm(region_ppt_res ~ region_ppt, data = resid_dat)$coefficients[2], intercept = lm(region_ppt_res ~ region_ppt, data = resid_dat)$coefficients[1]) +
  theme(legend.box = "horizontal", text = element_text(size = 8)) + labs(colour = NULL)

grobs <- ggplotGrob(px)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]


pdf(file = "C:/Users/sarah/Documents/PhD/monsoon_paper/figs/mlr_fig.pdf", width = 15/2.54, height = 12/2.54)
grid.arrange(
  p,p2,p3,p4,legend,
  widths = c(3,3,1),
  layout_matrix = rbind(c(1,2,NA),
                        c(1,2,5),
                        c(3,4,5),
                        c(3,4,NA)),
  left = textGrob(expression(paste("f(", delta^{18}, "O" ["precip"], ")")), rot = 90, gp = gpar(fontsize = 8, fontface  = "bold")))
dev.off()

