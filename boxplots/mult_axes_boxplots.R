#Make Multi-axis boxplots for each region

library(ggplot2)
library(gridExtra)
library(ggplotify)
library(cowplot)
library(dplyr)

setwd("C:/Users/sarah/Documents/")

# load SISAL and ECHAM boxplot data
dat_all <- read.csv("PhD/analyses/data/3_degs_boxplot_dat.csv") # cols = site_id, name, lat, lon, region, time_slice, value, dat_type (spel, ppt, temp, echam d18O)

# order variables for plotting
dat_all$time_slice <- factor(dat_all$time_slice, c("MH","LGM","LIG"))
dat_all$dat_type <- factor(dat_all$dat_type, c("spel_d18O","Echam_d18O","Echam_pre","Echam_temp"))

# format axis labels
scaleFUN <- function(x) sprintf("%.1f", x)

# set colour palette
palette <- c("#F5793A",
             "#A95AA1",
             "#85C0F9",
             "#0F2080") 

# ISM plot
# transform Echam d18O to match sisal d18O
# transform Echam ppt to match sisal d18O
dat_all[which(dat_all$region == "ISM" & dat_all$dat_type=="Echam_pre"),]$val <- (dat_all[which(dat_all$region == "ISM" & dat_all$dat_type=="Echam_pre"),]$val)/-10
dat_all[which(dat_all$region == "ISM" & dat_all$dat_type=="Echam_temp"),]$val <- ((dat_all[which(dat_all$region == "ISM" & dat_all$dat_type=="Echam_temp"),]$val)*-1)-2

p1a <- ggplot(filter(dat_all, region == "ISM"), aes(x = time_slice, y = val, fill = dat_type)) +
  geom_boxplot(outlier.colour = "grey", lwd = 0.1, outlier.size = 0.5) + scale_fill_manual(values = palette) +
  scale_y_reverse(sec.axis = sec_axis(trans = ~.*-10, name = "precipitation"), labels = scaleFUN) + 
  theme_classic() +
  theme(legend.position = "none", axis.line.y.left = element_line(colour = palette[1], size = 0.5), axis.line.y.right = element_line(colour = palette[3], size = 0.5), axis.line.x.bottom = element_line(size = 0.2),
        axis.title = element_blank(), text = element_text(size = 8)) 

p1b <- ggplot(filter(dat_all, region == "ISM"), aes(x = time_slice, y= val)) +
  scale_y_reverse(sec.axis = sec_axis(trans = ~ (.+2)/(-1), name = "",labels = scaleFUN)) +
  labs(x = "", y = "") +
  theme(axis.text.y.right = element_text(color = "black", size = 6), 
        axis.ticks.y.right = element_line(color = "black"),
        axis.line.y.right = element_line(color = palette[4], size = 0.5),
        axis.title.y.right = element_text(color = "black") ,
        
        axis.text.y.left = element_blank(),
        axis.line.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        
        panel.background = element_blank(),
        plot.margin = margin(1, 0, 1, 1, "mm"))

p1c <- ggplot(filter(dat_all, region == "ISM"), aes(x = time_slice, y = val)) +
  scale_y_reverse(labels = scaleFUN) +
  labs(x = "", y = "") +
  theme(axis.line.y = element_line(color = palette[2], size = 0.5),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        
        panel.background = element_blank(),
        plot.margin = margin(1, 0, 1, 1, "mm"),
        text = element_text(size = 8))

p1 <- ggplotGrob(p1a)
q1 <- ggplotGrob(p1b)
u1 <- ggplotGrob(p1c)
q1$heights <- p1$heights
u1$heights <- p1$heights
plot_1 <- grid.arrange(u1,p1,q1,widths=c(1,9,1))
plot_1 <- as.ggplot(plot_1)


# EAM
# transform Echam d18O to match sisal d18O
dat_all[which(dat_all$region == "EAM" & dat_all$dat_type=="Echam_d18O"),]$val <- (dat_all[which(dat_all$region == "EAM" & dat_all$dat_type=="Echam_d18O"),]$val)-3
dat_all[which(dat_all$region == "EAM" & dat_all$dat_type=="Echam_pre"),]$val <- ((dat_all[which(dat_all$region == "EAM" & dat_all$dat_type=="Echam_pre"),]$val)/-7)-4
dat_all[which(dat_all$region == "EAM" & dat_all$dat_type=="Echam_temp"),]$val <- ((dat_all[which(dat_all$region == "EAM" & dat_all$dat_type=="Echam_temp"),]$val)*-1)-3

p2a <- ggplot(filter(dat_all, region == "EAM"), aes(x = time_slice, y = val, fill = dat_type)) +
  geom_boxplot(outlier.colour = "grey", lwd = 0.1, outlier.size = 0.5) + scale_fill_manual(values = palette) +
  scale_y_reverse(sec.axis = sec_axis(trans = ~(.+4)*(-7), name = "precipitation"), labels = scaleFUN) + 
  theme_classic() +
  theme(legend.position = "none", axis.line.y.left = element_line(colour = palette[1], size = 0.5), axis.line.y.right = element_line(colour = palette[3], size = 0.5), axis.line.x.bottom = element_line(size = 0.5),
        axis.title = element_blank(), text = element_text(size = 8)) 

p2b <- ggplot(filter(dat_all, region == "EAM"), aes(x = time_slice, y= val)) +
  scale_y_reverse(sec.axis = sec_axis(trans = ~ (.+3)/(-1), name = "", labels = scaleFUN)) +
  labs(x = "", y = "") +
  theme(axis.text.y.right = element_text(color = "black", size = 6), 
        axis.ticks.y.right = element_line(color = "black"),
        axis.line.y.right = element_line(color = palette[4], size = 0.5),
        axis.title.y.right = element_text(color = "black") ,
        
        axis.text.y.left = element_blank(),
        axis.line.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        
        panel.background = element_blank(),
        plot.margin = margin(1, 0, 1, 1, "mm"))

p2c <- ggplot(filter(dat_all, region == "EAM"), aes(x = time_slice, y = val+3)) +
  scale_y_reverse(labels = scaleFUN) +
  labs(x = "", y = "") +
  theme(axis.line.y = element_line(color = palette[2], size = 0.5),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        
        panel.background = element_blank(),
        plot.margin = margin(1, 0, 1, 1, "mm"),
        text = element_text(size = 8))


p2 <- ggplotGrob(p2a)
q2 <- ggplotGrob(p2b)
u2 <- ggplotGrob(p2c)
q2$heights <- p2$heights
u2$heights <- p2$heights
plot_2 <- grid.arrange(u2,p2,q2,widths=c(1,9,1))
plot_2 <- as.ggplot(plot_2)


# IAM
# transform Echam d18O to match sisal d18O
dat_all[which(dat_all$region == "IAM" & dat_all$dat_type=="Echam_d18O"),]$val <- (dat_all[which(dat_all$region == "IAM" & dat_all$dat_type=="Echam_d18O"),]$val)-2
dat_all[which(dat_all$region == "IAM" & dat_all$dat_type=="Echam_pre"),]$val <- ((dat_all[which(dat_all$region == "IAM" & dat_all$dat_type=="Echam_pre"),]$val)/-10)-2
dat_all[which(dat_all$region == "IAM" & dat_all$dat_type=="Echam_temp"),]$val <- ((dat_all[which(dat_all$region == "IAM" & dat_all$dat_type=="Echam_temp"),]$val)*-1)-1


p3a <- ggplot(filter(dat_all, region == "IAM"), aes(x = time_slice, y = val, fill = dat_type)) +
  geom_boxplot(outlier.colour = "grey", lwd = 0.1, outlier.size = 0.5) + scale_fill_manual(values = palette) +
  scale_y_reverse(sec.axis = sec_axis(trans = ~(.+2)*(-10), name = "precipitation"), labels = scaleFUN) + 
  theme_classic() +
  theme(legend.position = "none", axis.line.y.left = element_line(colour = palette[1], size = 0.5), axis.line.y.right = element_line(colour = palette[3], size = 0.5), axis.line.x.bottom = element_line(size = 0.5),
        axis.title = element_blank(), text = element_text(size = 8)) 

p3b <- ggplot(filter(dat_all, region == "IAM"), aes(x = time_slice, y= val)) +
  scale_y_reverse(sec.axis = sec_axis(trans = ~ (.+1)/(-1), name = "", labels = scaleFUN)) +
  labs(x = "", y = "") +
  theme(axis.text.y.right = element_text(color = "black", size = 6), 
        axis.ticks.y.right = element_line(color = "black"),
        axis.line.y.right = element_line(color = palette[4], size = 0.5),
        axis.title.y.right = element_text(color = "black") ,
        
        axis.text.y.left = element_blank(),
        axis.line.y.left = element_blank(),
        axis.ticks.y.left = element_blank(),
        
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        
        panel.background = element_blank(),
        plot.margin = margin(1, 0, 1, 1, "mm"))

p3c <- ggplot(filter(dat_all, region == "IAM"), aes(x = time_slice, y = val+2)) +
  scale_y_reverse(labels = scaleFUN) +
  labs(x = "", y = "") +
  theme(axis.line.y = element_line(color = palette[2], size = 0.5),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        
        panel.background = element_blank(),
        plot.margin = margin(1, 0, 1, 1, "mm"),
        text = element_text(size = 8))

p3 <- ggplotGrob(p3a)
q3 <- ggplotGrob(p3b)
u3 <- ggplotGrob(p3c)
q3$heights <- p3$heights
u3$heights <- p3$heights
plot_3 <- grid.arrange(u3,p3,q3,widths=c(1,9,1))
plot_3 <- as.ggplot(plot_3)


# export multiplot as pdf 
pdf("PhD/monsoon_paper/figs/boxplot_fig.pdf", width = 3.54, height = 5.91)
plot_grid(plot_2, plot_1, plot_3, ncol = 1, align = "v")
dev.off()


