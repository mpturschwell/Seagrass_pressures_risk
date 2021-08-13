
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(brms)
library(patchwork)
library(ggplot2)
library(ggthemes)
mytime <- format(Sys.time(), "%Y-%m-%d")

# ----------------------------------------------------------------------------------------------------
# FIGURE 1
# ----------------------------------------------------------------------------------------------------
#read in spatial data
sg_dat <- st_read("Data/dat2000-2010_spatial.gpkg")
data("World")
World.mercator <- st_transform(World,"+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") 

# Cropping poles 
World_clip_bb <- st_bbox(c(xmin= -18601620, ymin = -6520048, xmax= 19601620, ymax = 8951339))

# make global plot of sites and trajectories 
world_clip <- tm_shape(World.mercator, bbox =World_clip_bb) + 
  tm_polygons() +
  tm_shape(sg_dat) + 
  tm_dots(
    col = "IRtrend_X50.",
    size = 0.3,
    shape = 1, border.lwd = 2.5,     # open circles 
    title = "Trend for 2000's",
    legend.show = TRUE,
    midpoint = NA,
    palette  = "RdYlBu", n = 5, contrast = c(0.2, 1),
    breaks=c(-Inf, -0.025, -0.005, 0.005, 0.025, Inf),
    labels = c("Rapidly declining", "Slowly declining", "Stable", "Slowly increasing", "Rapidly increasing"))+
  tm_credits("A", size = 1.5, position = c(0.04, 0.9))+
  tm_layout(legend.position = c(0.05, 0.12))+
  tm_scale_bar(breaks = c(0, 2500),
               position = c(0.03, 0.01),
               text.size = 1)  
world_clip
tmap_save(world_clip, filename = paste0("Figures/", mytime, "_Figure_1-A.png"), dpi = 400,
          width = 10, height = 8)

# Figure 1B - inset of North America
nth_am <- st_bbox(c(xmin= -8601620, ymin = 3000048, xmax= -5601620, ymax = 5551339))
Fig1B <- tm_shape(World, bbox = nth_am) + #tmaptools::bb(matrix(c(0,1000000,2000000,2000000),2,2))) +
  tm_polygons() +
  tm_shape(sg_dat) + 
  tm_dots(
    col = "IRtrend_X50.",
    size = 0.5,
    shape = 1, border.lwd = 2.5,
    title = "Trend for 2000's",
    legend.show = FALSE,
    midpoint = NA,
    palette  = "RdYlBu", n = 5, contrast = c(0.2, 1),
    breaks=c(-Inf, -0.025, -0.005, 0.005, 0.025, Inf),
    labels = c("Rapidly declining", "Slowly declining", "Stable", "Slowly increasing", "Rapidly increasing"))+
  tm_scale_bar(breaks = c(0, 1000),
               position = c(0.05, 0.01),
               text.size = 1)  +
  tm_credits("B", size = 1.5, position = c(0.04, 0.94))
Fig1B

# Figure 1C - inset of Europe and Med
med <- st_bbox(c(xmin= -1501620, ymin = 4000048, xmax= 2501620, ymax = 7551339))
Fig1C <- tm_shape(World.mercator, bbox =med) + 
  tm_polygons() +
  tm_shape(sg_dat) + 
  tm_dots(
    col = "IRtrend_X50.",
    size = 0.5,
    shape = 1, border.lwd = 2.5,
    title = "Trend for 2000's",
    legend.show = FALSE,
    midpoint = NA,
    palette  = "RdYlBu", n = 5, contrast = c(0.2, 1),
    breaks=c(-Inf, -0.025, -0.005, 0.005, 0.025, Inf),
    labels = c("Rapidly declining", "Slowly declining", "Stable", "Slowly increasing", "Rapidly increasing"))+
  tm_scale_bar(breaks = c(0, 1000),
               position = c(0.05, 0.01),
               text.size = 1)  +
  tm_credits("C", size = 1.5, position = c(0.04, 0.94))
Fig1C

Fig1BC <- tmap_arrange(Fig1B, Fig1C, nrow = 1)
Fig1BC
tmap_save(Fig1BC, filename = paste0("Figures/",mytime,"_Figure_1-BC.png"), width = 9, height = 6, dpi = 400)

# ----------------------------------------------------------------------------------------------------
# FIGURE 2 
# ----------------------------------------------------------------------------------------------------
dat2000 <- read_csv("Data/Data_seagrass_trends_pressures_2000-2010.csv", col_types = cols())
load(file = "Data/All_model_combinations.rda")

# Calculate one-sided probability of an effect
probs <- p_direction(
  dist_mod_100k_noSLR_DH,
  effects = c("all"),
  component = c("all"),
  parameters = NULL,
  method = "direct")

# Get probabilities for parameters to add to parameter plot 
addtoplot <- data.frame(probs$pd[5:15])

#plotting
par_plot <- mcmc_plot(dist_mod_100k_noSLR_DH, pars = c("_100k$", "_disc_"))+#, "^lh*")+# type = "areas")+
  theme_fivethirtyeight() +
  theme(axis.text.y = element_text(hjust = 0))+
  geom_vline(xintercept = 0, col = 'black', linetype = "dashed")+
  theme_bw()+
  xlim(-5,5)+
  xlab("Parameter estimate")
par_plot

plotme <- par_plot$data
plotme <- cbind(plotme, addtoplot)
colnames(plotme)[10] <- c("Probability")
plotme$Probability <- sprintf("%.2f", plotme$Probability)
levels(plotme$parameter) <- c("Destructive demersal fishing",
                              "Nutrient pollution",          
                              "Ocean acidification",         
                              "Organic chemical pollution",  
                              "Shipping",                    
                              "Extreme sea surface temperature events",
                              "Turbidity (mean)", 
                              "Turbidity (CV)",
                              "Opportunistic", 
                              "Mixed",
                              "Persistent")

fig2a <- plotme[1:8,] %>%
  arrange(desc(m)) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(parameter=factor(parameter, levels=parameter)) %>%   # This trick update the factor levels
  ggplot(aes(x=m, y=parameter)) +
  geom_pointrange(aes(xmin = ll, xmax = hh), color='grey27', shape=21, fatten = 2, size = 0.5) +
  geom_pointrange(aes(xmin = l, xmax = h), fill='dodgerblue1', color='grey52', shape=21, fatten = 2, size = 2) +
  geom_point(size=3, color="dodgerblue1") +
  theme_bw() +
  xlab("")+
  ylab("")+
  ggtitle("Pressures")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(-5,5)+
  geom_text(aes(x = 4, label = Probability))+
  geom_vline(aes(xintercept=0), linetype = 'dotted')
fig2a

fig2b <- plotme[9:11,] %>%
  arrange(desc(m)) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(parameter=factor(parameter, levels=parameter)) %>%   # This trick update the factor levels
  ggplot(aes(x=m, y=parameter)) +
  geom_pointrange(aes(xmin = ll, xmax = hh), color='grey27', shape=21, fatten = 2, size = 0.5) +
  geom_pointrange(aes(xmin = l, xmax = h), fill='blue', color='grey52', shape=21, fatten = 2, size = 2) +
  geom_point(size=3, color="dodgerblue1") +
  theme_bw() +
  xlab("Parameter estimate")+
  ylab("")+
  ggtitle("Life history")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(-5,5)+
  geom_text(aes(x = 4, label = Probability))+
  geom_vline(aes(xintercept=0), linetype = 'dotted')+
  coord_fixed(ratio = 1.3)
fig2b

fig2 <- fig2a/fig2b +  plot_annotation(tag_levels = 'A')
fig2
ggsave(path = "Figures/", filename = paste0(mytime,"_Figure_2.png"), 
       fig2, width = 8, height = 8, units = c("in"), dpi = 500)

# ----------------------------------------------------------------------------------------------------
# FIGURE 3
# ----------------------------------------------------------------------------------------------------
# Turbidity Variability panel 
fig3a_cv <- conditional_effects(dist_mod_100k_noSLR_DH, effects = "Turbidity.Coefficient.of.Variance_100k", categorical = TRUE, prob = c(0.75)) %>%  
  plot(theme = theme(panel.grid = element_blank()))

fig3a <- fig3a_cv$`Turbidity.Coefficient.of.Variance_100k:cats__`+ 
  scale_fill_brewer(palette = "RdYlBu")+ 
  scale_color_brewer(palette = "RdYlBu",  name = "Trajectory",labels = c("Rapidly declining", "Slowly declining", "Stable", "Slowly increasing", "Rapidly increasing"))+
  guides(fill = FALSE)+
  theme_clean()+
  theme(panel.grid = element_blank())+
  theme(legend.position = "none")+
  xlab("Turbidity variability")+
  ylab("Probability")+
  theme(axis.text.x = element_text(color = "grey20", size = 12),
        axis.text.y = element_text(color = "grey20", size = 12),
        axis.title.x = element_text(color = "grey20", size = 16),
        axis.title.y = element_text(color = "grey20", size = 16))
fig3a

# Demersal destructive fishing panel 
fig3b_dd <- conditional_effects(dist_mod_100k_noSLR_DH, effects = "mean_Dd_100k", categorical = TRUE, prob = c(0.75)) %>%  
  plot(theme = theme(panel.grid = element_blank()))
fig3b <- fig3b_dd$`mean_Dd_100k:cats__`+ 
  scale_fill_brewer(palette = "RdYlBu")+ 
  scale_color_brewer(palette = "RdYlBu",  name = "Trajectory",labels = c("Rapidly declining", "Slowly declining", "Stable", "Slowly increasing", "Rapidly increasing"))+
  guides(fill = FALSE)+
  theme_clean()+
  xlab("Destructive demersal fishing")+
  ylab("Probability")+
  theme(axis.text.x = element_text(color = "grey20", size = 12),
        axis.text.y = element_text(color = "grey20", size = 12),
        axis.title.x = element_text(color = "grey20", size = 16),
        axis.title.y = element_text(color = "grey20", size = 16))
fig3b

#Life history plot
fig3c_LH <- conditional_effects(dist_mod_100k_noSLR_DH, effects = "lh", categorical = TRUE, prob = c(0.75)) %>%  
  plot(theme = theme(panel.grid = element_blank()))

plot_dat <- fig3c_LH$`lh:cats__`$data

grey_cols <- brewer.pal(n = 6, "Greys")[3:6]
fig3c <- ggplot(plot_dat, aes(x = cats__, y = estimate__, colour = effect1__))+
  geom_point(size = 7, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin = lower__, ymax = upper__),
                width = 0.1,
                position=position_dodge(width=0.7))   +
  scale_colour_manual(values = grey_cols, name = "Life history",labels = c("Colonising", "Opportunistic", "Mixed", "Persistent"))+
  guides(fill = FALSE)+
  scale_x_discrete(labels = c("Rapidly declining", "Slowly declining", "Stable", "Slowly increasing", "Rapidly increasing"))+
  theme(legend.text=element_text(size=20))+
  theme_clean()+
  theme(panel.grid = element_blank())+
  xlab("Trajectory")+
  ylab("Probability")+
  theme(axis.text.x = element_text(color = "grey20", size = 12),
        axis.text.y = element_text(color = "grey20", size = 12),
        axis.title.x = element_text(color = "grey20", size = 16),
        axis.title.y = element_text(color = "grey20", size = 16))
fig3c

fig3 <- (fig3a | fig3b) / fig3c 
fig3 <- fig3 + plot_annotation(tag_levels = 'A')
fig3 

ggsave(path = "Figures/", filename = paste0(mytime, "_Figure_3-ABC.png"), 
       fig3, width = 10, height = 8, units = c("in"), dpi = 300)

# -----------------------------------------------------------------------------------
# FIGURE 4
# -----------------------------------------------------------------------------------
rprep <- raster('Data/2021-08-12_risk_predictions.tif')
data("World")
World <- st_transform(World, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") 
World_clip_bb <- st_bbox(c(xmin= -12601620, ymin = -6520048, xmax= 17601620, ymax = 8751339))
crs(rpred) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

tmap_mode("plot")
world_clip_hotspots <- tm_shape(World, bbox = World_clip_bb) + 
  tm_polygons() +
  tm_shape(rpred, raster.warp = FALSE) + 
  tm_raster(palette = "Reds", n = 4,#, contrast = c(0.4,1),
            breaks=c(0, 0.25, 0.50, 0.75, 1),
            legend.reverse = TRUE,
            title = "Risk of rapid decline")+
  tm_layout(scale = 1.1, legend.position = c(0.82, 0.75))#+

world_clip_hotspots
tmap_save(world_clip_hotspots, filename =  paste0("Figures/",mytime, "_Figure_4.png"), dpi = 400,
          width = 10, height = 10)

# -----------------------------------------------------------------------------------
#         FIGURE S1 - PROPORTIONAL BARCHART 
# -----------------------------------------------------------------------------------
dat2000 <- read_csv("Data/Data_seagrass_trends_pressures_2000-2010.csv", col_types = cols())

dat2000 <- dat2000 %>% 
  group_by(bioregion) %>% 
  mutate(bioregion_count = n()) %>% 
  ungroup() %>% 
  mutate(bioregion_updated = paste0(bioregion, ": n=", bioregion_count))

prop2000 <- ggplot(dat2000, aes(fill = category, x = decade))+
  geom_bar(position = position_fill(reverse = TRUE))+
  xlab("")+
  ylab("Proportion")+
  theme_clean()+
  scale_fill_brewer(name = "Trajectory", labels = c("Rapid decline", "Slow decline", "Stable", "Slow increase", "Rapid increase"), 
                    palette = "RdYlBu")+
  facet_wrap(~bioregion_updated)+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank())

prop2000
ggsave(filename = paste0("Figures/", mytime, "_Figure_S1_Bioregion_proportion_trajectory_barchart.png"), 
       prop2000, width = 10, height = 8, units = c("in"), dpi = 300)

# -----------------------------------------------------------------------------
# FIGURE S2 - insets of pressure maps
# -----------------------------------------------------------------------------
dd_rast_ag <- raster("Data/Dd_100km.tif")
dd_rast_ag[dd_rast_ag<=0] <- NA
plot(dd_rast_ag)
crs(dd_rast_ag) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

turb_cv_ag <- raster("Data/Turb_CV_100km.tif")
turb_cv_ag[turb_cv_ag<=0] <- NA 
plot(turb_cv_ag)
crs(turb_cv_ag) <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

#raster to match 
tmap_mode("plot") # 'view' is leaflet integrative viewer
data("World")
World <- st_transform(World, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs") 

# clip world to limit poles 
World_clip_bb <- st_bbox(c(xmin= -13801620, ymin = -6520048, xmax= 18601620, ymax = 8751339))

# RASTER MAPS 
europe <- st_bbox(c(xmin= -1501620, ymin = 3800048, xmax= 2501620, ymax = 8551339))
tm_shape(World, bbox =europe) +
  tm_polygons()

med_threat <-  tm_shape(dd_rast_ag,bbox =europe) + 
  tm_raster(palette = "YlOrRd", n = 4, contrast = c(0.4,1),
            #         breaks=c(0, 0.2, 0.6, Inf),
            labels = c("Low", "Medium", "High", "Very High"),
            legend.reverse = TRUE,
            title = "Destructive demersal\n fishing pressure")+
  tm_shape(World, bbox =europe) +
  tm_polygons()+
  tm_layout(legend.format = list(scientific = TRUE),
            legend.position = c("right","center"),
            legend.title.size = 1,
            legend.text.size = 0.8,
            legend.bg.color = "white",
            legend.bg.alpha = 1)+
  tm_scale_bar(breaks = c(0, 1000),
               position = c(0.3, 0.01),
               text.size = 1)+
  tm_credits("(A)", size = 2, position = c(0.05, 0.80))
med_threat


med_turb <-  tm_shape(turb_cv_ag,bbox =europe) + 
  tm_raster(palette = "YlOrRd", n = 4, contrast = c(0.4,1),
            #         breaks=c(0, 0.2, 0.6, Inf),
            labels = c("Low", "Medium", "High", "Very High", "Extreme"),
            legend.reverse = TRUE,
            title = "Turbidity variability")+
  tm_shape(World, bbox =europe) +
  tm_polygons()+
  tm_layout(legend.format = list(scientific = TRUE),
            legend.position = c("right","center"),
            legend.title.size = 1,
            legend.text.size = 0.8,
            legend.bg.color = "white",
            legend.bg.alpha = 1)+
  tm_scale_bar(breaks = c(0, 1000),
               position = c(0.3, 0.01),
               text.size = 1)+
  tm_credits("(B)", size = 2, position = c(0.05, 0.80))
med_turb

# USA RASTER MAPS 
usa <- st_bbox(c(xmin= -7601620, ymin = 4000048, xmax= -3601620, ymax = 7051339))
tm_shape(World, bbox =usa) +
  tm_polygons()

usa_threat <-  tm_shape(dd_rast_ag,bbox =usa) + 
  tm_raster(palette = "YlOrRd", n = 4, contrast = c(0.4,1),
            #         breaks=c(0, 0.2, 0.6, Inf),
            labels = c("Low", "Medium", "High", "Very High"),
            legend.reverse = TRUE,
            title = "Destructive demersal\n fishing pressure")+
  tm_shape(World, bbox =usa) +
  tm_polygons()+
  tm_layout(legend.format = list(scientific = TRUE),
            legend.position = c("left","center"),
            legend.title.size = 1,
            legend.text.size = 0.8,
            legend.bg.color = "white",
            legend.bg.alpha = 1)+
  tm_scale_bar(breaks = c(0, 1000),
               position = c(0.3, 0.01),
               text.size = 1)+
  tm_credits("(C)", size = 2, position = c(0.05, 0.80))
usa_threat

usa_turb <-  tm_shape(turb_cv_ag,bbox =usa) + 
  tm_raster(palette = "YlOrRd", n = 4, contrast = c(0.4,1),
            #         breaks=c(0, 0.2, 0.6, Inf),
            labels = c("Low", "Medium", "High", "Very High", "Extreme"),
            legend.reverse = TRUE,
            title = "Turbidity variability")+
  tm_shape(World, bbox =usa) +
  tm_polygons()+
  tm_layout(legend.format = list(scientific = TRUE),
            legend.position = c("left","center"),
            legend.title.size = 1,
            legend.text.size = 0.8,
            legend.bg.color = "white",
            legend.bg.alpha = 1)+
  tm_scale_bar(breaks = c(0, 1000),
               position = c(0.3, 0.01),
               text.size = 1)+
  tm_credits("(D)", size = 2, position = c(0.05, 0.80))
usa_turb

ALL <- tmap_arrange(med_threat, med_turb, usa_threat, usa_turb, nrow = 2, ncol = 2)
tmap_save(ALL, filename = paste0("Figures/", mytime, "_Figure_S2.png"), width = 8, height = 8, dpi = 400)

#--------------------------------------------------------#
# FIGURE S3 - parameters with only GAM data (>2 obs)
#--------------------------------------------------------#
rm(list=ls())
mytime <- format(Sys.time(), "%Y-%m-%d")
load(file = "Data/All_models-onlyGAMdat.rda")

probs <- p_direction(
  dist_mod_100k_noSLR_DH,
  effects = c("all"),
  component = c("all"),
  parameters = NULL,
  method = "direct")

addtoplot <- data.frame(probs$pd[5:15])

par_plot <- mcmc_plot(dist_mod_100k_noSLR_DH, pars = c("_100k$", "_disc_"))+
  theme_fivethirtyeight() +
  theme(axis.text.y = element_text(hjust = 0))+
  geom_vline(xintercept = 0, col = 'black', linetype = "dashed")+
  theme_bw()+
  xlim(-5,5)+
  xlab("Parameter estimate")
par_plot

plotme <- par_plot$data
plotme <- cbind(plotme, addtoplot)
colnames(plotme)[10] <- c("Probability")
plotme$Probability <- sprintf("%.2f", plotme$Probability)
levels(plotme$parameter) <- c("Destructive demersal fishing",
                              "Nutrient pollution",          
                              "Ocean acidification",         
                              "Organic chemical pollution",  
                              "Shipping",                    
                              "Extreme sea surface temperature events",
                              "Turbidity (mean)", 
                              "Turbidity (CV)",
                              "Opportunistic", 
                              "Mixed",
                              "Persistent")
fS3a <- plotme[1:8,] %>%
  arrange(desc(m)) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(parameter=factor(parameter, levels=parameter)) %>%   # This trick update the factor levels
  ggplot(aes(x=m, y=parameter)) +
  geom_pointrange(aes(xmin = ll, xmax = hh), color='grey27', shape=21, fatten = 2, size = 0.5) +
  geom_pointrange(aes(xmin = l, xmax = h), fill='dodgerblue1', color='grey52', shape=21, fatten = 2, size = 2) +
  geom_point(size=3, color="dodgerblue1") +
  theme_bw() +
  xlab("")+
  ylab("")+
  ggtitle("Pressures")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(-5,5)+
  geom_text(aes(x = 4, label = Probability))+
  geom_vline(aes(xintercept=0), linetype = 'dotted')
fS3a

fS3b <- plotme[9:11,] %>%
  arrange(desc(m)) %>%    # First sort by val. This sort the dataframe but NOT the factor levels
  mutate(parameter=factor(parameter, levels=parameter)) %>%   # This trick update the factor levels
  ggplot(aes(x=m, y=parameter)) +
  geom_pointrange(aes(xmin = ll, xmax = hh), color='grey27', shape=21, fatten = 2, size = 0.5) +
  geom_pointrange(aes(xmin = l, xmax = h), fill='blue', color='grey52', shape=21, fatten = 2, size = 2) +
  geom_point(size=3, color="dodgerblue1") +
  theme_bw() +
  xlab("Parameter estimate")+
  ylab("")+
  ggtitle("Life history")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlim(-5,5)+
  geom_text(aes(x = 4, label = Probability))+
  geom_vline(aes(xintercept=0), linetype = 'dotted')+
  coord_fixed(ratio = 1.3)
fS3b
fS3 <- fS3a/fS3b +  plot_annotation(tag_levels = 'A')
fS3
ggsave(path = "Figures/", filename = paste0(mytime,"_Figure_S3-GAM-only-parameter-plot.png"), 
       fS3, width = 8, height = 8, units = c("in"), dpi = 500)

#--------------------------------------------------------#
# FIGURE S7 - Frequency histograms of eight pressures 
#--------------------------------------------------------#
model_dat <- dist_mod_100k_noSLR_DH$data
covariates <- model_dat[2:9]
colnames(covariates) <- c("Demersal destructive fishing", "Nutrient pollution", "Ocean acidification",
                          "Organic chemical pollution", "Shipping", "Extreme sea surface temperature events",
                          "Turbidity (mean)", "Turbidity (Variability)")
longcovs <- pivot_longer(covariates, cols = everything())
longcovs$name <- as.factor(longcovs$name)

hist_pressures <-ggplot(data = longcovs, aes(x = value))+
  geom_histogram(bins = 50)+
  facet_wrap(~name)+
  xlab("Pressure value")+
  ylab("Count")+
  theme_clean()

ggsave(filename = paste0("Figures/", mytime, "_Figure_S7.png"), 
       hist_pressures, width = 10, height = 8, units = c("in"), dpi = 300)

#--------------------------------------------------------#
# FIGURE - S8 observed vs predicted risk map trajectories
#--------------------------------------------------------#
#import raster of risk ranking 
risk_prob <- rast("Data/2021-08-12_risk_predictions.tif")
plot(risk_prob)
# import seagrass sites
sg_dat <- st_read("Data/dat2000-2010_spatial.gpkg")

sg_proj <- st_transform(sg_dat, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
sg_proj_vec <- as(sg_proj, "Spatial")
v <- vect(sg_proj_vec)

out <- terra::extract(risk_prob, v, method = "bilinear")
hello <- cbind(sg_dat, out)

plot_dat <- hello %>% 
  as_tibble() %>% 
  dplyr::select('category', 'X2021.08.12_risk_predictions')

figS8 <-ggplot(plot_dat)+
  aes(x = factor(category, level = c('rapid_decline', 'slow_decline', 'stable',
                                     'slow_increase', 'rapid_increase')), y = X2021.08.12_risk_predictions)+
  geom_boxplot()+
  theme_clean()+
  ylab("Probability")+
  ylim(0,1)+
  xlab("Trajectory")+
  scale_x_discrete(labels = c('Rapid decline', 'Slow decline', 'Stable',
                              'Slow increase', 'Rapid increase'))
figS8
ggsave(filename = paste0("Figures/", mytime,"_Figure_S8.png"), 
       figS8, width = 8, height = 8, units = c("in"), dpi = 500)
