# Example of non linear transform of MCMC samples
library(brms)
library(raster)
library(tmap)
library(tidyverse)
library(dplyr)
library(sf)
mytime <- format(Sys.time(), "%Y-%m-%d")

# load model 
load(file = "Data/All_model_combinations.rda")

#select best model
m1 <- dist_mod_100k_noSLR_DH

#Now extract all MCMC samples
samples <- posterior_samples(m1)
names(samples)
nmcmc <- nrow(samples)

# import prediction data set with all predictors resampled to 100km 
pa.raster <- raster('Data/Seagrass_PresenceAbsence_Raster.tif')
crs(pa.raster) <-  "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

#import prediction data set (pressure values aggregated at 100km) and transform to Mollewiede projection
prediction_data <- st_read("Data/2021-06-04_global_prediction_data.gpkg")

prediction_data <- st_transform(prediction_data, crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

names(prediction_data)[1:9] <- c("bioregion", "mean_Dd_100k", "mean_Np_100k",  "mean_Oa_100k",
                                 "mean_Ocp_100k", "mean_Ship_100k", "mean_Sst_100k", "Turbidity.Mean_100k",
                                 "Turbidity.Coefficient.of.Variance_100k") 
prediction_data$ID <- as.factor(row.names(prediction_data))
prediction_data$lh <- rep(as.factor("P")) # repeating a single life-history but unused as there is no interaction with life history and predictors in the final model  

# replace NA with 0 because NA means no pressure data extracted so no pressure in that area 
prediction_data <- prediction_data %>%
  replace(is.na(.), 0)
sg_coords <- sf::st_coordinates(prediction_data$geom)

#Predict to every location in 'r' and calculate the number of times that location is in the top 10% of samples
#predict for every posterior samples
preds <- posterior_linpred(m1, newdata = prediction_data)
dim(preds)

#Now find top q% for each row (each sample)
myfun <- function(x,q){
  qx <- quantile(x, probs = q)
  (x>qx) + 0
}

#apply function to every row to get top 10% of locations for each MCMC sample (90th percentile)
predquant <- t(apply(preds, 1, myfun, 0.1))
#Now sum columns and divide by nmcmc to get proportion of times that location is in top 10%
prop_top_10 <- colSums(predquant)/nmcmc

# do 1-prop top get lowest 10 % i.e. most rapidly declining because of ordinal regression linear predictor (i.e. most at risk)
prop_top <- 1-prop_top_10

#Add values back into our raster
rpred <- pa.raster
rpred[rpred == 1] <- prop_top
rpred[rpred<=0] <- NA

plot(rpred)
writeRaster(rpred, filename = paste0("Data/", mytime,"_risk_predictions.tif"))
