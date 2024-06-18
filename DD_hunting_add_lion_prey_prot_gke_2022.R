rm(list = ls())
setwd("C:/Users/g23b661/Desktop/DD_analysis/Kafue_DD_covs_2022")

library(raster)
library(rgdal)
library(sp)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(viridis)
library(broom)
library(cowplot)
library(sf)
library(move)

#########################################################
# PREY --------------------------------------------------
#########################################################

#covariates for HDS model as multi-band geotiff from GEE
gkecovs <- stack("gke_area_covariates.tif")
gkecovs[[7]] <- reclassify(gkecovs[[7]], cbind(NA, 50000))

#puku--------------------------------------------------------
# out = jags output from HDS model fitting, one species at a time
load("GKE_puku300_2024-04-20_1236_hds_out.RData")

# raster math, using output from jagsUI as coefficients in the linear model
# of herd density, then multiplying by group size + 1

#this format if jagsUI is not loaded
gke_puku_densities<-(exp(out1$mean$beta0 +
                        (log(gkecovs[[7]]/1000)*out1$mean$beta1) +
                        (gkecovs[[2]]/100 * out1$mean$beta2) +
                        (gkecovs[[1]] * out1$mean$beta3) +
                        (log(gkecovs[[8]]/1000)*out1$mean$beta4)) *
                        (out1$mean$lambda.group +1))
  
#gke_puku_densities[gke_puku_densities>100]<-100
gke_puku__df  <- as.data.frame(gke_puku_densities, xy = TRUE)

#impala----------------------------------------------

#be careful, because this will load into THE SAME OBJECT out that was 
#used previously with other species

load("GKE_impala300_2024-04-20_127_hds_out.RData")

gke_imp_densities<-(exp(out1$mean$beta0 +
                           (log(gkecovs[[7]]/1000)*out1$mean$beta1) +
                           (gkecovs[[2]]/100 * out1$mean$beta2) +
                           (gkecovs[[1]] * out1$mean$beta3) +
                           (log(gkecovs[[8]]/1000)*out1$mean$beta4)) *
                       (out1$mean$lambda.group +1))

#gke_imp_densities[gke_imp_densities>100]<-100
gke_imp__df  <- as.data.frame(gke_imp_densities, xy = TRUE)

#small antelope (duiker, grysbok, bushbuck------------------------------

load("GKE_smallant300_2024-04-20_947_hds_out.RData")

gke_smallant_densities<-(exp(out1$mean$beta0 +
                          (log(gkecovs[[7]]/1000)*out1$mean$beta1) +
                          (gkecovs[[2]]/100 * out1$mean$beta2) +
                          (gkecovs[[1]] * out1$mean$beta3) +
                          (log(gkecovs[[8]]/1000)*out1$mean$beta4)) *
                      (out1$mean$lambda.group +1))

#gke_smallant_densities[gke_smallant_densities>10]<-10
gke_smallant__df  <- as.data.frame(gke_smallant_densities, xy = TRUE)

# combined prey density and biomass---------------------------

gke_awd_prey_density <- gke_imp_densities + gke_puku_densities + gke_smallant_densities
gke_awd_prey_biomass <- (gke_imp_densities*31.9 + 
                        gke_puku_densities*37 + 
                        gke_smallant_densities * 14)

gke_awd_prey_density__df  <- as.data.frame(gke_awd_prey_density, xy = TRUE)
gke_awd_prey_biomass__df  <- as.data.frame(gke_awd_prey_biomass, xy = TRUE)
writeRaster(x=gke_awd_prey_biomass, filename = 'gke_prey_biomass.tif', format="GTiff", overwrite=TRUE)


#########################################################
#  DOGS -------------------------------------------------
#########################################################

DD_dog_locs <- read.csv("Minutely_WD_2022_with_HUNTS.csv")
summary(DD_dog_locs)
# filter out a small number w/o lat/long.  Only need to filter on one of these 
DD_dog_locs <- DD_dog_locs %>% filter(!is.na(DD_dog_locs$DR.Long))

DD_dog_spatial <- SpatialPointsDataFrame(coords=DD_dog_locs[,22:23], data = DD_dog_locs)
crs(DD_dog_spatial) <- "+proj=longlat +datum=WGS84 +no_defs"             

#########################################################
# PROTECTION --------------------------------------------
#########################################################


# park and gma are binary rasters, 1/0 for values inside and outside that protection level
park <- gkecovs[[9]]*2
gma <- gkecovs[[10]]

plot(park)
plot(gma)
protection <- park + gma
plot(protection)

#########################################################
# LIONS --------------------------------------------------
#########################################################

# load workspace Lion_UD_6Apr24 to avoid running dBBMM code
# from Goodheart et al.'s published analysis
# make naming format consistent
load("Lion_UD_6Apr24.RData")
LI2022 <- LI_2022_WM_Raster
LI2022_Dry <- LI_2022Dry_WM_Raster

# Rasters with dry season or full year lion UD for 2022 and 2023 lion data

# LI2022 
# LI2022_Dry will use these for DD tag analysis for best temporal match
# LI2023
# LI2023_Dry will use these for DD tag analysis for best temporal match

#remove extreme outliers on each end, less than 0.1% of the data
LI2022_Dry <- reclassify(LI2022_Dry, cbind(0, 1.000327e-18 , NA))
LI2022_Dry <- reclassify(LI2022_Dry, cbind(0.030021615,0.5, 0.030021615))
summary(LI2022_Dry)
lion_2022_dry_df<- as.data.frame(LI2022_Dry, xy = TRUE)

# Extract the new covariates at each dog location
# add to existing data 
# and write to csv

pointCoordinates <- DD_dog_locs[,23:24]
coordinates(pointCoordinates)= ~ DR.Long + DR.Lat
lion_extract <- raster::extract(LI2022_Dry, pointCoordinates)
prot_extract <- raster::extract(protection, pointCoordinates)
prey_biomass_extract <- raster::extract(gke_awd_prey_biomass, pointCoordinates)
prey_density_extract <- raster::extract(gke_awd_prey_density, pointCoordinates)
puku_extract <- raster::extract(gke_puku_densities, pointCoordinates)
impala_extract <- raster::extract(gke_imp_densities, pointCoordinates)
smallant_extract <- raster::extract(gke_smallant_densities, pointCoordinates)

DD_dog_locs <- cbind(DD_dog_locs, lion_extract, prot_extract, prey_biomass_extract, 
                     prey_density_extract, puku_extract, impala_extract, smallant_extract)


#########################################################
# PACK SIZE, COMPOSITION, DENNING------------------------
#########################################################

# create pack variable from filename variable
DD_dog_locs <- DD_dog_locs %>% mutate(Pack = case_when(str_detect(Filename, "Eden") ~ "Eden", 
                                      str_detect(Filename, "Lumbeya") ~ "Lumbeya",
                                      str_detect(Filename, "Mayukuyuku") ~ "Mayukuyuku",
                                      str_detect(Filename, "Musekwe") ~ "Musekwe",
                                      str_detect(Filename, "Shishamba") ~ "Shishamba",
                                      str_detect(Filename, "Kakumbi") ~ "Kakumbi",
                                      str_detect(Filename, "Manzi") ~ "Manzi",
                                      str_detect(Filename, "Liuwa") ~ "Liuwa"))

#check result
table(DD_dog_locs$Pack)

# read pack size, composition, denning data file (includes 2022 and 2023)
pack_props <- read.csv("pack_size_comp_den.csv")

DD_dog_locs <- left_join(DD_dog_locs, pack_props, by = 'Pack')

write.table(DD_dog_locs, file='DD_dogs_with_all_covs_gke_2022.csv', append=FALSE, sep= ',', row.names = FALSE, col.names=TRUE)

# examine autocorrelation structure
ar(DD_dog_locs$VeDBA.smoothed)

# examine raw relationship VeDBA v prey biomass
library(ggsci)
ggplot(data = DD_dog_locs, aes(x = prey_biomass_extract, y = VeDBA.smoothed, 
                               group = 1, colour = as.factor(Pack))) +
  geom_point(alpha = 0.2) +
  guides(colour=FALSE, alpha = FALSE) +
  theme_classic() +
  scale_colour_jco() +
  xlab(bquote('Prey Biomass ( '*'kg'~ km^-2*')')) +
  ylab("VeDBA")

ggsave("vedba_gke_2022.tif", device= 'tiff' )


