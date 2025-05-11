## MAIN ANALYSES AND PLOTS ##

## Author: Philippe Fernandez-Fournier, 2025

#setwd("D:/sfuvault/SFU/MooersLab/Ch3_WinnersLosers/Analysis/")
rm(list=ls())

# Add package versions
library(sf)
library(ape)
library(scales)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(colorRamps)
library(RColorBrewer)
library(ggridges)
library(lme4)
library(lmerTest)
library(brms)
library(bayesplot)
library(mgcv)
library(picante)
library(data.table)
library(phytools)
library(ggtree)
library(ggtreeExtra)
library(patchwork)
library(ggeffects)
library(marginaleffects)
library(splines)

## Load data 
assemb.data <- read.csv("./data/data_assemb_final.csv")
assemb.data$study <- as.factor(as.character(assemb.data$study))

comm.data <- read.csv("./data/data_pres_final.csv") 

comm.data <- comm.data %>% 
  dplyr::select(-taxon, -species_ass_int_all, -contrib_mpd)

write.csv(comm.data, "./data/data_pres_final.csv", row.names = FALSE)

species.data <- read.csv("./data/data_species_final.csv") 
species.data$study <- as.factor(species.data$study)

time_data <- read.csv("./data/data_time_final.csv")

phylo_fish  <- read.tree("./data/phylo_fish_PFF.tre")

## Load results
post_table <- read.csv("../Graphs/posteriors_table.csv")
gbif_obs <- read.csv("./results/gbif_occurrences_clean.csv")
slopes_summary_df_MPD <- read.csv("./results/slopes_summary_MPD.csv")
slopes_summary_df_CTI <- read.csv("./results/slopes_summary_CTI.csv")
species_sti <- read.csv("./results/species_STI_final.csv")

## Colors
col_fish <- "#32658EFF"
col_scale_sti <- colorRampPalette(c("royalblue3", "firebrick3"))
col_scale_mpd <- rev(brewer.pal(n = 11, name = "RdBu"))
col_scale_bpres <- c("#762A83", "#b179cb", "grey70", "#65c26d", "#1c7837")
col_biomes <- c("darkorange2", "darkorange2", "mediumseagreen", "mediumseagreen", "mediumseagreen", "mediumseagreen", "steelblue4")

# Variable colors
col_delta_mpd   <- col_scale_mpd[10]
col_delta_bpres <- col_scale_bpres[5]
col_delta_sti   <- col_scale_sti(11)[6] 

# Set visual limits for betas (max of 5/95% quantile)
hist(assemb.data$median_betaMPD, breaks = 100)
(mpd_lim <- round(max(quantile(assemb.data$median_betaMPD, probs = c(0.05, 0.95), na.rm = TRUE)), 3))

hist(assemb.data$betaSST_lm, breaks = 100)
(sst_lim <- round(max(quantile(assemb.data$betaSST_lm, probs = c(0.05, 0.95), na.rm = TRUE)), 1))

hist(assemb.data$median_betaCTI, breaks = 100)
(cti_lim <- round(max(quantile(assemb.data$median_betaCTI, probs = c(0.01, 0.99), na.rm = TRUE)), 2))

# Useful sites
# https://uclpg-msc-sgds.github.io/GEOG0125/bayesian-generalised-additive-models-gams.html
# https://r.qcbs.ca/workshop08/book-en/gam-with-interaction-terms.html


# NEW DATA ----------------------------------------------------------------

# Add slope posterior SDs (for weighting)
assemb.data$median_betaMPD_sd <- slopes_summary_df_MPD$sd[match(assemb.data$ass.name, slopes_summary_df_MPD$rarefyID)]
assemb.data$median_betaCTI_sd <- slopes_summary_df_CTI$sd[match(assemb.data$ass.name, slopes_summary_df_CTI$ass.name)]

# Biomes
# Split West and East Northern Atlantic
assemb.data$biome2 <- with(assemb.data, ifelse(biome == "Temperate_Northern_Atlantic" & longitude < -30, 
                                               "Temperate_Northwest_Atlantic",
                                               ifelse(biome == "Temperate_Northern_Atlantic" & longitude >= -30, 
                                                      "Temperate_Northeast_Atlantic", 
                                                      biome)))

# New dataset and combine Indo-Pacific and split Atlantic
assemb.data <- assemb.data %>% 
  mutate(biome_abbr = case_when(
    biome2 == "Eastern_Indo_Pacific" ~ "Indo-Pacific",
    biome2 == "Central_Indo_Pacific" ~ "Indo-Pacific",
    biome2 == "Western_Indo_Pacific" ~ "Indo-Pacific",
    biome2 == "Tropical_Atlantic" ~ "Tropical Atlantic",
    biome2 == "Temperate_Northwest_Atlantic" ~ "Temperate NW Atlantic",
    biome2 == "Temperate_Northeast_Atlantic" ~ "Temperate NE Atlantic",
    biome2 == "Temperate_Australasia" ~ "Temperate Australasia",
    biome2 == "Temperate_Northern_Pacific" ~ "Temperate N. Pacific",
    biome2 == "Arctic" ~ "Arctic",
    TRUE ~ as.character(biome2)  # In case any biome is missing
  ))

assemb.data$biome_abbr <- factor(assemb.data$biome_abbr, levels = c("Indo-Pacific", "Tropical Atlantic", 
                                                                    "Temperate N. Pacific", "Temperate Australasia", 
                                                                    "Temperate NW Atlantic", "Temperate NE Atlantic", 
                                                                    "Arctic"))
assemb.data <- assemb.data %>%
  mutate(biome_atl = case_when(biome_abbr %in% c("Temperate NW Atlantic", "Temperate NE Atlantic") ~ biome_abbr,
                               TRUE ~ "Other"))

assemb.data$biome_atl <- factor(assemb.data$biome_atl,
                                levels = c("Temperate NW Atlantic", "Temperate NE Atlantic", "Other"))



# Assemblage betaMPD weight relative to inverse of assemblage slope SD
assemb.data <- assemb.data %>% 
  mutate(weight_mpd = 1/median_betaMPD_sd)

assemb.data$w_mpd <- assemb.data$weight_mpd/mean(assemb.data$weight_mpd)
hist(assemb.data$w_mpd, breaks = 50)

# Assemblage betaCTI weight relative to inverse of assemblage slope SD
assemb.data <- assemb.data %>% 
  mutate(weight_cti = 1/median_betaCTI_sd)

assemb.data$w_cti <- assemb.data$weight_cti/mean(assemb.data$weight_cti)
hist(assemb.data$w_cti, breaks = 50)

# Remove betaSST outliers and save new df
hist(assemb.data$betaSST_lm)

bsst_q1 <- quantile(assemb.data$betaSST_lm, 0.25, na.rm = TRUE)
bsst_q3 <- quantile(assemb.data$betaSST_lm, 0.75, na.rm = TRUE)
bsst_iqr <- bsst_q3 - bsst_q1
bsst_low <- bsst_q1 - 3 * bsst_iqr
bsst_upp <- bsst_q3 + 3 * bsst_iqr

assemb_data_bsst <- subset(assemb.data, betaSST_lm >= bsst_low & betaSST_lm <= bsst_upp)
hist(assemb_data_bsst$betaSST_lm)

# Scale variables
assemb.data$MPDchange_scaled <- as.numeric(scale(assemb.data$median_betaMPD))
assemb.data$SRchange_scaled <- as.numeric(scale(assemb.data$median_betaSR))
assemb.data$betaJacc_scaled  <- as.numeric(scale(assemb.data$median_betaJacc))
assemb.data$betaPhyloSor_scaled  <- as.numeric(scale(assemb.data$betaPhyloSor_lm))
assemb.data$latitude_scaled <- as.numeric(scale(assemb.data$abs_latitude))
assemb.data$CTIchange_scaled <- as.numeric(scale(assemb.data$median_betaCTI))
assemb.data$SSTchange_scaled <- as.numeric(scale(assemb.data$betaSST_lm))
assemb.data$SSTmean_scaled <- as.numeric(scale(assemb.data$SST_mean))


# DATA STATS --------------------------------------------------------------

# Number of observations (rows with presence in comm.data)
length(which(comm.data$presence == 1))

# Range of years of observation
range(time_data$YEAR)

# Number of observed years
ny_dat <- comm.data %>% 
  select(ass.name, n.years) %>% 
  distinct()

median(ny_dat$n.years)
mean(ny_dat$n.years)
hist(ny_dat$n.years)

quantile(ny_dat$n.years, 0.99)

# Number of assemblage per biome (modified)
table(assemb.data$biome_atl)

# Number of studies
length(unique(assemb.data$study))

# Number of taxonomic species in fish species
length(unique(species.data$sp.name))
length(unique(comm.data$sp.name))

# Number of taxonomic orders in fish species
length(unique(species.data$order))

# Number of taxonomic families in fish species
length(unique(species.data$family))

# GBIF observations
colnames(gbif_obs)
nrow(gbif_obs)

mean(gbif_obs$year, na.rm = TRUE)
median(gbif_obs$year, na.rm = TRUE)
quantile(gbif_obs$year, probs = 0.05, na.rm = TRUE)

hist(gbif_obs$year) 

# Types of STI estimates
species_sti2 <- species_sti[species_sti$sp.name %in% comm.data$sp.name,]
table(species_sti2$source)

# View assemb_SR by temperature
data_viz <- comm.data %>%
  group_by(ass.name) %>%
  mutate(assemb_SR = n_distinct(species)) %>% 
  summarize(mean_ED = mean(ed.std),
            mean_SR = mean(assemb_SR)) %>% 
  ungroup()

ggplot(data_viz, aes(x = log(mean_SR), y = mean_ED)) +
  geom_point() +
  theme_bw()

lm(mean_ED ~ log(mean_SR), data = data_viz)

# Range of latitudes
range(assemb.data$abs_latitude)

# Histogram of mean temperatures
plot_temps <- ggplot(assemb.data, aes(x = SST_mean)) +
  geom_histogram(bins = 50, fill = col_fish) +
  theme_classic(); plot_temps


## Proportion of assemblages with a positive slope
# MPD ~ time
length(which(assemb.data$median_betaMPD > 0))
nrow(assemb.data)
length(which(assemb.data$median_betaMPD > 0)) / nrow(assemb.data)

# sesMPD ~ time
length(which(assemb.data$median_betaSESMPD > 0))
nrow(assemb.data)
length(which(assemb.data$median_betaSESMPD > 0)) / nrow(assemb.data)

# MNTD ~ time
length(which(assemb.data$median_betaMNTD > 0))
nrow(assemb.data)
length(which(assemb.data$median_betaMNTD > 0)) / nrow(assemb.data)

# CTI ~ time
length(which(assemb.data$median_betaCTI > 0))
nrow(assemb.data)
length(which(assemb.data$median_betaCTI > 0)) / nrow(assemb.data)

# count signs of betaCTI and betaSST combination
CTI_pos <- assemb.data$median_betaCTI > 0
SST_pos <- assemb.data$betaSST_lm > 0
table(CTI_positive = CTI_pos, SST_positive = SST_pos)

# count signs of betaCTI and betaSST combination
MPD_pos <- assemb.data$median_betaMPD > 0
SST_pos <- assemb.data$betaSST_lm > 0
table(MPD_positive = MPD_pos, SST_positive = SST_pos)


# Density plot of sampling year vs. GBIF observations
plot_years <- ggplot() +
  geom_density(data = gbif_obs,  aes(x = year), 
               color = "#F46D43", fill = "#FDAE61", linewidth = 1, alpha = 0.5) +
  geom_density(data = time_data, aes(x = YEAR), 
                            color = "#313695", fill = "#4575B4", linewidth = 1, alpha = 0.5) +
  theme_classic(); plot_years

ggsave(filename = '../Graphs/Supp_GBIF_OBS_year.jpeg', plot_years,
       dpi = 300, width = 5, height = 4, units = "in")



# MAPS -------------------------------------------------------------
#load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true"))
#save(NE_box, NE_countries, NE_graticules, lbl.X, lbl.Y, file = "./data/data_naturalearth.RData")
load("./data/data_naturalearth.RData")

# Read in MEOW regions from TNC source
meow <-read_sf("./MEOW/MEOW.shp") #downloaded from https://hub.arcgis.com/datasets/74b6ac5c8fc24dcb8abaad6428a5dfa4_0/

# Create spatial data (countries)
world <- st_as_sf(rnaturalearth::ne_countries(scale = "large", returnclass = "sf"))

# Convert assemblage-level points data
assemb_point <- st_as_sf(assemb.data, coords = c("longitude", "latitude"), crs = 4326) 
assemb_point$long <- st_coordinates(assemb_point)[, 1]
assemb_point$lat <- st_coordinates(assemb_point)[, 2]

## WORLD MAP
# Set up mapping
PROJ <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 
meow_rob <- st_transform(meow, crs = PROJ)

world_rob <- st_as_sf(NE_countries) %>% st_transform(crs = PROJ)
world_graticules_rob <- st_as_sf(NE_graticules) %>% st_transform(crs = PROJ)
world_box_rob   <- st_as_sf(NE_box) %>% st_transform(crs = PROJ)

# Transforming lbl.Y coordinates
lbl.Y$lon <- c(-160, -165, -170, -170, -170, -170, -170, -170,  # N 80-10
               -175,                                     # 0 
               -170, -170, -170, -170, -170, -170, -165, -160,  # S 10-80
               160, 165, 170, 170, 170, 170, 170, 170,   # N 80-10
               175,                                      # 0
               170, 170, 170, 170, 170, 170, 165, 160)  # S 10-80
lbl.Y_sf <- st_as_sf(lbl.Y, coords = c("lon", "lat"), crs = 4326)  
lbl.Y_prj <- st_transform(lbl.Y_sf, crs = PROJ)
lbl.Y.prj <- cbind(st_coordinates(lbl.Y_prj), lbl.Y)  
names(lbl.Y.prj)[1:2] <- c("X.prj", "Y.prj")

# Transforming lbl.X coordinates
lbl.X_sf <- st_as_sf(lbl.X, coords = c("lon", "lat"), crs = 4326)  
lbl.X_prj <- st_transform(lbl.X_sf, crs = PROJ)
lbl.X.prj <- cbind(st_coordinates(lbl.X_prj), lbl.X)  
names(lbl.X.prj)[1:2] <- c("X.prj", "Y.prj")

assemb_point_sf <- st_as_sf(assemb.data, coords = c("longitude", "latitude"), crs = 4326)  
assemb_point_rob <- st_transform(assemb_point_sf, crs = PROJ)
assemb_point_rob <- data.frame(long = st_coordinates(assemb_point_rob)[, 1], 
                               lat = st_coordinates(assemb_point_rob)[, 2],
                               ass.name = assemb_point_rob$ass.name)

assemb_point_rob <- left_join(assemb_point_rob, assemb.data, by = join_by(ass.name))


# World map
plot_assembs <- ggplot() +
  geom_sf(data=world_box_rob, color = "black", fill = "transparent", size = 0.25) +
  geom_sf(data=world_graticules_rob, linetype="dotted", color = "grey80", size = 0.1) +
  geom_point(data = assemb_point_rob, aes(x = long, y = lat), 
               size = 2.5, color = col_fish, alpha = 0.7) +
  geom_sf(data=world_rob, color = "black", fill = "grey80", size = 0.1) +
  geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color = "grey50", size = 3) +
  #geom_text(data = lbl.X.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size = 2) +
  theme_void() +
  theme(legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  guides(alpha = "none"); plot_assembs

ggsave("../Graphs/Fig01a_map_assemblages.jpeg", plot_assembs,
     dpi = 300, width = 10, height = 5, units = "in")


## ADD DENSITY PLOT OF ASSEMBLAGES ON RIGHT
plot_lats <- ggplot(assemb.data, aes(x = latitude)) +
  geom_density(fill = col_fish, alpha = 0.8, adjust = 2) +
  coord_flip() +
  scale_x_continuous(breaks = seq(-80, 80, by = 10), limits = c(-90, 90), expand = c(0, 0)) + 
  theme_test() +
  theme(axis.title.y = element_blank(),   
        axis.text.y = element_blank(),    
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(),    
        axis.ticks.x = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_line(color = NA),
        plot.margin = unit(c(1, 0, 1, 0), "in")); plot_lats

plot_f1 <- plot_assembs + plot_lats + plot_layout(widths = c(9, 1))

ggsave("../Graphs/Fig01_map_assemblages.jpeg", plot_f1,
       dpi = 300, width = 10, height = 7, units = "in")


# COLORED BY betaMPD
assemb_point_rob <- assemb_point_rob %>% 
  arrange(abs(median_betaMPD))# Arange by magnitude of betaMPD

map_mpd <- ggplot() +
  geom_sf(data=world_box_rob, color = "black", fill = "transparent", size = 0.25) +
  geom_sf(data=world_graticules_rob, linetype="dotted", color = "grey80", size = 0.1) +
  geom_point(data = assemb_point_rob, aes(x = long, y = lat, color = median_betaMPD), 
             size = 0.5, alpha = 0.5, shape = 19) +
  geom_sf(data=world_rob, color = "black", fill = "grey80", size = 0.1) +
  geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color = "grey50", size = 3) +
  scale_color_gradientn(colors = col_scale_mpd, name = expression(paste(beta["MPD"])),
                        limits = c(-mpd_lim, mpd_lim),
                        oob = scales::squish) + # squish any values outside the limits to the extreme colors
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  guides(alpha = "none"); map_mpd

ggsave("../Graphs/map_betaMPD.jpeg", map_mpd,
       dpi = 300, width = 10.5, height = 5, units = "in")

# COLORED BY betaCTI
assemb_point_rob <- assemb_point_rob %>% 
  arrange(abs(median_betaCTI))# Arange by magnitude of betaMPD

map_cti <- ggplot() +
  geom_sf(data=world_box_rob, color = "black", fill = "transparent", size = 0.25) +
  geom_sf(data=world_graticules_rob, linetype="dotted", color = "grey80", size = 0.1) +
  geom_point(data = assemb_point_rob, aes(x = long, y = lat, color = median_betaCTI), 
             size = 0.5, alpha = 0.7, shape = 19) +
  geom_sf(data=world_rob, color = "black", fill = "grey80", size = 0.1) +
  geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color = "grey50", size = 3) +
  scale_color_gradientn(colors = col_scale_mpd, name = expression(paste(beta["CTI"])),
                        limits = c(-cti_lim, cti_lim),
                        oob = scales::squish) + # squish any values outside the limits to the extreme colors
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  guides(alpha = "none"); map_cti

ggsave("../Graphs/map_betaCTI.jpeg", map_cti, dpi = 300, width = 10.5, height = 5, units = "in")

# COLORED BY betaSST
assemb_point_rob <- assemb_point_rob %>% 
  arrange(abs(betaSST_lm))# Arange by magnitude of betaMPD

map_sst <- ggplot() +
  geom_sf(data=world_box_rob, color = "black", fill = "transparent", size = 0.25) +
  geom_sf(data=world_graticules_rob, linetype="dotted", color = "grey80", size = 0.1) +
  geom_point(data = assemb_point_rob, aes(x = long, y = lat, color = betaSST_lm), 
             size = 0.5, alpha = 0.5, shape = 19) +
  geom_sf(data=world_rob, color = "black", fill = "grey80", size = 0.1) +
  geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color = "grey50", size = 2) +
  scale_color_gradientn(colors = col_scale_mpd,
                        name = "Change in SST\n (°C/decade)",
                        limits = c(-1.6, 1.6),
                        breaks = c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5),
                        labels = c("-1.5", "-1.0", "-0.5", "0", "+0.5", "+1.0", "+1.5"),
                        oob = scales::squish) +    
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 10),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  guides(alpha = "none"); map_sst

ggsave("../Graphs/map_betaSST_short.jpeg", map_sst, dpi = 300, width = 7.5, height = 4.2, units = "in")

# COLORED BY TEMPERATURE
assemb_point_rob <- assemb_point_rob %>% 
  arrange(abs(SST_mean))# Arange by magnitude of betaMPD

map_temp <- ggplot() +
  geom_sf(data=world_box_rob, color = "black", fill = "transparent", size = 0.25) +
  geom_sf(data=world_graticules_rob, linetype="dotted", color = "grey80", size = 0.1) +
  geom_point(data = assemb_point_rob, aes(x = long, y = lat, color = SST_mean), 
             size = 0.5, alpha = 0.8, shape = 19) +
  geom_sf(data=world_rob, color = "black", fill = "grey80", size = 0.1) +
  geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color = "grey50", size = 2) +
  scale_color_gradientn(colors = col_scale_sti(11), 
                        name = "Mean SST",
                        limits = c(-1, 31),
                        breaks = c(0, 5, 10, 15, 20, 25, 30),
                        labels = c("0°C", "5°C", "10°C", "15°C", "20°C", "25°C", "30°C"),
                        oob = scales::squish) + # squish any values outside the limits to the extreme colors
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 10),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  guides(alpha = "none"); map_temp

ggsave("../Graphs/map_temp.jpeg", map_temp, dpi = 300, width = 7.5, height = 4, units = "in")


# COLORED BY betaSR
assemb_point_rob <- assemb_point_rob %>% 
  arrange(abs(median_betaSR))# Arange by magnitude of betaMPD

map_temp <- ggplot() +
  geom_sf(data=world_box_rob, color = "black", fill = "transparent", size = 0.25) +
  geom_sf(data=world_graticules_rob, linetype="dotted", color = "grey80", size = 0.1) +
  geom_point(data = assemb_point_rob, aes(x = long, y = lat, color = median_betaSR), 
             size = 0.5, alpha = 0.7, shape = 19) +
  geom_sf(data=world_rob, color = "black", fill = "grey80", size = 0.1) +
  geom_text(data = lbl.Y.prj, aes(x = X.prj, y = Y.prj, label = lbl), color = "grey50", size = 3) +
  #geom_text(data = lbl.X.prj, aes(x = X.prj, y = Y.prj, label = lbl), color="grey50", size = 2) +
  scale_color_gradientn(colors = col_scale_mpd, name = expression(paste(beta["SR"])),
                        limits = c(-0.03, 0.03),
                        oob = scales::squish) + # squish any values outside the limits to the extreme colors
  theme_void() +
  theme(legend.position = "right",
        legend.text = element_text(size = 12),
        plot.margin = unit(c(0, 0.5, 0, 0), "cm")) +
  guides(alpha = "none"); map_temp

ggsave("../Graphs/map_betaSR.jpeg", map_temp, dpi = 300, width = 10.5, height = 5, units = "in")

# betaCTI ~ TEMPERATURE ---------------------------------------------------------------------
# Community Temperature Index

## Compare models with different k
models_gam_cti <- lapply(3:15, function(k) {
  gam(median_betaCTI ~ s(SSTmean_scaled, k = k), 
      weights = w_cti, # Weight relative to inverse of assemblage slope posteriors SD
      data = assemb.data, method = "REML")
})

# Name models
names(models_gam_cti) <- paste0("k_", 3:15)

# Compare k checks and AICs
lapply(models_gam_cti, k.check)
(cti_gam_AICs <- sapply(models_gam_cti, AIC))
plot(cti_gam_AICs, xaxt = "n", xlab = "Model", ylab = "AIC")
axis(1, at = 1:13, labels = names(cti_gam_AICs), las = 2)

## GAM (Generalized Additive Model)
gam_cti_temp <- gam(median_betaCTI ~ s(SSTmean_scaled, k = 7), weights = w_cti,
                   data = assemb.data, method = "REML")

summary(gam_cti_temp)

# Check diagnostics
par(mfrow = c(2,2))
gam.check(gam_cti_temp)
par(mfrow = c(1,1))

# RESULTS: convergence and stable fit.

# Store means and SDs for back-transformation
temp_mean <- mean(assemb.data$SST_mean, na.rm = TRUE)
temp_sd <- sd(assemb.data$SST_mean, na.rm = TRUE)

# Create new data for predictions
newdata_cti <- data.frame(SSTmean_scaled = seq((0 - temp_mean)/temp_sd, # starting at 10 degrees abs latitude
                                                max(assemb.data$SSTmean_scaled, na.rm = TRUE), length.out = 200))

# Predict smooth effect of Temperature, excluding the SR effect and random effect
pred_gam_cti <- predict.gam(gam_cti_temp, newdata = newdata_cti, se.fit = TRUE, type = "response")

# Store predictions
newdata_cti$fit <- pred_gam_cti$fit
newdata_cti$upper <- newdata_cti$fit + 1.96 * pred_gam_cti$se.fit
newdata_cti$lower <- newdata_cti$fit - 1.96 * pred_gam_cti$se.fit

# Back-transform Temperature for plotting
newdata_cti$SST_mean <- newdata_cti$SSTmean_scaled * temp_sd + temp_mean

# Plot grey
plot_ctitemp_gam2 <- ggplot(assemb.data, aes(x = SST_mean, y = median_betaCTI)) +
  geom_point(color = "grey60", size = assemb.data$w_cti, alpha = 0.2) +
  geom_hline(yintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = newdata_cti, aes(x = SST_mean, y = fit), color = "black", linewidth = 1) +
  geom_ribbon(data = newdata_cti, aes(x = SST_mean, y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  ylim(c(-cti_lim, cti_lim)) +
  scale_x_continuous(breaks = seq(0, 25, by = 5), trans = "reverse") + 
  ylab("Rate of change in CTI") +
  xlab("Assemblage temperature (°C)") +
  theme_bw() +
  coord_flip(); plot_ctitemp_gam2

ggsave(filename = '../Graphs/betaCTI_temp_gam.jpeg', plot_ctitemp_gam2,
       dpi = 300, width = 4, height = 7, units = "in")


# betaMPD ~ TEMPERATURE ---------------------------------------------

## Compare models with different k
models_gam_mpd <- lapply(3:15, function(k) {
  gam(median_betaMPD ~ s(SSTmean_scaled, k = k) + s(SRchange_scaled, k = 10), weights = w_mpd, 
      data = assemb.data, method = "REML")
})

# Name models
names(models_gam_mpd) <- paste0("k_", 3:15)

# Compare k checks and AICs
lapply(models_gam_mpd, k.check)
sapply(models_gam_mpd, AIC)
plot(sapply(models_gam_mpd, AIC), xaxt = "n", xlab = "Model", ylab = "AIC")
axis(1, at = 1:13, labels = names(sapply(models_gam_mpd, AIC)), las = 2)

## GAM (Generalized Additive Model)
gam_mpd_temp <- gam(median_betaMPD ~ s(SSTmean_scaled, k = 7) + s(SRchange_scaled, k = 10), weights = w_mpd, 
                   data = assemb.data, method = "REML")

summary(gam_mpd_temp)

# Check diagnostics
par(mfrow = (c(2,2)))
gam.check(gam_mpd_temp)
par(mfrow = (c(1,1)))

# Store means and SDs for back-transformation
temp_mean <- mean(assemb.data$SST_mean, na.rm = TRUE)
temp_sd <- sd(assemb.data$SST_mean, na.rm = TRUE)
mpd_mean <- mean(assemb.data$median_betaMPD, na.rm = TRUE)
mpd_sd <- sd(assemb.data$median_betaMPD, na.rm = TRUE)

# Create new data for predictions
newdata_mpd <- data.frame(SSTmean_scaled = seq((0 - temp_mean)/temp_sd, # starting at 0 degrees C
                                                max(assemb.data$SSTmean_scaled, na.rm = TRUE), length.out = 200),
                          SRchange_scaled = 0)  # Set SR change to zero to parse effect of temperature

# Predict smooth effect of temperature, excluding the SR effect and random effect
pred_gam_mpd <- predict.gam(gam_mpd_temp, newdata = newdata_mpd, se.fit = TRUE, type = "response", exclude = "s(SRchange_scaled)")

# Store predictions
newdata_mpd$fit <- pred_gam_mpd$fit
newdata_mpd$upper <- newdata_mpd$fit + 1.96 * pred_gam_mpd$se.fit
newdata_mpd$lower <- newdata_mpd$fit - 1.96 * pred_gam_mpd$se.fit

# Back-transform temperature and MPDchange for plotting
newdata_mpd$SST_mean <- newdata_mpd$SSTmean_scaled * temp_sd + temp_mean

# Plot grey
plot_mpdtemp_gam2 <- ggplot(assemb.data, aes(x = SST_mean, y = median_betaMPD)) +
  geom_point(color = "grey60", alpha = 0.2, size = assemb.data$w_mpd) +
  geom_hline(yintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = newdata_mpd, aes(x = SST_mean, y = fit), color = "black", linewidth = 1) +
  geom_ribbon(data = newdata_mpd, aes(x = SST_mean, y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  ylim(c(-0.003, 0.005)) +
  scale_x_continuous(breaks = seq(0, 25, by = 5), trans = "reverse") + 
  ylab("Rate of change in MPD") +
  xlab("Assemblage temperature (°C)") +
  theme_bw() +
  coord_flip(); plot_mpdtemp_gam2

ggsave(filename = '../Graphs/betaMPD_temp_gam.jpeg', plot_mpdtemp_gam2,
       dpi = 300, width = 4, height = 7, units = "in")


# betaCTI ~ betaSST ------------------------------------------------------

# Fit weighted linear model
lm_cti_sst <- lm(median_betaCTI ~ betaSST_lm, data = assemb_data_bsst, weights = w_cti)

# Create new data for predictions
newdata_cti_sst <- data.frame(betaSST_lm = seq(min(assemb_data_bsst$betaSST_lm, na.rm = TRUE),
                                               max(assemb_data_bsst$betaSST_lm, na.rm = TRUE),
                                               length.out = 200))

# Get predicted values and 95% confidence intervals
pred_cti_sst <- predict(lm_cti_sst, newdata_cti_sst, interval = "confidence")
newdata_cti_sst$fit <- pred_cti_sst[, "fit"]
newdata_cti_sst$lwr <- pred_cti_sst[, "lwr"]
newdata_cti_sst$upr <- pred_cti_sst[, "upr"]

# Plot: CTI vs SST trend with CI
plot_cti_vs_sst <- ggplot() +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "Positive/Positive"), 
            color = "transparent", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "Negative/Negative"), 
            color = "transparent", alpha = 0.2) +
  scale_fill_manual(values = c("Positive/Positive" = col_scale_mpd[10], 
                               "Negative/Negative" = col_scale_mpd[2])) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  geom_point(data = assemb_data_bsst, aes(x = betaSST_lm, y = median_betaCTI), 
             size = assemb_data_bsst$w_mpd, alpha = 0.2, color = "grey30") +
  #geom_ribbon(data = newdata_cti_sst, aes(x = betaSST_lm, ymin = lwr, ymax = upr), 
  #            inherit.aes = FALSE, alpha = 0.2) +
  #geom_line(data = newdata_cti_sst, aes(x = betaSST_lm, y = fit), 
  #          inherit.aes = FALSE, color = "black", linewidth = 1) +
  ylim(-0.07, 0.1) +
  xlab("Change in SST (°C/decade)") +
  ylab("Change in CTI") +
  theme_bw() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"));plot_cti_vs_sst


ggsave(filename = '../Graphs/betaCTI_betaSST.jpeg', plot_cti_vs_sst,
       dpi = 300, width = 5, height = 4.5, units = "in")



# betaMPD ~ betaSST ------------------------------------------------------

# Fit weighted linear model
lm_mpd_sst <- lm(median_betaMPD ~ betaSST_lm, data = assemb_data_bsst, weights = w_mpd)

summary(lm_mpd_sst)

# Create new data for predictions
newdata_mpd_sst <- data.frame(betaSST_lm = seq(min(assemb_data_bsst$betaSST_lm, na.rm = TRUE),
                                               max(assemb_data_bsst$betaSST_lm, na.rm = TRUE),
                                               length.out = 200))

# Get predicted values and 95% confidence intervals
pred_mpd_sst <- predict(lm_mpd_sst, newdata_mpd_sst, interval = "confidence")
newdata_mpd_sst$fit <- pred_mpd_sst[, "fit"]
newdata_mpd_sst$lwr <- pred_mpd_sst[, "lwr"]
newdata_mpd_sst$upr <- pred_mpd_sst[, "upr"]

# Create a new variable for quadrant coloring based on the signs of betaSST_lm and median_betaMPD
assemb_data_bsst$quadrant <- ifelse(assemb_data_bsst$betaSST_lm > 0 & assemb_data_bsst$median_betaMPD > 0, "Positive/Positive",
                                    ifelse(assemb_data_bsst$betaSST_lm < 0 & assemb_data_bsst$median_betaMPD < 0, "Negative/Negative", 
                                           NA))
assemb_quad <- subset(assemb_data_bsst, !is.na(quadrant))
rect_dat <- data.frame(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf)

# Plot: MPD vs SST trend with CI
plot_mpd_vs_sst <- ggplot() +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "Positive/Positive"), 
            color = "transparent", alpha = 0.2) +
  geom_rect(aes(xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = "Negative/Negative"), 
            color = "transparent", alpha = 0.2) +
  scale_fill_manual(values = c("Positive/Positive" = col_scale_mpd[10], 
                               "Negative/Negative" = col_scale_mpd[2])) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  geom_point(data = assemb_data_bsst, aes(x = betaSST_lm, y = median_betaMPD), 
             size = assemb_data_bsst$w_mpd, alpha = 0.2, color = "grey30") +
  #geom_ribbon(data = newdata_mpd_sst, aes(x = betaSST_lm, ymin = lwr, ymax = upr), 
  #            inherit.aes = FALSE, alpha = 0.2) +
  #geom_line(data = newdata_mpd_sst, aes(x = betaSST_lm, y = fit), 
  #          inherit.aes = FALSE, color = "black", linewidth = 1) +
  ylim(-0.005, 0.01) +
  xlab("Change in SST (°C/decade)") +
  ylab("Change in MPD") +
  theme_bw() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"));plot_mpd_vs_sst

ggsave(filename = '../Graphs/betaMPD_betaSST.jpeg', plot_mpd_vs_sst,
       dpi = 300, width = 5, height = 4.5, units = "in")



# SAVE COMBINED PLOT betaSST -------------------------------------------------------
plot_cti_vs_sst_nox <- plot_cti_vs_sst + xlab(NULL)

plot_betaSST_comb <- plot_cti_vs_sst_nox / plot_mpd_vs_sst +
  plot_annotation(tag_levels = 'a', theme = theme(plot.tag = element_text(face = "bold")))

ggsave(filename = '../Graphs/Fig05_betaCTI-MPD_betaSST.jpeg', plot_betaSST_comb,
       dpi = 300, width = 4.5, height = 7, units = "in")


# betaSST ~ TEMPERATURE ---------------------------------------------
## Compare models with different k
models_gam_sst <- lapply(3:15, function(k) {
  gam(betaSST_lm ~ s(SSTmean_scaled, k = k),
      data = assemb_data_bsst, method = "REML")
})

# Name models
names(models_gam_sst) <- paste0("k_", 3:15)

# Compare k checks and AICs
lapply(models_gam_sst, k.check)
sapply(models_gam_sst, AIC)
plot(sapply(models_gam_sst, AIC))

## GAM (Generalized Additive Model)
gam_sst_temp <- gam(betaSST_lm ~ s(SSTmean_scaled, k = 7), data = assemb_data_bsst, method = "REML")

summary(gam_sst_temp)

# Check diagnostics
par(mfrow = c(2,2))
gam.check(gam_sst_temp)
par(mfrow = c(1,1))

# Store means and SDs for back-transformation
temp_mean <- mean(assemb_data_bsst$SST_mean, na.rm = TRUE)
temp_sd <- sd(assemb_data_bsst$SST_mean, na.rm = TRUE)

# Create new data for predictions
newdata_sst <- data.frame(SSTmean_scaled = seq((0 - temp_mean)/temp_sd, 
                                                max(assemb_data_bsst$SSTmean_scaled, na.rm = TRUE), length.out = 200)) 

# Predict smooth effect of temperature, excluding the SR effect and random effect
pred_gam_sst <- predict.gam(gam_sst_temp, newdata = newdata_sst, se.fit = TRUE, type = "response")

# Store predictions
newdata_sst$fit <- pred_gam_sst$fit
newdata_sst$upper <- newdata_sst$fit + 1.96 * pred_gam_sst$se.fit
newdata_sst$lower <- newdata_sst$fit - 1.96 * pred_gam_sst$se.fit

# Back-transform temperature and SSTchange for plotting
newdata_sst$SST_mean <- newdata_sst$SSTmean_scaled * temp_sd + temp_mean

# Plot grey
plot_ssttemp_gam2 <- ggplot(assemb_data_bsst, aes(x = SST_mean, y = betaSST_lm)) +
  geom_point(color = "grey60", alpha = 0.2) +
  geom_hline(yintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = newdata_sst, aes(x = SST_mean, y = fit), color = "black", linewidth = 1) +
  geom_ribbon(data = newdata_sst, aes(x = SST_mean, y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  #geom_smooth(method = "lm", color = "black") +
  xlim(c(28.2, 0)) +
  ylim(c(-0.3, 2.5)) +
  ylab("Change in SST") +
  xlab("Assemblage temperature (degrees C)") +
  theme_bw() +
  coord_flip(); plot_ssttemp_gam2

ggsave(filename = '../Graphs/betaSST_temp_gam.jpeg', plot_ssttemp_gam2,
       dpi = 300, width = 4, height = 7, units = "in")

# Plot biomes
plot_ssttemp_gam3 <- ggplot(assemb_data_bsst, aes(x = SST_mean, y = betaSST_lm)) +
  geom_point(aes(color = biome_atl), alpha = 0.3) +
  geom_hline(yintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = newdata_sst, aes(x = SST_mean, y = fit), color = "black", linewidth = 1) +
  geom_ribbon(data = newdata_sst, aes(x = SST_mean, y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  scale_color_manual(name = "Biogeographic region:",
                     values = c("Temperate NW Atlantic" = "#26828EFF",
                                "Temperate NE Atlantic" = "#35B779FF",
                                "Other" = "grey60")) +
  xlim(c(28.2, 0)) +
  ylim(c(-0.3, 2.5)) +
  ylab("Change in SST") +
  xlab("Assemblage temperature (degrees C)") +
  theme_bw() +
  theme(legend.position = c(0.77, 0.09),
        legend.background = element_rect(fill = alpha('white', 0.6), color = "grey80")) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  coord_flip(); plot_ssttemp_gam3

ggsave(filename = '../Graphs/betaSST_temp_gam_biomes.jpeg', plot_ssttemp_gam3,
       dpi = 300, width = 4.5, height = 7, units = "in")






#### *** APPENDIX FIGURES *** ####

# STUDIES METADATA ------------------------------------------------------
meta <- read.csv("./data/bioTIMEmetadataSept18.csv") # Available at: https://biotime.st-andrews.ac.uk/download.php
assemb.data2 <- read.csv("./data/data_assemb_evol.csv")
assemb.data2$study <- as.factor(assemb.data2$study)

meta2 <- meta[meta$STUDY_ID %in% assemb.data2$study,
              c("STUDY_ID", "CLIMATE", "TITLE", "START_YEAR", "END_YEAR")]

write.csv(meta2, "./results/Supp_TableS1_studies.csv", row.names = FALSE)



# betaMPD ~ deltaSR -----------------------------------------------------

## LM
lm_mpd_sr <- lm(median_betaMPD ~ median_betaSR, data = assemb.data)

summary(lm_mpd_sr)

## Plot
plot_mpd_sr <- ggplot(data = assemb.data, aes(x = median_betaSR, y = median_betaMPD)) +
  geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey20", linetype = "dashed") +
  geom_point(alpha = 0.4, color = col_fish, size = 2) +
  geom_smooth(method = "lm", color = "black", linewidth = 1, se = FALSE) +
  xlab("Change in SR") +
  ylab("Change in MPD") +
  theme_bw(); plot_mpd_sr

ggsave(filename = '../Graphs/betaMPD_betaSR.jpeg', plot_mpd_sr,
       dpi = 300, width = 6, height = 5, units = "in")



# MPD ~ SR SIMULATIONS ----------------------------------------------------
spp_unique <- unique(species.data$sp.name)
all(spp_unique %in% phylo_fish[[1]]$tip.label) # check if they are all in phylo

sim.delta <- data.frame()

for(i in 1:10000){
  
  # tracker
  cat("\r", paste0(round(i/10000*100, 1), '%     '))
  
  
  # SR and species assemblage
  sr.i <- sample(5:100, 1) # sample between 5 and 100 species
  comm.i <- sample(spp_unique, sr.i, replace = FALSE)
  
  # phylo of community (sample tree from multiPhylo at random each time)
  phy_i <- keep.tip(phylo_fish[[sample(1:length(phylo_fish), 1)]], tip = comm.i)
  #plot(phy_i, cex = 0.7)
  
  # Calculate MPD
  coph.i <- cophenetic(phy_i)
  diag(coph.i) <- NA
  mpd.i <- mean(coph.i[lower.tri(coph.i)])
  
  data.i <- data.frame(SR = sr.i,
                       MPD = mpd.i)
  
  # Append to data
  sim.delta <- bind_rows(sim.delta, data.i)
  
}


jpeg(filename = paste0("../Graphs/sim_MPD_SR_correlation.jpeg"),
     res = 300, width = 6, height = 5, units = "in")

ggplot(data = sim.delta, aes(x = SR, y = MPD)) +
  geom_point(color = col_fish, alpha = 0.3) + 
  geom_hline(yintercept = mean(sim.delta$MPD), linetype = "dashed", linewidth = 1) +
  #geom_smooth(method = "lm") +
  theme_bw() +
  theme(legend.position = 'none', 
        axis.text = element_text(size = 10),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

dev.off()



# LINEAR betaMPD ~ TEMPERATURE (BRM) ---------------------------------------------
# study as covariate (no interaction)
formula_brm_mpd_lat <- as.formula('median_betaMPD ~ SSTmean_scaled + SRchange_scaled + (1 | study)')

## Bayesian model
brm_mpd_temp <- brm(formula_brm_mpd_lat, data = assemb.data,
                   family = gaussian(),                            # Likelihood (default for continuous outcomes)
                   prior = c(prior(normal(0, 1), class = "b"),     # Prior for fixed effects
                             prior(cauchy(0, 1), class = "sd")),   # Prior for random effects
                   control = list(adapt_delta = 0.95,
                                  max_treedepth = 15,
                                  stepsize = 0.01),
                   iter = 4000,                                    # Number of iterations (including warmup)
                   warmup = 1000,
                   chains = 4,
                   cores = 4,
                   refresh = 100)

# The estimate for SSTmean_scaled quantifies how much median_betaMPD is expected to change for a 1-SD increase 
# in SST, holding SRchange_scaled constant and after accounting for the random intercept variability across study.

save(brm_mpd_temp, file = "./results/brm_fit_MPD-temp.RData")
#load("./results/brm_fit_MPD-temp.RData")

(summ_brm_mpd_temp <- summary(brm_mpd_lat))


## DIAGNOSTICS

# Posterior Predictive Checks
pp_check(brm_mpd_temp)

# Residual Diagnostics
brm_mpd_temp_res <- residuals(brm_mpd_temp)

# Check the distribution of residuals
hist(brm_mpd_temp_res[, "Estimate"], breaks = 100) # Histogram of Residuals
qqnorm(brm_mpd_temp_res[, "Estimate"])                # QQ plot for normality
qqline(brm_mpd_temp_res[, "Estimate"], col = "red")   # A bit heavy-tailed...

# Compute Bayesian R-squared to summarize model fit
bayes_R2(brm_mpd_temp)

## PLOT
# Extract means of posteriors
alpha_mean <- summ_brm_mpd_temp$fixed[1, "Estimate"]      # Intercept (first row in 'fixed' part)
beta_mean <- summ_brm_mpd_temp$fixed[2, "Estimate"]       # Slope (second row in 'fixed' part)

## Plot the posterior distribution of the slope
post_brm_mpd_temp <- as_draws_df(brm_mpd_temp)
slope_df_mpd_temp <- data.frame(slope = post_brm_mpd_temp$b_SSTmean_scaled)
median(slope_df_mpd_temp$slope)
quantile(slope_df_mpd_temp$slope, probs = c(0.025, 0.975))

# deltaSR slope
median(post_brm_mpd_temp$b_SRchange_scaled)
quantile(post_brm_mpd_temp$b_SRchange_scaled, probs = c(0.025, 0.975))

## Plot posterior distribution of slope
color_scheme_set("darkgray")
plot_b_temp <- mcmc_areas(post_brm, pars = "b_SSTmean_scaled", prob = 0.95) +
  geom_vline(aes(xintercept = 0), color = "grey20", linetype = "dashed") +  
  xlim(c(-0.00001, 0.0001)) +
  labs(title = "Posterior distribution of the slope of betaMPD ~ temperature",
       x = "Slope estimate (shaded 95% CI)") +
  coord_cartesian(ylim = c(0.9, 1.00)) +  # Adjust the y-axis limits
  theme_classic(); plot_b_temp

ggsave(filename = '../Graphs/betaMPD_temp_beta.jpeg', plot_b_temp,
       dpi = 300, width = 6, height = 3.5, units = "in")


## Plot study random intercept estimates (ordered by posterior means)
post_brm_mpd_temp_r <- post_brm_mpd_temp[startsWith(colnames(post_brm_mpd_temp), "r_study")]

# Reorder columns in post_brm_mpd_temp_r
post_brm_mpd_temp_r <- post_brm_mpd_temp_r[, order(colMeans(post_brm_mpd_temp_r))]

study_num <- str_extract(colnames(post_brm_mpd_temp_r), "\\d+") # clean up the column names to have only the study number

color_scheme_set("darkgray")
plot_r_study <- mcmc_areas(post_brm_mpd_temp_r, prob = 0.8, prob_outer = 0.99) +
  geom_vline(aes(xintercept = 0), color = "grey20", linetype = "dashed") +  # Mean line
  labs(x = "Intercept estimate",
       y = "Study ID") +
  scale_y_discrete(labels = study_num) +  # Rename y-axis tick marks
  theme_bw(); plot_r_study

ggsave(filename = '../Graphs/betaMPD_temp_r_study.jpeg', plot_r_study,
       dpi = 300, width = 6, height = 10, units = "in")


# Extract the mean and standard deviation used for scaling
temp_mean <- attr(scale(assemb.data$SST_mean), "scaled:center")
temp_sd <- attr(scale(assemb.data$SST_mean), "scaled:scale")
temp_labels <- c(0, 5, 10, 15, 20 , 25)
temp_scaled_values <- (temp_labels - temp_mean) / temp_sd # Convert labels to scaled values

## Plot BRM fit
plot_mpdtemp_brm <- ggplot(assemb.data, aes(x = SSTmean_scaled, y = median_betaMPD)) +
  geom_point(color = "grey60", alpha = 0.2) +
  geom_hline(yintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_abline(intercept = alpha_mean, slope = beta_mean, 
              linewidth = 1, linetype = "solid", color = "black") +
  #ylim(c(-0.0075, 0.01)) +
  ylim(c(-0.005, 0.01)) +
  ylab("Rate of change in MPD") +
  xlab("Assemblage mean SST (°C)") +
  scale_x_continuous(breaks = temp_scaled_values, 
                     labels = temp_labels) +
  theme_bw() +
  theme(legend.position = "none"); plot_mpdtemp_brm

ggsave(filename = '../Graphs/betaMPD_temp_brm.jpeg', plot_mpdtemp_brm,
       dpi = 300, width = 6.5, height = 5, units = "in")



# GAM betaJacc ~ TEMPERATURE ------------------------------------------------------

## Compare models with different k
models_gam_betaJacc <- lapply(3:15, function(k) {
  gam(median_betaJacc ~ s(SSTmean_scaled, k = k), 
      data = assemb.data, method = "REML")
})

# Name models
names(models_gam_betaJacc) <- paste0("k_", 3:15)

# Compare k checks and AICs
lapply(models_gam_betaJacc, k.check)
(betaJacc_gam_AICs <- sapply(models_gam_betaJacc, AIC))
plot(betaJacc_gam_AICs, xaxt = "n", xlab = "Model", ylab = "AIC")
axis(1, at = 1:13, labels = names(betaJacc_gam_AICs), las = 2)

## GAM (Generalized Additive Model)
gam_betaJacc_temp <- gam(median_betaJacc ~ s(SSTmean_scaled, k = 4), data = assemb.data, method = "REML")

summary(gam_betaJacc_temp)

# Check diagnostics
par(mfrow = c(2,2))
gam.check(gam_betaJacc_temp)
par(mfrow = c(1,1))

# RESULTS: convergence and stable fit.

# Store means and SDs for back-transformation
temp_mean <- mean(assemb.data$SST_mean, na.rm = TRUE)
temp_sd <- sd(assemb.data$SST_mean, na.rm = TRUE)

# Create new data for predictions
newdata_betaJacc <- data.frame(SSTmean_scaled = seq((0 - temp_mean)/temp_sd, # starting at 10 degrees abs latitude
                                               max(assemb.data$SSTmean_scaled, na.rm = TRUE), length.out = 200))

# Predict smooth effect of Temperature, excluding the SR effect and random effect
pred_gam_betaJacc <- predict.gam(gam_betaJacc_temp, newdata = newdata_betaJacc, se.fit = TRUE, type = "response")

# Store predictions
newdata_betaJacc$fit <- pred_gam_betaJacc$fit
newdata_betaJacc$upper <- newdata_betaJacc$fit + 1.96 * pred_gam_betaJacc$se.fit
newdata_betaJacc$lower <- newdata_betaJacc$fit - 1.96 * pred_gam_betaJacc$se.fit

# Back-transform Temperature for plotting
newdata_betaJacc$SST_mean <- newdata_betaJacc$SSTmean_scaled * temp_sd + temp_mean

# plot with SEs
plot_jacctemp_gam <- ggplot(assemb.data, aes(x = SST_mean, y = median_betaJacc)) +
  geom_point(color = "grey60", alpha = 0.2) +
  geom_hline(yintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = newdata_betaJacc, aes(x = SST_mean, y = fit), color = "black", linewidth = 1) +
  geom_ribbon(data = newdata_betaJacc, aes(x = SST_mean, y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  ylim(c(-0.03, 0.06)) +
  scale_x_continuous(breaks = seq(0, 25, by = 5), trans = "reverse") + 
  xlab("Assemblage temperature (°C)") +
  ylab("Rate of change in Jaccard index") +
  theme_bw()  +
  coord_flip(); plot_jacctemp_gam

# Save the plot
ggsave(filename = '../Graphs/betaJacc_temp_gam.jpeg', plot_jacctemp_gam,
       dpi = 300, width = 4, height = 7, units = "in")


# GAM betaPhyloSor ~ TEMPERATURE ------------------------------------------------------

## Compare models with different k
models_gam_betaPhyloSor <- lapply(3:15, function(k) {
  gam(median_betaPhyloSor ~ s(SSTmean_scaled, k = k), 
      data = assemb.data, method = "REML")
})

# Name models
names(models_gam_betaPhyloSor) <- paste0("k_", 3:15)

# Compare k checks and AICs
lapply(models_gam_betaPhyloSor, k.check)
(betaPhyloSor_gam_AICs <- sapply(models_gam_betaPhyloSor, AIC))
plot(betaPhyloSor_gam_AICs, xaxt = "n", xlab = "Model", ylab = "AIC")
axis(1, at = 1:13, labels = names(betaPhyloSor_gam_AICs), las = 2)

## GAM (Generalized Additive Model)
gam_betaPhyloSor_temp <- gam(median_betaPhyloSor ~ s(SSTmean_scaled, k = 5), data = assemb.data, method = "REML")

summary(gam_betaPhyloSor_temp)

# Check diagnostics
par(mfrow = c(2,2))
gam.check(gam_betaPhyloSor_temp)
par(mfrow = c(1,1))

# RESULTS: convergence and stable fit.

# Store means and SDs for back-transformation
temp_mean <- mean(assemb.data$SST_mean, na.rm = TRUE)
temp_sd <- sd(assemb.data$SST_mean, na.rm = TRUE)

# Create new data for predictions
newdata_betaPhyloSor <- data.frame(SSTmean_scaled = seq((0 - temp_mean)/temp_sd, # starting at 10 degrees abs latitude
                                                    max(assemb.data$SSTmean_scaled, na.rm = TRUE), length.out = 200))

# Predict smooth effect of Temperature, excluding the SR effect and random effect
pred_gam_betaPhyloSor <- predict.gam(gam_betaPhyloSor_temp, newdata = newdata_betaPhyloSor, se.fit = TRUE, type = "response")

# Store predictions
newdata_betaPhyloSor$fit <- pred_gam_betaPhyloSor$fit
newdata_betaPhyloSor$upper <- newdata_betaPhyloSor$fit + 1.96 * pred_gam_betaPhyloSor$se.fit
newdata_betaPhyloSor$lower <- newdata_betaPhyloSor$fit - 1.96 * pred_gam_betaPhyloSor$se.fit

# Back-transform Temperature for plotting
newdata_betaPhyloSor$SST_mean <- newdata_betaPhyloSor$SSTmean_scaled * temp_sd + temp_mean

# plot with SEs
plot_phylosortemp_gam <- ggplot(assemb.data, aes(x = SST_mean, y = median_betaPhyloSor)) +
  geom_point(color = "grey60", alpha = 0.2) +
  geom_hline(yintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = newdata_betaPhyloSor, aes(x = SST_mean, y = fit), color = "black", linewidth = 1) +
  geom_ribbon(data = newdata_betaPhyloSor, aes(x = SST_mean, y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  ylim(c(-0.03, 0.06)) +
  scale_x_continuous(breaks = seq(0, 25, by = 5), trans = "reverse") + 
  #xlab("Assemblage temperature (°C)") +
  xlab(NULL) +
  ylab("Rate of change in PhyloSor") +
  theme_bw()  +
  coord_flip(); plot_phylosortemp_gam

# Save the plot
ggsave(filename = '../Graphs/betaPhyloSor_temp_gam.jpeg', plot_phylosortemp_gam,
       dpi = 300, width = 4, height = 7, units = "in")


## betaJacc and deltaPhyloSor plotted and saved together

Jacc_PS_plot <- (plot_jacctemp_gam * plot_phylosortemp_gam) +
  plot_annotation(tag_levels = 'a')  # Add (a) and (b) labels

# Save the combined plot
ggsave("../Graphs/Jacc_PhyloSor_temp_gam.jpeg", Jacc_PS_plot, width = 8, height = 7, dpi = 300)

# GAM betaSESMPD ~ TEMPERATURE ------------------------------------------------------

## Compare models with different k
models_gam_betaSESMPD <- lapply(3:15, function(k) {
  gam(median_betaSESMPD ~ s(SSTmean_scaled, k = k), 
      data = assemb.data, method = "REML")
})

# Name models
names(models_gam_betaSESMPD) <- paste0("k_", 3:15)

# Compare k checks and AICs
lapply(models_gam_betaSESMPD, k.check)
(betaSESMPD_gam_AICs <- sapply(models_gam_betaSESMPD, AIC))
plot(betaSESMPD_gam_AICs, xaxt = "n", xlab = "Model", ylab = "AIC")
axis(1, at = 1:13, labels = names(betaSESMPD_gam_AICs), las = 2)

## GAM (Generalized Additive Model)
gam_betaSESMPD_temp <- gam(median_betaSESMPD ~ s(SSTmean_scaled, k = 7), data = assemb.data, method = "REML")

summary(gam_betaSESMPD_temp)

# Check diagnostics
par(mfrow = c(2,2))
gam.check(gam_betaSESMPD_temp)
par(mfrow = c(1,1))

# RESULTS: convergence and stable fit.

# Store means and SDs for back-transformation
temp_mean <- mean(assemb.data$SST_mean, na.rm = TRUE)
temp_sd <- sd(assemb.data$SST_mean, na.rm = TRUE)

# Create new data for predictions
newdata_betaSESMPD <- data.frame(SSTmean_scaled = seq((0 - temp_mean)/temp_sd, # starting at 10 degrees abs latitude
                                                        max(assemb.data$SSTmean_scaled, na.rm = TRUE), length.out = 200))

# Predict smooth effect of Temperature, excluding the SR effect and random effect
pred_gam_betaSESMPD <- predict.gam(gam_betaSESMPD_temp, newdata = newdata_betaSESMPD, se.fit = TRUE, type = "response")

# Store predictions
newdata_betaSESMPD$fit <- pred_gam_betaSESMPD$fit
newdata_betaSESMPD$upper <- newdata_betaSESMPD$fit + 1.96 * pred_gam_betaSESMPD$se.fit
newdata_betaSESMPD$lower <- newdata_betaSESMPD$fit - 1.96 * pred_gam_betaSESMPD$se.fit

# Back-transform Temperature for plotting
newdata_betaSESMPD$SST_mean <- newdata_betaSESMPD$SSTmean_scaled * temp_sd + temp_mean

# plot with SEs
plot_SESMPD_temp_gam <- ggplot(assemb.data, aes(x = SST_mean, y = median_betaSESMPD)) +
  geom_point(color = "grey60", alpha = 0.2) +
  geom_hline(yintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = newdata_betaSESMPD, aes(x = SST_mean, y = fit), color = "black", linewidth = 1) +
  geom_ribbon(data = newdata_betaSESMPD, aes(x = SST_mean, y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  ylim(c(-0.01, 0.01)) +
  scale_x_continuous(breaks = seq(0, 25, by = 5), trans = "reverse") + 
  xlab("Assemblage temperature (°C)") +
  #xlab(NULL) +
  ylab("Rate of change in sesMPD") +
  theme_bw()  +
  coord_flip(); plot_SESMPD_temp_gam

# Save the plot
ggsave(filename = '../Graphs/betaSESMPD_temp_gam.jpeg', plot_SESMPD_temp_gam,
       dpi = 300, width = 4, height = 7, units = "in")


# GAM betaSR ~ TEMPERATURE ------------------------------------------------------

## Compare models with different k
models_gam_betaSR <- lapply(3:15, function(k) {
  gam(median_betaSR ~ s(SSTmean_scaled, k = k), 
      data = assemb.data, method = "REML")
})

# Name models
names(models_gam_betaSR) <- paste0("k_", 3:15)

# Compare k checks and AICs
lapply(models_gam_betaSR, k.check)
(betaSR_gam_AICs <- sapply(models_gam_betaSR, AIC))
plot(betaSR_gam_AICs, xaxt = "n", xlab = "Model", ylab = "AIC")
axis(1, at = 1:13, labels = names(betaSR_gam_AICs), las = 2)

## GAM (Generalized Additive Model)
gam_betaSR_temp <- gam(median_betaSR ~ s(SSTmean_scaled, k = 5), data = assemb.data, method = "REML")

summary(gam_betaSR_temp)

# Check diagnostics
par(mfrow = c(2,2))
gam.check(gam_betaSR_temp)
par(mfrow = c(1,1))

# RESULTS: convergence and stable fit.

# Store means and SDs for back-transformation
temp_mean <- mean(assemb.data$SST_mean, na.rm = TRUE)
temp_sd <- sd(assemb.data$SST_mean, na.rm = TRUE)

# Create new data for predictions
newdata_betaSR <- data.frame(SSTmean_scaled = seq((0 - temp_mean)/temp_sd, # starting at 10 degrees abs latitude
                                                        max(assemb.data$SSTmean_scaled, na.rm = TRUE), length.out = 200))

# Predict smooth effect of Temperature, excluding the SR effect and random effect
pred_gam_betaSR <- predict.gam(gam_betaSR_temp, newdata = newdata_betaSR, se.fit = TRUE, type = "response")

# Store predictions
newdata_betaSR$fit <- pred_gam_betaSR$fit
newdata_betaSR$upper <- newdata_betaSR$fit + 1.96 * pred_gam_betaSR$se.fit
newdata_betaSR$lower <- newdata_betaSR$fit - 1.96 * pred_gam_betaSR$se.fit

# Back-transform Temperature for plotting
newdata_betaSR$SST_mean <- newdata_betaSR$SSTmean_scaled * temp_sd + temp_mean

# plot with SEs
plot_SRtemp_gam <- ggplot(assemb.data, aes(x = SST_mean, y = median_betaSR)) +
  geom_point(color = "grey60", alpha = 0.2) +
  geom_hline(yintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = newdata_betaSR, aes(x = SST_mean, y = fit), color = "black", linewidth = 1) +
  geom_ribbon(data = newdata_betaSR, aes(x = SST_mean, y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  ylim(c(-0.02, 0.03)) +
  scale_x_continuous(breaks = seq(0, 25, by = 5), trans = "reverse") + 
  #xlab("Assemblage temperature (°C)") +
  xlab(NULL) +
  ylab("Rate of change in SR") +
  theme_bw()  +
  coord_flip(); plot_SRtemp_gam

# Save the plot
ggsave(filename = '../Graphs/betaSR_temp_gam.jpeg', plot_SRtemp_gam,
       dpi = 300, width = 4, height = 7, units = "in")


## betaJacc and deltaPhyloSor plotted and saved together

sesMPD_SR_plot <- (plot_SESMPD_temp_gam * plot_SRtemp_gam) +
  plot_annotation(tag_levels = 'a')  # Add (a) and (b) labels

# Save the combined plot
ggsave("../Graphs/Supp_sesMPD-SR_temp_gam.jpeg", sesMPD_SR_plot, width = 8, height = 7, dpi = 300)




# CTI ~ TIME --------------------------------------------------------------
## From LME
b_cti <- 0.003319641
int_cti <- 14.268313082

# Calculate assemblage means of STI and shifts through time
assemb.data <- comm.data %>%
  filter(!is.na(sti)) %>%
  filter(presence == 1) %>% # only present species
  group_by(ass.name, study, year) %>%
  summarise(cti = mean(sti, na.rm = TRUE)) %>% # Average STI of a given assemblage = CTI
  ungroup() %>%
  group_by(ass.name) %>%
  mutate(cYEAR = year - mean(year)) %>% # re-center year by assemblage
  ungroup()

colnames(assemb.data)

# Plot
plot_cti_time <- ggplot(assemb.data, aes(x = cYEAR, y = cti)) +
  geom_line(stat="smooth", method = "lm", aes(group = ass.name), se = FALSE, fullrange = FALSE,
            color = "black", alpha = 0.05, linewidth = 0.6) +
  #geom_smooth(method = "lm", color = col_scale_sti(11)[10], linewidth = 2) +
  geom_abline(slope = b_cti, intercept = int_cti, 
              color = col_scale_sti(11)[10], linewidth = 2) +
  labs(x = "Year of observation (centered)",
       y = "Assemblage CTI") +  
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 0.5, 0.5), "cm")); plot_cti_time # Adjust margins (top, right, bottom, left)

ggsave(filename = '../Graphs/CTI_time.jpeg', plot_cti_time,
       dpi = 300, width = 5, height = 5.5, units = "in")




# betaCTI ~ betaSR ------------------------------------------------------
## LM
lm_cti_sr <- lm(median_betaCTI ~ median_betaSR, data = assemb.data)

summary(lm_cti_sr)

# Plot
plot_cti_sr <- ggplot(data = assemb.data, aes(x = median_betaSR, y = median_betaCTI)) +
  geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey20", linetype = "dashed") +
  geom_point(alpha = 0.2, color = col_fish, size = 2) +
  geom_smooth(method = "lm", color = "black", linewidth = 1, se = FALSE) +
  xlab("Change in SR") +
  ylab("Change in CTI") +
  theme_bw(); plot_cti_sr

ggsave(filename = '../Graphs/betaCTI_betaSR.jpeg', plot_cti_sr,
       dpi = 300, width = 6, height = 5, units = "in")



# MPD ~ TIME --------------------------------------------------------------
# BRM slope and intercept (median of post dist)
b_mpd   <- post_table$b_median[post_table$metric == "MPD"]
int_mpd <- post_table$Int_median[post_table$metric == "MPD"]

# Plot
(pretty_mpd       <- pretty(time_data$MPD))
pretty_mpd_transf <- log(pretty_mpd)

plot_mpd <- ggplot(time_data, aes(x = cYEAR, y = log(MPD))) +
  geom_line(stat="smooth", method = "lm", aes(group = rarefyID), se = FALSE, fullrange = FALSE,
            color = "black", alpha = 0.05, linewidth = 0.6) +
  #geom_smooth(method = "lm", color = col_delta_mpd, linewidth = 2) +
  geom_abline(slope = b_mpd, intercept = int_mpd, 
              color = col_delta_mpd, linewidth = 2) +
  scale_y_continuous(breaks = pretty_mpd_transf, labels = pretty_mpd) +
  labs(x = "Year of observation (centered)",
       y = "Assemblage MPD (log scale)") +  
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 0.5, 0.5), "cm")); plot_mpd # Adjust margins (top, right, bottom, left)

ggsave(filename = '../Graphs/MPD_time.jpeg', plot_mpd,
       dpi = 300, width = 6, height = 5, units = "in")



# betaMPD ~ LAT colored by betaSST --------------------------------------------------

assemb_sst <- assemb.data %>% 
  arrange(betaSST_lm)

# Plot colored by deltaSST
plot_mpdlat_gam_sst <- ggplot(assemb_sst, aes(x = abs_latitude, y = median_betaMPD)) +
  geom_point(aes(color = betaSST_lm), alpha = 0.2) +
  geom_hline(yintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = newdata_mpd, aes(x = latitude, y = fit), color = "black", linewidth = 1) +
  geom_ribbon(data = newdata_mpd, aes(x = latitude, y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  xlim(c(10,62)) +
  ylab("Rate of change in MPD (per year)") +
  xlab("Absolute latitude") +
  scale_color_gradientn(name = expression(paste(beta["SST"])),
                        colors = col_scale_sst, 
                        limits = c(-0.6, 0.6),
                        oob = scales::squish) + # squish any values outside the limits to the extreme colors
  theme_bw(); plot_mpdlat_gam_sst

ggsave(filename = '../Graphs/betaMPD_lat_gam_betaSST.jpeg', plot_mpdlat_gam_sst,
       dpi = 300, width = 6.5, height = 5, units = "in")

# Plot colored by deltaSST FLIPPED
plot_mpdlat_gam_sst2 <- ggplot(assemb_sst, aes(x = abs_latitude, y = median_betaMPD)) +
  geom_point(aes(color = betaSST_lm), alpha = 0.2) +
  geom_hline(yintercept = 0, color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_line(data = newdata_mpd, aes(x = latitude, y = fit), color = "black", linewidth = 1) +
  geom_ribbon(data = newdata_mpd, aes(x = latitude, y = fit, ymin = lower, ymax = upper), alpha = 0.2) +
  ylim(c(-0.003, 0.005)) +
  xlim(c(10,62)) +
  ylab("Rate of change in MPD") +
  xlab("Absolute latitude") +
  scale_color_gradientn(name = expression(paste(beta["SST"])),
                        colors = col_scale_sst, 
                        limits = c(-0.6, 0.6),
                        oob = scales::squish) + # squish any values outside the limits to the extreme colors
  theme_bw()  +
  coord_flip(); plot_mpdlat_gam_sst2

ggsave(filename = '../Graphs/betaMPD_lat_gam_betaSST_FLIPPED.jpeg', plot_mpdlat_gam_sst2,
       dpi = 300, width = 4.5, height = 7, units = "in")


# betaCTI ~ betaMPD -----------------------------------------------------

## Plot
plot_mpd_cti <- ggplot(data = assemb.data, aes(x = median_betaMPD, y = median_betaCTI)) +
  geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey20", linetype = "dashed") +
  geom_point(alpha = 0.4, color = col_fish, size = 2) +
  geom_smooth(method = "lm", color = "black", linewidth = 1, se = FALSE) +
  #ylim(c(-0.06, 0.06)) +
  xlab("Change in MPD") +
  ylab("Change in CTI") +
  theme_bw(); plot_mpd_cti

ggsave(filename = '../Graphs/betaCTI_betaMPD.jpeg', plot_mpd_cti,
       dpi = 300, width = 6, height = 5, units = "in")

## Plot biomes
plot_mpd_cti2 <- ggplot(data = assemb.data, aes(x = median_betaMPD, y = median_betaCTI, group = biome_atl, color = biome_atl)) +
  geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "grey20", linetype = "dashed") +
  geom_point(alpha = 0.2, size = 1) +
  scale_color_manual(name = "",
                     values = c("Temperate NW Atlantic" = "#26828EFF",
                                "Temperate NE Atlantic" = "#35B779FF",
                                "Other" = "grey60")) +
  geom_smooth(method = "lm", linewidth = 1, se = FALSE, fullrange = FALSE) +
  #ylim(c(-0.06, 0.06)) +
  xlab("Change in MPD") +
  ylab("Change in CTI") +
  theme_bw() +
  theme(legend.position = "top"); plot_mpd_cti2

ggsave(filename = '../Graphs/betaCTI_betaMPD_biome.jpeg', plot_mpd_cti2,
       dpi = 300, width = 5, height = 4.5, units = "in")

# ED ~ STI --------------------------------------------------------
colnames(species.data)

# Create data, center STI
species_sti <- species.data %>%
  mutate(ed_trans = sqrt(ed.std)) %>% # transformed ed
  group_by(ass.name) %>%
  mutate(ed_trans_cent = (ed_trans - mean(ed_trans))) %>% # centered ed
  mutate(sti_cent = (sti - mean(sti))) %>% # centered STI
  ungroup() %>% 
  mutate(sti_cent_scaled = as.numeric(scale(sti_cent))) %>%  # scaled STI
  mutate(ed_cent_scaled = as.numeric(scale(ed_trans_cent))) # scaled STI

## LME
lme_ed_sti <- lmer(ed_trans_cent ~ sti_cent_scaled + (1 + sti_cent_scaled | ass.name), data = species_sti)

summary(lme_ed_sti)


# Plot Assemblage trends (BIOMES!)
plot_ed_sti2 <- ggplot(species_sti, aes(x = sti_cent, y = ed.std)) +
  geom_line(stat="smooth", method = "lm", aes(group = ass.name), se = FALSE, fullrange = FALSE,
            color = "black", alpha = 0.05, linewidth = 0.6) +
  geom_smooth(method = "lm", color = "royalblue3", linewidth = 2) +
  scale_y_continuous(breaks = c(0.01, 0.05, 0.10, 0.20), trans = scales::sqrt_trans()) +
  labs(y = "Local ED (square root scale)",
       x = "STI (centered)") +  
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 0.5, 0.5), "cm")); plot_ed_sti2 # Adjust margins (top, right, bottom, left)

ggsave(filename = '../Graphs/ed_sti_assemb.jpeg', plot_ed_sti2,
       dpi = 300, width = 6, height = 5, units = "in")



# New data 
species_tax_sti <- species.data %>%
  group_by(sp.name, sti) %>%
  summarize(mean_ed = mean(ed.std)) %>%
  ungroup

## Plot taxonomic species averages
plot_ed_sti <- ggplot(species_tax_sti, aes(x = sti, y = mean_ed)) +
  geom_point(alpha = 0.4, color = col_fish) + 
  geom_line(stat="smooth", method = "lm", se = FALSE, 
            color = "black", alpha = 1, linewidth = 1) +
  scale_y_continuous(breaks = c(0.01, 0.05, 0.10, 0.20), trans = scales::sqrt_trans()) +
  labs(x = "Species thermal affinity (STI)",
       y = "Average local ED (square root scale)") +
  theme_classic(); plot_ed_sti

ggsave(filename = '../Graphs/ed_sti_species.jpeg', plot_ed_sti,
       dpi = 300, width = 6, height = 5, units = "in")



# CTI ~ TIME * TEMP ---------------------------------------------------------------
colnames(time_data)

# Create data, center SST
time_cti <- time_data %>%
  filter(rarefyID %in% assemb.data$ass.name) %>% 
  group_by(rarefyID) %>% 
  mutate(sst_cent = sst - mean(sst)) %>% 
  mutate(yr_cent  = YEAR - mean(YEAR)) %>% 
  ungroup() %>% 
  mutate(sst_scaled = (sst - mean(sst_cent)) / sd(sst_cent)) %>% 
  mutate(yr_scaled  = (yr_cent - mean(yr_cent)) / sd(yr_cent))

# LMER with nonlinear effect of SST
lme_cti_sst3 <- lmer(cti ~ yr_scaled + sst_scaled + yr_scaled:ns(sst_scaled, df = 3) 
                     + (1 + yr_scaled | STUDY_ID) + (1 + yr_scaled | rarefyID), 
                     data = time_cti, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(lme_cti_sst3)

# for back-transformation
sst_mean <- mean(time_cti$sst_cent)
sst_sd <- sd(time_cti$sst_cent)
yr_mean <- mean(time_cti$yr_cent)
yr_sd <- sd(time_cti$yr_cent)

# SST values to plot
#sst_quants <- quantile(time_cti$sst_scaled, probs = c(0.01, 0.25, 0.50, 0.75, 0.99), na.rm = TRUE)
#(sst_values <- (sst_quants * sst_sd) + sst_mean)
sst_values <- c(0, 5, 10, 15, 20, 25, 30)
sst_groups <- (sst_values - sst_mean) / sst_sd

sst_labels <- setNames(paste0(round(sst_values, 0) , "°C"), sst_groups)

# year to plot 
range(time_cti$yr_cent)
yr_values <- c(-30, -20, -10, 0, 10, 20, 30)
yr_breaks <- (yr_values - yr_mean) / yr_sd

# Plot predicted values
cti_sst_pred <- ggpredict(lme_cti_sst3, terms = c("yr_scaled [all]", "sst_scaled [sst_groups]"))

plot_cti_sst <- ggplot(data = time_cti, aes(x = yr_scaled, y = cti)) +
  geom_point(color = "grey50", alpha = 0.1) +
  geom_line(data = cti_sst_pred,inherit.aes = FALSE, aes(x = x, y = predicted, color = group), 
            linewidth = 1) +  
  #geom_ribbon(data = cti_sst_pred, inherit.aes = FALSE, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), 
  #            alpha = 0.3, color = NA) +
  scale_x_continuous(breaks = yr_breaks, labels = yr_values) +
  scale_y_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30), labels = c("0°C", "5°C", "10°C", "15°C", "20°C", "25°C", "30°C")) +
  scale_color_manual(values = col_scale_sti(length(sst_labels)), labels = sst_labels) +
  scale_fill_manual(values = col_scale_sti(length(sst_labels)), labels = sst_labels) +
  guides(fill = "none") +
  labs(x = "Observation year (centered)",
       y = "Assemblage CTI",
       color = "Assemblage\nmean SST") + 
  theme_bw(); plot_cti_sst

ggsave(filename = '../Graphs/CTI_time_SST.jpeg', plot_cti_sst,
       dpi = 300, width = 5.5, height = 4.5, units = "in")


# MPD ~ TIME * TEMP ---------------------------------------------------------------
colnames(time_data)

# Create data, center SST
time_mpd <- time_data %>%
  filter(rarefyID %in% assemb.data$ass.name) %>% 
  group_by(rarefyID) %>% 
  mutate(sst_cent = sst - mean(sst)) %>% 
  mutate(yr_cent  = YEAR - mean(YEAR)) %>% 
  ungroup() %>% 
  mutate(sst_scaled = (sst - mean(sst_cent)) / sd(sst_cent)) %>% 
  mutate(yr_scaled  = (yr_cent - mean(yr_cent)) / sd(yr_cent))

# LMER with nonlinear effect of SST
lme_mpd_sst3 <- lmer(MPD ~ yr_scaled + sst_scaled + yr_scaled:ns(sst_scaled, df = 3) 
                     + (1 + yr_scaled | STUDY_ID) + (1 + yr_scaled | rarefyID), 
                     data = time_mpd, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(lme_mpd_sst3)

# for back-transformation
sst_mean <- mean(time_mpd$sst_cent)
sst_sd <- sd(time_mpd$sst_cent)
yr_mean <- mean(time_mpd$yr_cent)
yr_sd <- sd(time_mpd$yr_cent)

# SST values to plot
#sst_quants <- quantile(time_cti$sst_scaled, probs = c(0.01, 0.25, 0.50, 0.75, 0.99), na.rm = TRUE)
#(sst_values <- (sst_quants * sst_sd) + sst_mean)
range(time_mpd$sst)
sst_values <- c(0, 5, 10, 15, 20, 25, 30)
sst_groups <- (sst_values - sst_mean) / sst_sd
sst_labels <- setNames(paste0(round(sst_values, 0) , "°C"), sst_groups)

# year to plot 
yr_values <- c(-30, -20, -10, 0, 10, 20, 30)
yr_breaks <- (yr_values - yr_mean) / yr_sd

# Plot predicted values
mpd_sst_pred <- ggpredict(lme_mpd_sst3, terms = c("yr_scaled [all]", "sst_scaled [sst_groups]"))

plot_mpd_sst <- ggplot(data = time_mpd, aes(x = yr_scaled, y = MPD)) +
  geom_point(color = "grey50", alpha = 0.1) +
  geom_line(data = mpd_sst_pred,inherit.aes = FALSE, aes(x = x, y = predicted, color = group), 
            linewidth = 1) +  
  #geom_ribbon(data = mpd_sst_pred, inherit.aes = FALSE, aes(x = x, ymin = conf.low, ymax = conf.high, fill = group), 
  #            alpha = 0.3, color = NA) +
  scale_x_continuous(breaks = yr_breaks, labels = yr_values) +
  scale_color_manual(values = col_scale_sti(length(sst_labels)), labels = sst_labels) +
  scale_fill_manual(values = col_scale_sti(length(sst_labels)), labels = sst_labels) +
  guides(fill = "none") +
  labs(x = "Calendar year",
       y = "Assemblage MPD (My)",
       color = "Assemblage\nmean SST") + 
  theme_bw(); plot_mpd_sst

ggsave(filename = '../Graphs/MPD_time_SST.jpeg', plot_mpd_sst,
       dpi = 300, width = 5.5, height = 4.5, units = "in")

