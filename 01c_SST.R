#*##################################################################################*#
## RETRIEVE SST DATA FOR EACH GBIF OBSERVATION & CALCULATE THERMAL AFFINITY METRICS ##
#*##################################################################################*#

## Author: Philippe Fernandez-Fournier, 2025

#setwd("D:/sfuvault/SFU/MooersLab/Ch3_WinnersLosers/Analysis/")
rm(list=ls())

library(raster)
library(biooracler)
library(dplyr)
library(tidyr)
library(CoordinateCleaner)
library(countrycode)
library(ncdf4)
library(lubridate)
library(rfishbase)
library(broom)
library(ggplot2)

## Load data 
species.data <- read.csv("./data/data_species_evol.csv")

comm.data <- read.csv("./data/data_pres_evol.csv") 

assemb.data <- read.csv("./data/data_assemb.csv")
assemb.data$study <- as.factor(as.character(assemb.data$study))

time_data <- read.csv("./data/data_time_final.csv")

gbif_occdata_clean <- read.csv("./results/gbif_occurrences_clean.csv")


# ERA5 data ---------------------------------------------------------------
#https://cds.climate.copernicus.eu/

# Open the NetCDF file for SST data
nc <- nc_open("./data/ERA5/era5_sst.nc")

print(nc)

# Extract necessary variables from the NetCDF file
sst <- ncvar_get(nc, "sst")
time <- ncvar_get(nc, "valid_time")
lat <- ncvar_get(nc, "latitude")
lon <- ncvar_get(nc, "longitude")

# Convert the time to year-month format
time <- as.POSIXct(time, origin = "1970-01-01", tz = "UTC")
sst_year_month <- paste(year(time), month(time), sep = "-")



# SPECIES STI: GBIF -------------------------------------------------------------

# Function to extract SST value for given coordinates and year-month
get_monthly_sst <- function(lat_val, lon_val, year_month, search_radius = 1) {
  #lat_val <- gbif_sst$decimalLatitude[1]
  #lon_val <- gbif_sst$decimalLongitude_ERA5[1]
  #year_month <- gbif_sst$year_month[1]
  
  # Find the closest latitude and longitude index
  lat_idx <- which.min(abs(lat - lat_val))
  lon_idx <- which.min(abs(lon - lon_val))
  
  # Find the corresponding year-month
  time_idx <- which(sst_year_month == year_month)
  
  # Ensure the indices are within bounds
  lon_min <- max(lon_idx - search_radius, 1)  # Lower bound
  lon_max <- min(lon_idx + search_radius, length(lon))  # Upper bound
  lat_min <- max(lat_idx - search_radius, 1)  # Lower bound
  lat_max <- min(lat_idx + search_radius, length(lat))  # Upper bound
  
  # Extract the SST value
  sst_value <- sst[lon_idx, lat_idx, time_idx]
  
  # If SST is NaN (probably on land) extract values from the nearby grid points within bounds
  if (is.nan(sst_value)) {
    sst_value_rad <- sst[lon_min:lon_max, lat_min:lat_max, time_idx]
    
    # Compute the mean of the SST values in the search radius
    sst_value <- mean(sst_value_rad, na.rm = TRUE)
  }
  
  # Convert to Celsius
  sst_celsius <- sst_value - 273.15
  
  return(sst_celsius)
}


# Filter data
gbif_sst <- gbif_occdata_clean %>%
  filter(!is.na(month)) %>% 
  filter(basisOfRecord == "HUMAN_OBSERVATION") %>% 
  mutate(year_month = paste(year, month, sep = "-")) %>% 
  filter(year <= 2024)

# Convert longitude to degrees_east like in ERA5
# "longitude values in the range [0:360] referenced to the Greenwich Prime Meridian."
gbif_sst$decimalLongitude_ERA5 <- ifelse(gbif_sst$decimalLongitude < 0, 
                                         gbif_sst$decimalLongitude + 360, 
                                         gbif_sst$decimalLongitude) 

# get SST values for each observation
gbif_sst <- gbif_sst %>%
  rowwise() %>%
  mutate(sst_value = get_monthly_sst(decimalLatitude, decimalLongitude_ERA5, year_month))

# Summarize Species Temperature Index (STI, ie thermal affinity) for each species
# Mid-point of the 5th and 95th quantiles (Day et al. 2018)
spp_sti <- gbif_sst %>%
  group_by(sp.name) %>%
  filter(sst_value >= quantile(sst_value, 0.05, na.rm = TRUE),
         sst_value <= quantile(sst_value, 0.95, na.rm = TRUE)) %>%
  summarise(sti = median(sst_value, na.rm = TRUE),
            sti_lower = min(sst_value, na.rm = TRUE),
            sti_upper = max(sst_value, na.rm = TRUE))


write.csv(spp_sti, "./results/species_STI_ERA5.csv", row.names = FALSE)
#spp_sti <- read.csv("./results/species_STI_ERA5.csv")


# SPECIES STI: FISHBASE ---------------------------------------------------
#species.data <- read.csv("./data/data_species_final.csv") 

# Retrieve temperature data from FishBase
species_list <- gsub("_", " ", unique(species.data$sp.name))
sp_stocks_fb <- stocks(species_list = species_list)
sp_stocks <- sp_stocks_fb %>% 
  select(Species, TempPref50, TempPreferred, EnvTemp)
sp_stocks$sp.name <- gsub(" ", "_", sp_stocks$Species)

plot(sp_stocks$TempPref50, sp_stocks$TempPreferred)
abline(a = 0, b = 1, col = "red")
     
length(which(!is.na(sp_stocks$TempPref50)))
length(which(!is.na(sp_stocks$TempPreferred)))

# Compare GBIF STI vs. FishBase Temp Pref
spp_sti <- left_join(sp_stocks, spp_sti, by = join_by(sp.name))

plot(spp_sti$TempPreferred, spp_sti$sti,
     xlab = "FishBase STI", 
     ylab = "GBIF STI")
abline(a = 0, b = 1, col = "red")

# Replace GBIF STI by FishBase Temp Pref when available
spp_sti <- spp_sti %>%
  mutate(sti = if_else(!is.na(TempPreferred), TempPreferred, sti))

# Fill missing by checking valid name
(spp_miss <- unique(species.data$sp.name[!species.data$sp.name %in% spp_sti$sp.name]))

spp_sti2 <- spp_sti

for (sp in spp_miss) {
  # sp = spp_miss[2]
  cat("\r", paste0(which(spp_miss == sp), "/", length(spp_miss)))
  sp_query <- gsub("_", " ", sp)
  sp_valid <- tryCatch(validate_names(sp_query), error = function(e) NULL)
  
  if (!is.null(sp_valid) && length(sp_valid) > 0) {
    valid_name <- sp_valid[1]  # take the first matched/accepted name
    sp_fishbase <- tryCatch(stocks(species_list = valid_name), error = function(e) NULL)
    
    if (!is.null(sp_fishbase) && !is.na(sp_fishbase$TempPreferred[1])) {
      spp_sti2 <- rbind(spp_sti2, tibble(
        Species = valid_name,
        TempPref50 = NA,
        TempPreferred = sp_fishbase$TempPreferred[1],
        EnvTemp = sp_fishbase$EnvTemp[1],
        sp.name = sp,
        sti = sp_fishbase$TempPreferred[1],
        sti_lower = NA,
        sti_upper = NA
      ))
    }
  }
}

# FILL IN LEFTOVER MISSING SPECIES MANUALLY 
# From mean preferred temp on FishBase website, saved here:
spp_miss_df <- read.csv("./data/species_sti_fishbase_synonyms.csv")

spp_sti2_nona <- spp_sti2[!is.na(spp_sti2$sti),]

spp_sti3 <- bind_rows(spp_sti2_nona, spp_miss_df)
spp_sti3 <- spp_sti3[-which(duplicated(spp_sti3$sp.name)),]

all(unique(species.data$sp.name) %in% spp_sti3$sp.name) # yes!

write.csv(spp_sti3, "./results/species_STI_final.csv", row.names = FALSE)
#spp_sti3 <- read.csv("./results/species_STI_final.csv")


# MODIFY SPECIES DATA -----------------------------------------------------
#species.data <- read.csv("./data/data_species_final.csv") 

# Add STI
species.data$sti <- spp_sti3$sti[match(species.data$sp.name, spp_sti3$sp.name)]

# CHECK: Number of species without STI (should be zero)
length(unique(species.data$sp.name[which(is.na(species.data$sti))]))

# Add EnvTemp to species.data
species.data$EnvTemp <- spp_sti3$EnvTemp[match(species.data$sp.name, spp_sti3$Species)]

# Save
write.csv(species.data, "./data/data_species_final.csv", row.names = FALSE)

# MODIFY COMM DATA -----------------------------------------------------
#comm.data <- read.csv("./data/data_pres_final.csv") 

comm.data$sti <- spp_sti3$sti[match(comm.data$sp.name, spp_sti3$sp.name)]

# CHECK: Number of species without STI (should be zero)
length(which(is.na(comm.data$sti)))

write.csv(comm.data, "./data/data_pres_final.csv", row.names = FALSE)

# ADD CTI TO RAREFIED_MEDIANS -----------------------------------------------
# Community Temperature Index
#time_data <- read.csv("./data/data_time_final.csv")

# Calculate assemblage means of STI and shifts through time
assemb_cti <- comm.data %>%
  filter(presence == 1) %>% # only present species
  group_by(ass.name, year) %>%
  summarise(cti = mean(sti, na.rm = TRUE), .groups = "keep") %>% # Average STI of a given assemblage = CTI
  ungroup() %>%
  rename(rarefyID = ass.name,
         YEAR = year)

#time_data <- time_data %>%  select(-cti)

time_data2 <- left_join(time_data, assemb_cti, by = join_by(rarefyID, YEAR))

write.csv(time_data2, "./data/data_time_final.csv", row.names = FALSE)


# ASSEMB SST YEARLY -----------------------------------------------------------
#time_data <- read.csv("./data/data_time_final.csv")

# Expand data to all years between min and max of each assemblage
time_data_full <- time_data %>%
  dplyr::select(rarefyID, STUDY_ID, YEAR, lat.grid.center, lon.grid.center_ERA5) %>%
  group_by(rarefyID) %>%
  complete(YEAR = seq(min(YEAR), max(YEAR))) %>% # Generate all years from min to max per rarefyID
  mutate(STUDY_ID = first(na.omit(STUDY_ID)), # Fill other columns based on first value for this assemblage
         lat.grid.center = first(na.omit(lat.grid.center)),
         lon.grid.center_ERA5 = first(na.omit(lon.grid.center_ERA5))) %>%
  ungroup()

# get SST values for each (full) observed year
time_data_full <- time_data_full %>%
  rowwise() %>%
  mutate(sst = get_yearly_sst(lat.grid.center, lon.grid.center_ERA5, YEAR))

# Add new centered YEAR
time_data_full <- time_data_full %>%
  group_by(rarefyID) %>%
  mutate(cYEAR = YEAR - mean(YEAR)) %>% # center year (by assemblage)
  ungroup()

write.csv(time_data_full, "./data/data_time_allyears_final.csv", row.names = FALSE)

# ASSEMB SST MONTHLY -----------------------------------------------------------
#time_data <- read.csv("./data/data_time_final.csv")
#time_data_full <- read.csv("./data/data_time_allyears_final.csv")
#assemb.data <- read.csv("./data/data_assemb_final.csv")

# Expand data to all months (01-12) for each year between min and max for each assemblage
time_data_monthly <- time_data %>%
  dplyr::select(rarefyID, YEAR) %>%
  group_by(rarefyID) %>%
  complete(YEAR = seq(min(YEAR), max(YEAR))) %>%
  ungroup()

time_data_monthly <- time_data_monthly %>%
  tidyr::crossing(MONTH = 1:12) %>% 
  arrange(rarefyID, YEAR, MONTH) %>%
  mutate(year_month = paste(YEAR, MONTH, sep = "-"))


# Add columns
time_data_monthly <- left_join(time_data_monthly, 
                               time_data %>% 
                                 dplyr::select(rarefyID, STUDY_ID, lat.grid.center, lon.grid.center_ERA5) %>%
                                 distinct(), by = "rarefyID")


# Get SST values for each (full) observed year-month combination
time_data_monthly <- time_data_monthly %>%
  rowwise() %>%
  mutate(sst = get_monthly_sst(lat_val = lat.grid.center, lon_val = lon.grid.center_ERA5, year_month = year_month))

# Add new centered YEAR
time_data_monthly <- time_data_monthly %>%
  group_by(rarefyID) %>%
  mutate(cYEAR = YEAR - mean(YEAR),
         obs_month = row_number()) %>%  # Create the sequential observation month for each rarefyID
  ungroup()

# Example assemblages
set.seed(163)  
samp_assemb <- sample(time_data_monthly$rarefyID, 25) 

# Filter and plot
plot_sst_assemb <- time_data_monthly %>%
  filter(rarefyID %in% samp_assemb) %>%
  ggplot(aes(x = obs_month, y = sst)) +
  geom_line() +
  geom_smooth(method = "lm", se = FALSE, color = "royalblue3", linewidth = 1) +
  labs(title = "SST per month in 25 sampled assemblages", 
       x = "Month of observation", 
       y = "SST (°C)") +
  ylim(c(-5,30)) +
  facet_wrap(~ rarefyID, ncol = 5, scales = "free_x") +
  theme_bw(); plot_sst_assemb


ggsave("../Graphs/assemb_samp_monthlySST2.jpeg", plot_sst_assemb,
       dpi = 300, width = 8, height = 7, units = "in")



## SAVE
write.csv(time_data_monthly, "./data/data_time_allmonths_final.csv", row.names = FALSE)
#time_data_monthly <- read.csv("./data/data_time_allmonths_final.csv")

# Add assemblage temperature (mean monthly SST)
sst_mean <- time_data_monthly %>% 
  group_by(rarefyID) %>% 
  summarize(SST_mean = mean(sst)) %>% 
  rename(ass.name = rarefyID)

assemb.data <- left_join(assemb.data, sst_mean, by = join_by(ass.name))


# betaSST PER ASSEMBLAGE --------------------------------------------------

# Run individual LMs for each rarefyID and extract the results
sst_time_lm <- time_data_monthly %>%
  filter(!is.na(sst)) %>% 
  group_by(rarefyID) %>%
  do({
    # Fit the linear model for each rarefyID
    model <- lm(sst ~ obs_month, data = .)
    
    # Get the model summary and tidy the results
    tidy_model <- tidy(model)
    
    # Get R2 and RMSE
    r2 <- summary(model)$r.squared
    rmse <- sqrt(mean(residuals(model)^2))
    
    # Extract slope and intercept
    slope <- tidy_model$estimate[2] 
    intercept <- tidy_model$estimate[1]
    
    # Return the results
    tibble(
      rarefyID = unique(.$rarefyID),
      r_squared = r2,
      p_value = tidy_model$p.value[2],  # p-value of the slope
      rmse = rmse,
      slope = slope,
      intercept = intercept
    )
  }) %>%
  ungroup()

# Slope in degrees C per decade
sst_time_lm$slope_decade <- sst_time_lm$slope * 12 * 10

# MODIFY ASSEMB.DATA ------------------------------------------------------

assemb.data$betaSST_lm <- sst_time_lm$slope_decade[match(assemb.data$ass.name, sst_time_lm$rarefyID)]

write.csv(assemb.data, "./data/data_assemb_final.csv", row.names = FALSE)






