#*##############################*#
## RETRIEVE AND CLEAN GBIF DATA ##
#*##############################*#

## Author: Philippe Fernandez-Fournier, 2025

#setwd("D:/sfuvault/SFU/MooersLab/Ch3_WinnersLosers/Analysis/")
rm(list=ls())

library(dplyr)
library(ggplot2)
library(rgbif)

## Load data 
species.data <- read.csv("./data/data_species_evol.csv")

# RETRIEVE DATA FROM GBIF -----------------------------------

species_list <- unique(species.data$sp.name)

# Set up GBIF fields to query and empty dataframe
fields_gbif <- c("name", "order", "family",              # Taxonomy
                 "decimalLatitude", "decimalLongitude",  # Coordinates
                 "year", "month", "day",                 # Date
                 "countryCode", "country",               # Place
                 "basisOfRecord", "issues")              # Record details


# loop to retrieve GBIF observations and save each species as a file
for(ii in 1:length(species_list)){
  
  species_name <- gsub("_", " ", species_list[ii])
  
  # Track progress
  cat("\r", paste0(ii, '/', length(species_list), '    '))
  
  tryCatch({
    # Fetch occurrence data for the species
    occurrences <- occ_search(scientificName = species_name, 
                              limit = 10000, # list first n occurrences (ordered latest to oldest)
                              fields = fields_gbif)$data
    
    # If no data or latitude returned, skip to next
    if (is.null(occurrences) || nrow(occurrences) == 0 || !"decimalLatitude" %in% colnames(occurrences)) {
      cat("\nNo data for:", species_name, "\n")
      return(NULL)
    }
    
    # Remove duplicate observations
    occurrences <- occurrences %>%
      distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)
    
    # Add species name to occurrences
    occurrences$sp.name <- species_list[ii]
    
    # Save the occurrences data frame for this species
    write.csv(occurrences, paste0("./results/gbif/", species_list[ii], ".csv"), row.names = FALSE)
    
  }, 
  
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}


#  CHECK MISSINGNESS AND RELOOP -------------------------------------------

# list files saved and species
gbif_files <- list.files(path = "./results/gbif/", pattern = "\\.csv$", full.names = FALSE)
gbif_spp <- gsub(".csv", "", gbif_files)

# Run loop again with these species
species_mism <- species_list[!species_list %in% gbif_spp]


# CHECK TAXONOMY AND RELOOP ------------------------------------------------

# Function to get the latest accepted taxonomy
get_species_names <- function(species_list) {
  # Empty data frame for results
  species_data <- data.frame(
    sp_orig = character(),
    sp_new = character(),
    stringsAsFactors = FALSE
  )
  
  # Loop through the species list
  for (species in species_list) {
    tryCatch({
      # Query the GBIF backbone taxonomy
      response <- name_backbone(name = species)
      
      # Append result to species_data
      species_data <- rbind(
        species_data,
        data.frame(
          sp_orig = species,
          sp_new = response$species,
          stringsAsFactors = FALSE
        )
      )
    }, error = function(e) {
      # Append an entry for unresolved species
      species_data <- rbind(
        species_data,
        data.frame(
          sp_orig = species,
          sp_new = NA,
          stringsAsFactors = FALSE
        )
      )
    })
  }
  
  return(species_data)
}

sp_taxonomy <- get_species_names(species_mism)

# Run the new names in GBIF (Run querying loop above again)
species_list <- sp_taxonomy$sp_new

# revert to original names
gbif_files <- list.files(path = "./results/gbif/", pattern = "\\.csv$", full.names = FALSE)
gbif_spp <- gsub(".csv", "", gbif_files)

gbif_spp_new <- gbif_spp[gbif_spp %in% sp_taxonomy$sp_new]

# Rename files to original names by hand...
sp_taxonomy

## STILL 17 SPECIES (OUT OF 2002 TOTAL) UNMATCHED BUT NOT BAD


# COMBINE GBIF OBSERVATIONS -----------------------------------------------
gbif_files <- list.files(path = "./results/gbif/", pattern = "\\.csv$", full.names = TRUE)

gbif_occdata <- gbif_files %>%
  lapply(read.csv) %>%   # Read each CSV file into a list of data frames
  bind_rows()            # Combine all data frames into one

## Filter and clean
# Remove preserved specimen (as in Khaliq et al. 2024)
gbif_occdata <- gbif_occdata[!gbif_occdata$basisOfRecord == "PRESERVED_SPECIMEN",]

# Remove NAs
gbif_occdata <- gbif_occdata %>% 
  filter(!is.na(decimalLatitude)) 

# Save
write.csv(gbif_occdata, "./results/gbif_occurrences.csv", row.names = FALSE)
#gbif_occdata <- read.csv("./results/gbif_occurrences.csv")

# CLEAN COORDINATES -------------------------------------------------------

sp_names <- unique(gbif_occdata$sp.name)

gbif_occdata_clean <- data.frame()

for (i in 1:length(sp_names)){
  tryCatch({
    
    # Track progress
    cat("\r", paste0(i, '/', length(sp_names), '    '))
    
    spec_i <- gbif_occdata[gbif_occdata$sp.name == sp_names[i],]
    spec_i$countryCode <- countrycode(spec_i$countryCode, origin="iso2c", destination = "iso3c", 
                                      warn = FALSE) # Some countries are NA but they don't get flagged in cleaning, so it's OK
    trouble <- clean_coordinates(x = spec_i, lon = "decimalLongitude", lat = "decimalLatitude",
                                 countries = "countryCode", species = NULL,
                                 tests = c("capitals",   # tests if it is near a capital
                                           "centroids",  # tests a radius around country centroids
                                           "gbif"),      # tests a one-degree radius around the GBIF headquarters in Copenhagen, Denmark.
                                 capitals_rad = 5000,    # Radius around capital for "capitals" test 
                                 verbose = FALSE)
    
    spec_i_clean <- spec_i[trouble$.summary,] # filter to non-problematic coordinates
    
    # bind to new gbif occurrences dataset
    gbif_occdata_clean <- bind_rows(gbif_occdata_clean, spec_i_clean)
    
  },
  error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# Remove years that are NA or before 1980
gbif_occdata_clean <- gbif_occdata_clean %>%
  filter(!is.na(year)) %>%
  filter(year >= 1980)

# Save
write.csv(gbif_occdata_clean, "./results/gbif_occurrences_clean.csv", row.names = FALSE)
#gbif_occdata_clean <- read.csv("./results/gbif_occurrences_clean.csv")

# Save number of observation before and after cleaning
nobs_spp       <- table(gbif_occdata$sp.name)
nobs_spp_clean <- table(gbif_occdata_clean$sp.name)
nobs_spp_clean <- nobs_spp_clean[match(names(nobs_spp), names(nobs_spp_clean))] # reorder to match just in case

prop_cleaned_df <- data.frame(sp.name     = names(nobs_spp),
                              n_obs_orig  = as.numeric(nobs_spp),
                              n_obs_clean = as.numeric(nobs_spp_clean))

write.csv(prop_cleaned_df, "./results/gbif_nobs_cleaning.csv", row.names = FALSE)

