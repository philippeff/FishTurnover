#*####################################################*#
## Hierarchical model of JACCARD SIMILARITY THROUGH TIME ##
#*####################################################*#

## Author: Philippe Fernandez-Fournier, 2025

##  This script is designed to run on a Digital Alliance Canada cluster
##  See shell script for job submission to the cluster

#https://vuorre.com/posts/2019-02-18-analyze-analog-scale-ratings-with-zero-one-inflated-beta-models/

#setwd("D:/sfuvault/SFU/MooersLab/Ch3_WinnersLosers/Analysis")
rm(list=ls())

## packages	
library(brms)
library(dplyr)

## date
st=format(Sys.time(), "%Y-%m-%d")

# set the number of cores
ncores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
#ncores <- 4

## data
time_data <- read.csv("./data/data_time_final.csv")

# Create dataset of only columns needed
time_data_brm <- time_data %>%
  dplyr::select(STUDY_ID, rarefyID, cYEAR, jaccard_turnover)

#hist(time_data_brm$jaccard_turnover)

# STUDY_ID:   study ID
# rarefyID:   assemblage ID ("study_grid")
# cYEAR:      centered YEAR
# jaccard_turnover: Turnover component of Jaccard dissimilarity to the first year of observation within assemblage

## Formula
formula.Jacc <- as.formula('jaccard_turnover ~ cYEAR + (1+cYEAR|STUDY_ID) + (1+cYEAR|rarefyID)') # same as DDI paper

## Set priors
# check which parameters can have priors
#get_prior(formula.Jacc, data = time_data_brm, family = brmsfamily('beta'))

# Weakly informative priors        
hier_prior_Jacc <- c(
  set_prior("normal(0, 1)", class = "Intercept"),                    # Intercept for mean Jaccard
  set_prior("normal(1, 1)", class = "b", coef = "cYEAR"),            # Slope is positive
  
  set_prior("exponential(0.5)", class = "sd", group = "STUDY_ID"),   # Study-level intercept SD
  set_prior("exponential(1)", class = "sd", coef = "cYEAR", group = "STUDY_ID"), # Study-level slope SD
  
  set_prior("exponential(0.5)", class = "sd", group = "rarefyID"),   # Assemblage-level intercept SD
  set_prior("exponential(1)", class = "sd", coef = "cYEAR", group = "rarefyID") # Assemblage-level slope SD
)

## HIERARCHICAL BRM
brm_Jacc <- brm(formula.Jacc, 
               family  = brmsfamily('zero_one_inflated_beta'), # or gaussian??
               data    = time_data_brm,
               prior   = hier_prior_Jacc,
               warmup  = 1000, 
               iter    = 4000, 
               init    = 0,
               control = list(adapt_delta = 0.9, 
                              max_treedepth = 15),
               cores   = ncores, 
               chains  = 4,
               refresh = 100)

save(brm_Jacc, file = paste0('./results/brm_Jacc_turnover_', st, '.RData'))

