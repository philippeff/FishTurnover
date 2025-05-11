#*############################################*#
## Hierarchical model of SPECIES RICHNESS(SR) ##
#*############################################*#

## Author: Philippe Fernandez-Fournier, 2025

##  This script is designed to run on a Digital Alliance Canada cluster
##  See shell script for job submission to the cluster

setwd("D:/sfuvault/SFU/MooersLab/Ch3_WinnersLosers/Analysis")
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
  dplyr::select(STUDY_ID, rarefyID, cYEAR, S)

# STUDY_ID: study ID
# rarefyID: assemblage ID ("study_grid")
# cYEAR:    centered YEAR
# S:     Species richness (SR)


## Formula
formula.SR <- as.formula('S ~ cYEAR + (1+cYEAR|STUDY_ID) + (1+cYEAR|rarefyID)') # same as DDI paper

## Set priors
# check which parameters can have priors
#get_prior(formula.SR, data = time_data_brm, family = brmsfamily('gaussian'))

# Weakly informative priors        
hier_prior_SR <- c(
  set_prior(prior = 'exponential(0.5)', class = 'sd', coef = 'Intercept', group = 'rarefyID'), # sd for assemblage level int offsets
  set_prior(prior = 'exponential(1)',   class = 'sd', coef = 'cYEAR',     group = 'rarefyID'), # sd for assemblage level slope offsets
  
  set_prior(prior = 'exponential(0.5)', class = 'sd', coef = 'Intercept', group = 'STUDY_ID'), # sd for study level int offsets
  set_prior(prior = 'exponential(1)',   class = 'sd', coef = 'cYEAR',     group = 'STUDY_ID'), # sd for study level slope offsets
  
  set_prior(prior = 'exponential(0.25)', class = 'sigma'), # Residual sd
  
  set_prior(prior = 'lkj(2)', class = 'cor')) # parameter for covariance matrix


## HIERARCHICAL BRM
# When choosing gaussian vs. lognormal family:
# in the log-transformed model, the slope is relative (percentage change).
# in the untransformed model, the slope is absolute (direct change in value).

brm_SR <- brm(formula.SR, 
               family = brmsfamily('lognormal'), # or gaussian??
               data =   time_data_brm,
               prior =  hier_prior_SR,
               warmup = 1000, 
               iter = 4000, 
               init = 'random',
               control = list(adapt_delta = 0.9, 
                              max_treedepth = 15),
               cores = ncores, 
               chains = 4,
               refresh = 100)

save(brm_SR, file = paste0('./results/brm_SR_', st, '.RData'))

