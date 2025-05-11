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
  dplyr::select(STUDY_ID, rarefyID, cYEAR, PhyloSor)

# Change PhyloSor to DISsimilarity
time_data_brm$PhyloSor <- 1-time_data_brm$PhyloSor

#hist(time_data_brm$PhyloSor)

# STUDY_ID:   study ID
# rarefyID:   assemblage ID ("study_grid")
# cYEAR:      centered YEAR
# Phylosor:   phylogenetic Sorensen similarity, is a measure of phylogenetic community similarity 

## Formula
formula.PhyloSor <- as.formula('PhyloSor ~ cYEAR + (1+cYEAR|STUDY_ID) + (1+cYEAR|rarefyID)') # same as DDI paper

## Set priors
# check which parameters can have priors
#get_prior(formula.PhyloSor, data = time_data_brm, family = brmsfamily('beta'))

# Weakly informative priors        
hier_prior_PhyloSor <- c(
  set_prior("normal(0, 1)", class = "Intercept"),                    # Intercept for mean Phylosor (logit scale)
  set_prior("normal(1, 1)", class = "b", coef = "cYEAR"),            # Slope is positive
  
  set_prior("exponential(0.5)", class = "sd", group = "STUDY_ID"),   # Study-level intercept SD
  set_prior("exponential(1)", class = "sd", coef = "cYEAR", group = "STUDY_ID"), # Study-level slope SD
  
  set_prior("exponential(0.5)", class = "sd", group = "rarefyID"),   # Assemblage-level intercept SD
  set_prior("exponential(1)", class = "sd", coef = "cYEAR", group = "rarefyID") # Assemblage-level slope SD
)

## HIERARCHICAL BRM
brm_PhyloSor <- brm(formula.PhyloSor, 
               family  = brmsfamily('zero_one_inflated_beta'), # or gaussian??
               data    = time_data_brm,
               prior   = hier_prior_PhyloSor,
               warmup  = 1000, 
               iter    = 4000, 
               init    = 0,
               control = list(adapt_delta = 0.9, 
                              max_treedepth = 15),
               cores   = ncores, 
               chains  = 4,
               refresh = 100)

save(brm_PhyloSor, file = paste0('./results/brm_PhyloSor_diss_', st, '.RData'))

