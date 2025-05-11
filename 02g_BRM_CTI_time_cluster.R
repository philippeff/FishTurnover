#*#########################################################*#
## Hierarchical model of community temperature index (CTI) ##
#*#########################################################*#

## Author: Philippe Fernandez-Fournier, 2025

##  This script is designed to run on a Digital Alliance Canada cluster
##  See shell script for job submission to the cluster

#setwd("D:/sfuvault/SFU/MooersLab/Ch3_WinnersLosers/Analysis")
rm(list=ls())

library(dplyr)
library(brms)

## date
st=format(Sys.time(), "%Y-%m-%d")

## Load data 
time_data <- read.csv("./data/data_time_final.csv")

# Create dataset of only columns needed
time_data_brm <- time_data %>%
  filter(!is.na(cti)) %>% # remove NAs
  dplyr::select(STUDY_ID, rarefyID, cYEAR, cti)

# STUDY_ID: study ID
# rarefyID: assemblage ID ("study_grid")
# cYEAR:    centered YEAR
# cti:     Community Temperature Index

# Formula
form_brm_cti <- as.formula("cti ~ cYEAR + (1+cYEAR|STUDY_ID) + (1+cYEAR|rarefyID)")

## Bayesian model
brm_CTI <- brm(form_brm_cti, data = time_data_brm,
               family = gaussian(),                            # Likelihood (default for continuous outcomes)
               prior = c(prior(normal(0, 1), class = "b"),     # Prior for fixed effects
                         prior(cauchy(0, 1), class = "sd")),   # Prior for random effects
               control = list(adapt_delta = 0.90,
                              max_treedepth = 15),
               iter = 4000,                                    # Number of iterations (including warmup)
               warmup = 1000,
               chains = 4,
               cores = 4,
               refresh = 100)

save(brm_CTI, file = paste0('./results/brm_CTI_', st, '.RData'))
