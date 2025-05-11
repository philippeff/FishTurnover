## SUMMARIZE POSTERIOR FROM HIERARCHICAL BAYESIAN MODELS ##

## Author: Philippe Fernandez-Fournier, 2025

#setwd("D:/sfuvault/SFU/MooersLab/Ch3_WinnersLosers/Analysis")
rm(list=ls())

## Load packages
library(brms)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(bayesplot)
library(lme4)

## Load data
assemb.data <- read.csv("./data/data_assemb_final.csv")
assemb.data$study <- as.character(assemb.data$study)

time_data <- read.csv("./data/data_time_final.csv")

## Empty posterior table
post_table <- data.frame()
#post_table <- read.csv("../Graphs/posteriors_table.csv")

# BRM of MPD ~ TIME ----------------------------------------------------------
#rm(list=ls())

# load latest MPD BRM model (by date)
MPD_files <- list.files(pattern = 'brm_MPD', path = "./results/")
load(paste0('./results/', MPD_files[length(MPD_files)]))

brm_MPD$formula
brm_MPD$family

# Extract the posterior samples
post_MPD <- as_draws_df(brm_MPD)

# Extract the fixed effect for cYEAR
post_b_cYEAR_MPD <- post_MPD$b_cYEAR
hist(post_b_cYEAR_MPD, breaks = 100)
median(post_b_cYEAR_MPD)
quantile(post_b_cYEAR_MPD, probs = c(0.025, 0.975))

median(post_MPD$b_Intercept)

# Proportion of distribution in the same direction (sign) as the median
MPD_pdir <- length(which(sign(post_b_cYEAR_MPD) == sign(median(post_b_cYEAR_MPD)))) / length(post_b_cYEAR_MPD)

## Add to posteriors table
MPD_summ <- summary(brm_MPD)$fixed
post_table <- rbind(post_table, data.frame(metric = "MPD",
                                           b_median = median(post_b_cYEAR_MPD),
                                           Int_median = median(post_MPD$b_Intercept),
                                           CI_low = MPD_summ[2,"l-95% CI"],
                                           CI_high = MPD_summ[2,"u-95% CI"], 
                                           pdir = MPD_pdir, 
                                           Rhat = MPD_summ[2,"Rhat"], 
                                           ESS = MPD_summ[2,"Bulk_ESS"]))

# Plot b_cYEAR
# Convert to array and subset for plotting
post_MPD_array <- as.array(brm_MPD)
thin_factor <- 10 # How much thinning for subset? use 100 for testing
post_MPD_array_sub <- post_MPD_array[seq(1, dim(post_MPD_array)[1], thin_factor), , ]

color_scheme_set("darkgray")
plot_b_MPD <- mcmc_areas(post_MPD_array_sub, pars = c("b_cYEAR"), prob = 0.95) +
  geom_vline(aes(xintercept = 0), color = "grey20", linetype = "dashed") + 
  labs(title = "posterior distribution of change in MPD",
       x = "b_cYEAR (Overall slope of MPD ~ time)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  coord_cartesian(ylim = c(0.9, 0.95)); plot_b_MPD

ggsave(filename = '../Graphs/MPD_b_cYEAR.jpeg', plot_b_MPD,
       dpi = 300, width = 6, height = 5, units = "in")

# Extract the random slopes (offsets from b_cYEAR) at the rarefyID level
post_r_assemb_MPD <- post_MPD %>% 
  dplyr::select(starts_with("r_rarefyID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_assemb_MPD[1:10, 1:4]

# Extract the random slopes (offsets from b_cYEAR) at the STUDY_ID level
post_r_study_MPD <- post_MPD %>% 
  dplyr::select(starts_with("r_STUDY_ID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_study_MPD[1:10, 1:4]

# Combine the effects to calculate the posterior slopes for each rarefyID
assemb_id_MPD <- colnames(post_r_assemb_MPD) %>% gsub("r_rarefyID\\[|,cYEAR\\]", "", .)
study_id_MPD <- sub("_.*", "", assemb_id_MPD)
study_cols_MPD <- paste0("r_STUDY_ID[", study_id_MPD, ",cYEAR]")

post_b_assemb_MPD <- post_r_assemb_MPD + post_r_study_MPD[study_cols_MPD] + post_b_cYEAR_MPD

# Summarize the posterior distributions for each rarefyID
slopes_summary_MPD <- apply(post_b_assemb_MPD, 2, function(x) {
  c(median = median(x), sd = sd(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
})

slopes_summary_df_MPD <- as.data.frame(t(slopes_summary_MPD))
slopes_summary_df_MPD$rarefyID <- sub(".*\\[([0-9_]+),.*", "\\1", rownames(slopes_summary_df_MPD))

head(slopes_summary_df_MPD)

## Add to assemb.data
assemb.data$median_betaMPD <- slopes_summary_df_MPD$median[match(assemb.data$ass.name, slopes_summary_df_MPD$rarefyID)]

# Proportion of assemblages with a positive slope
length(which(assemb.data$median_betaMPD > 0))
nrow(assemb.data)
length(which(assemb.data$median_betaMPD > 0)) / nrow(assemb.data)

write.csv(slopes_summary_df_MPD, "./results/slopes_summary_MPD.csv", row.names = FALSE)

rm(list = ls(pattern = "MPD"))
gc()

# BRM of sesMPD ~ TIME ----------------------------------------------------------

# load latest MPD BRM model (by date)
sesMPD_files <- list.files(pattern = 'brm_sesMPD', path = "./results/")
load(paste0('./results/', sesMPD_files[length(sesMPD_files)]))

brm_sesMPD$formula
brm_sesMPD$family

# Extract the posterior samples
post_sesMPD <- as_draws_df(brm_sesMPD)

# Extract the fixed effect for cYEAR
post_b_cYEAR_sesMPD <- post_sesMPD$b_cYEAR
hist(post_b_cYEAR_sesMPD, breaks = 100)
median(post_b_cYEAR_sesMPD)
quantile(post_b_cYEAR_sesMPD, probs = c(0.025, 0.975))

# Proportion of distribution in the same direction (sign) as the median
sesMPD_pdir <- length(which(sign(post_b_cYEAR_sesMPD) == sign(median(post_b_cYEAR_sesMPD)))) / length(post_b_cYEAR_sesMPD)

## Add to posteriors table
sesMPD_summ <- summary(brm_sesMPD)$fixed
post_table <- rbind(post_table, data.frame(metric = "sesMPD",
                                           b_median = median(post_b_cYEAR_sesMPD),
                                           Int_median = median(post_sesMPD$b_Intercept),
                                           CI_low = sesMPD_summ[2,"l-95% CI"],
                                           CI_high = sesMPD_summ[2,"u-95% CI"], 
                                           pdir = sesMPD_pdir, 
                                           Rhat = sesMPD_summ[2,"Rhat"], 
                                           ESS = sesMPD_summ[2,"Bulk_ESS"]))

# Plot b_cYEAR
# Convert to array and subset for plotting
post_sesMPD_array <- as.array(brm_sesMPD)
thin_factor <- 10 # How much thinning for subset? use 100 for testing
post_sesMPD_array_sub <- post_sesMPD_array[seq(1, dim(post_sesMPD_array)[1], thin_factor), , ]

color_scheme_set("darkgray")
plot_b_sesMPD <- mcmc_areas(post_sesMPD_array_sub, pars = c("b_cYEAR"), prob = 0.95) +
  geom_vline(aes(xintercept = 0), color = "grey20", linetype = "dashed") + 
  labs(title = "posterior distribution of change in sesMPD",
       x = "b_cYEAR (Overall slope of sesMPD ~ time)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  coord_cartesian(ylim = c(0.9, 0.95)); plot_b_sesMPD

ggsave(filename = '../Graphs/sesMPD_b_cYEAR.jpeg', plot_b_sesMPD,
       dpi = 300, width = 6, height = 5, units = "in")

# Extract the random slopes (offsets from b_cYEAR) at the rarefyID level
post_r_assemb_sesMPD <- post_sesMPD %>% 
  dplyr::select(starts_with("r_rarefyID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_assemb_sesMPD[1:10, 1:4]

# Extract the random slopes (offsets from b_cYEAR) at the STUDY_ID level
post_r_study_sesMPD <- post_sesMPD %>% 
  dplyr::select(starts_with("r_STUDY_ID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_study_sesMPD[1:10, 1:4]

# Combine the effects to calculate the posterior slopes for each rarefyID
assemb_id_sesMPD <- colnames(post_r_assemb_sesMPD) %>% gsub("r_rarefyID\\[|,cYEAR\\]", "", .)
study_id_sesMPD <- sub("_.*", "", assemb_id_sesMPD)
study_cols_sesMPD <- paste0("r_STUDY_ID[", study_id_sesMPD, ",cYEAR]")

post_b_assemb_sesMPD <- post_r_assemb_sesMPD + post_r_study_sesMPD[study_cols_sesMPD] + post_b_cYEAR_sesMPD

# Summarize the posterior distributions for each rarefyID
slopes_summary_sesMPD <- apply(post_b_assemb_sesMPD, 2, function(x) {
  c(median = median(x), sd = sd(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
})

slopes_summary_df_sesMPD <- as.data.frame(t(slopes_summary_sesMPD))
slopes_summary_df_sesMPD$rarefyID <- sub(".*\\[([0-9_]+),.*", "\\1", rownames(slopes_summary_df_sesMPD))

head(slopes_summary_df_sesMPD)

## Add to assemb.data
assemb.data$median_betaSESMPD <- slopes_summary_df_sesMPD$median[match(assemb.data$ass.name, slopes_summary_df_sesMPD$rarefyID)]

# Proportion of assemblages with a positive slope
length(which(assemb.data$median_betaSESMPD > 0))
nrow(assemb.data)
length(which(assemb.data$median_betaSESMPD > 0)) / nrow(assemb.data)

write.csv(slopes_summary_df_sesMPD, "./results/slopes_summary_sesMPD.csv", row.names = FALSE)

rm(list = ls(pattern = "sesMPD"))
gc()

# BRM of MNTD ~ TIME ----------------------------------------------------------

# load latest MNTD BRM model (by date)
MNTD_files <- list.files(pattern = 'brm_MNTD', path = "./results/")
load(paste0('./results/', MNTD_files[length(MNTD_files)]))

brm_MNTD$formula
brm_MNTD$family

# Extract the posterior samples
post_MNTD <- as_draws_df(brm_MNTD)

# Extract the fixed effect for cYEAR
post_b_cYEAR_MNTD <- post_MNTD$b_cYEAR
hist(post_b_cYEAR_MNTD, breaks = 100)
median(post_b_cYEAR_MNTD)
quantile(post_b_cYEAR_MNTD, probs = c(0.025, 0.975))

# Proportion of distribution in the same direction (sign) as the median
MNTD_pdir <- length(which(sign(post_b_cYEAR_MNTD) == sign(median(post_b_cYEAR_MNTD)))) / length(post_b_cYEAR_MNTD)

## Add to posteriors table
MNTD_summ <- summary(brm_MNTD)$fixed
post_table <- rbind(post_table, data.frame(metric = "MNTD",
                                           b_median = median(post_b_cYEAR_MNTD),
                                           Int_median = median(post_MNTD$b_Intercept),
                                           CI_low = MNTD_summ[2,"l-95% CI"],
                                           CI_high = MNTD_summ[2,"u-95% CI"], 
                                           pdir = MNTD_pdir, 
                                           Rhat = MNTD_summ[2,"Rhat"], 
                                           ESS = MNTD_summ[2,"Bulk_ESS"]))

# Plot b_cYEAR
# Convert to array and subset for plotting
post_MNTD_array <- as.array(brm_MNTD)
thin_factor <- 10 # How much thinning for subset? use 100 for testing
post_MNTD_array_sub <- post_MNTD_array[seq(1, dim(post_MNTD_array)[1], thin_factor), , ]

color_scheme_set("darkgray")
plot_b_MNTD <- mcmc_areas(post_MNTD_array_sub, pars = c("b_cYEAR"), prob = 0.95) +
  geom_vline(aes(xintercept = 0), color = "grey20", linetype = "dashed") + 
  labs(title = "posterior distribution of change in MNTD",
       x = "b_cYEAR (Overall slope of MNTD ~ time)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  coord_cartesian(ylim = c(0.9, 0.95)); plot_b_MNTD

ggsave(filename = '../Graphs/MNTD_b_cYEAR.jpeg', plot_b_MNTD,
       dpi = 300, width = 6, height = 5, units = "in")

# Extract the random slopes (offsets from b_cYEAR) at the rarefyID level
post_r_assemb_MNTD <- post_MNTD %>% 
  dplyr::select(starts_with("r_rarefyID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_assemb_MNTD[1:10, 1:4]

# Extract the random slopes (offsets from b_cYEAR) at the STUDY_ID level
post_r_study_MNTD <- post_MNTD %>% 
  dplyr::select(starts_with("r_STUDY_ID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_study_MNTD[1:10, 1:4]

# Combine the effects to calculate the posterior slopes for each rarefyID
assemb_id_MNTD <- colnames(post_r_assemb_MNTD) %>% gsub("r_rarefyID\\[|,cYEAR\\]", "", .)
study_id_MNTD <- sub("_.*", "", assemb_id_MNTD)
study_cols_MNTD <- paste0("r_STUDY_ID[", study_id_MNTD, ",cYEAR]")

post_b_assemb_MNTD <- post_r_assemb_MNTD + post_r_study_MNTD[study_cols_MNTD] + post_b_cYEAR_MNTD

# Summarize the posterior distributions for each rarefyID
slopes_summary_MNTD <- apply(post_b_assemb_MNTD, 2, function(x) {
  c(median = median(x), sd = sd(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
})

slopes_summary_df_MNTD <- as.data.frame(t(slopes_summary_MNTD))
slopes_summary_df_MNTD$rarefyID <- sub(".*\\[([0-9_]+),.*", "\\1", rownames(slopes_summary_df_MNTD))

head(slopes_summary_df_MNTD)

## Add to assemb.data
assemb.data$median_betaMNTD <- slopes_summary_df_MNTD$median[match(assemb.data$ass.name, slopes_summary_df_MNTD$rarefyID)]

# Proportion of assemblages with a positive slope
length(which(assemb.data$median_betaMNTD > 0))
nrow(assemb.data)
length(which(assemb.data$median_betaMNTD > 0)) / nrow(assemb.data)

write.csv(slopes_summary_df_MNTD, "./results/slopes_summary_MNTD.csv", row.names = FALSE)

rm(list = ls(pattern = "MNTD"))
gc()

# BRM of SR ~ TIME ----------------------------------------------------------

# load latest SR BRM model (by date)
SR_files <- list.files(pattern = 'brm_SR', path = "./results/")
load(paste0('./results/', SR_files[length(SR_files)]))

brm_SR$formula
brm_SR$family

# Extract the posterior samples
post_SR <- as_draws_df(brm_SR)

# Extract the fixed effect for cYEAR
post_b_cYEAR_SR <- post_SR$b_cYEAR
hist(post_b_cYEAR_SR, breaks = 100)
median(post_b_cYEAR_SR)
quantile(post_b_cYEAR_SR, probs = c(0.025, 0.975))

# Proportion of distribution in the same direction (sign) as the median
SR_pdir <- length(which(sign(post_b_cYEAR_SR) == sign(median(post_b_cYEAR_SR)))) / length(post_b_cYEAR_SR)

## Add to posteriors table
SR_summ <- summary(brm_SR)$fixed
post_table <- rbind(post_table, data.frame(metric = "SR",
                                           b_median = median(post_b_cYEAR_SR),
                                           Int_median = median(post_SR$b_Intercept),
                                           CI_low = SR_summ[2,"l-95% CI"],
                                           CI_high = SR_summ[2,"u-95% CI"], 
                                           pdir = SR_pdir, 
                                           Rhat = SR_summ[2,"Rhat"], 
                                           ESS = SR_summ[2,"Bulk_ESS"]))

# Plot b_cYEAR
# Convert to array and subset for plotting
post_SR_array <- as.array(brm_SR)
thin_factor <- 10 # How much thinning for subset? use 100 for testing
post_SR_array_sub <- post_SR_array[seq(1, dim(post_SR_array)[1], thin_factor), , ]

color_scheme_set("darkgray")
plot_b_SR <- mcmc_areas(post_SR_array_sub, pars = c("b_cYEAR"), prob = 0.95) +
  geom_vline(aes(xintercept = 0), color = "grey20", linetype = "dashed") + 
  labs(title = "posterior distribution of change in SR",
       x = "b_cYEAR (Overall slope of SR ~ time)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  coord_cartesian(ylim = c(0.9, 0.95)); plot_b_SR

ggsave(filename = '../Graphs/SR_b_cYEAR.jpeg', plot_b_SR,
       dpi = 300, width = 6, height = 5, units = "in")

# Extract the random slopes (offsets from b_cYEAR) at the rarefyID level
post_r_assemb_SR <- post_SR %>% 
  dplyr::select(starts_with("r_rarefyID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_assemb_SR[1:10, 1:4]

# Extract the random slopes (offsets from b_cYEAR) at the STUDY_ID level
post_r_study_SR <- post_SR %>% 
  dplyr::select(starts_with("r_STUDY_ID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_study_SR[1:10, 1:4]

# Combine the effects to calculate the posterior slopes for each rarefyID
assemb_id_SR <- colnames(post_r_assemb_SR) %>% gsub("r_rarefyID\\[|,cYEAR\\]", "", .)
study_id_SR <- sub("_.*", "", assemb_id_SR)
study_cols_SR <- paste0("r_STUDY_ID[", study_id_SR, ",cYEAR]")

post_b_assemb_SR <- post_r_assemb_SR + post_r_study_SR[study_cols_SR] + post_b_cYEAR_SR

# Summarize the posterior distributions for each rarefyID
slopes_summary_SR <- apply(post_b_assemb_SR, 2, function(x) {
  c(median = median(x), sd = sd(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
})

slopes_summary_df_SR <- as.data.frame(t(slopes_summary_SR))
slopes_summary_df_SR$rarefyID <- sub(".*\\[([0-9_]+),.*", "\\1", rownames(slopes_summary_df_SR))

head(slopes_summary_df_SR)

## Add to assemb.data
assemb.data$median_betaSR <- slopes_summary_df_SR$median[match(assemb.data$ass.name, slopes_summary_df_SR$rarefyID)]

# Proportion of assemblages with a positive slope
length(which(assemb.data$median_betaSR > 0))
nrow(assemb.data)
length(which(assemb.data$median_betaSR > 0)) / nrow(assemb.data)

write.csv(slopes_summary_df_SR, "./results/slopes_summary_SR.csv", row.names = FALSE)

rm(list = ls(pattern = "SR"))
gc()

# BRM of Jacc ~ TIME ----------------------------------------------------------
#rm(list=ls())

# load latest Jacc BRM model (by date)
Jacc_files <- list.files(pattern = 'brm_Jacc', path = "./results/")
load(paste0('./results/', Jacc_files[length(Jacc_files)]))

brm_Jacc$formula
brm_Jacc$family

# Extract the posterior samples
post_Jacc <- as_draws_df(brm_Jacc)

# Extract the fixed effect for cYEAR
post_b_cYEAR_Jacc <- post_Jacc$b_cYEAR
hist(post_b_cYEAR_Jacc, breaks = 100)
median(post_b_cYEAR_Jacc)
quantile(post_b_cYEAR_Jacc, probs = c(0.025, 0.975))

# Proportion of distribution in the same direction (sign) as the median
Jacc_pdir <- length(which(sign(post_b_cYEAR_Jacc) == sign(median(post_b_cYEAR_Jacc)))) / length(post_b_cYEAR_Jacc)

## Add to posteriors table
Jacc_summ <- summary(brm_Jacc)$fixed
post_table <- rbind(post_table, data.frame(metric = "Jacc",
                                           b_median = median(post_b_cYEAR_Jacc),
                                           Int_median = median(post_Jacc$b_Intercept),
                                           CI_low = Jacc_summ[2,"l-95% CI"],
                                           CI_high = Jacc_summ[2,"u-95% CI"], 
                                           pdir = Jacc_pdir, 
                                           Rhat = Jacc_summ[2,"Rhat"], 
                                           ESS = Jacc_summ[2,"Bulk_ESS"]))

# Plot b_cYEAR
# Convert to array and subset for plotting
post_Jacc_array <- as.array(brm_Jacc)
thin_factor <- 10 # How much thinning for subset? use 100 for testing
post_Jacc_array_sub <- post_Jacc_array[seq(1, dim(post_Jacc_array)[1], thin_factor), , ]

color_scheme_set("darkgray")
plot_b_Jacc <- mcmc_areas(post_Jacc_array_sub, pars = c("b_cYEAR"), prob = 0.95) +
  geom_vline(aes(xintercept = 0), color = "grey20", linetype = "dashed") + 
  labs(title = "posterior distribution of change in Jacc",
       x = "b_cYEAR (Overall slope of Jacc ~ time)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  coord_cartesian(ylim = c(0.9, 0.95)); plot_b_Jacc

ggsave(filename = '../Graphs/Jacc_b_cYEAR.jpeg', plot_b_Jacc,
       dpi = 300, width = 6, height = 5, units = "in")

# Extract the random slopes (offsets from b_cYEAR) at the rarefyID level
post_r_assemb_Jacc <- post_Jacc %>% 
  dplyr::select(starts_with("r_rarefyID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_assemb_Jacc[1:10, 1:4]

# Extract the random slopes (offsets from b_cYEAR) at the STUDY_ID level
post_r_study_Jacc <- post_Jacc %>% 
  dplyr::select(starts_with("r_STUDY_ID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_study_Jacc[1:10, 1:4]

# Combine the effects to calculate the posterior slopes for each rarefyID
assemb_id_Jacc <- colnames(post_r_assemb_Jacc) %>% gsub("r_rarefyID\\[|,cYEAR\\]", "", .)
study_id_Jacc <- sub("_.*", "", assemb_id_Jacc)
study_cols_Jacc <- paste0("r_STUDY_ID[", study_id_Jacc, ",cYEAR]")

post_b_assemb_Jacc <- post_r_assemb_Jacc + post_r_study_Jacc[study_cols_Jacc] + post_b_cYEAR_Jacc

# Summarize the posterior distributions for each rarefyID
slopes_summary_Jacc <- apply(post_b_assemb_Jacc, 2, function(x) {
  c(median = median(x), sd = sd(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
})

slopes_summary_df_Jacc <- as.data.frame(t(slopes_summary_Jacc))
slopes_summary_df_Jacc$rarefyID <- sub(".*\\[([0-9_]+),.*", "\\1", rownames(slopes_summary_df_Jacc))

head(slopes_summary_df_Jacc)

## Add to assemb.data
assemb.data$median_betaJacc <- slopes_summary_df_Jacc$median[match(assemb.data$ass.name, slopes_summary_df_Jacc$rarefyID)]

# Proportion of assemblages with a positive slope
length(which(assemb.data$median_betaJacc > 0))
nrow(assemb.data)
length(which(assemb.data$median_betaJacc > 0)) / nrow(assemb.data)

write.csv(slopes_summary_df_Jacc, "./results/slopes_summary_Jacc.csv", row.names = FALSE)

rm(list = ls(pattern = "Jacc"))
gc()

# BRM of PhyloSor ~ TIME ----------------------------------------------------------
#rm(list=ls())

# load latest PhyloSor BRM model (by date)
PhyloSor_files <- list.files(pattern = 'brm_PhyloSor', path = "./results/")
load(paste0('./results/', PhyloSor_files[length(PhyloSor_files)]))

brm_PhyloSor$formula
brm_PhyloSor$family

# Extract the posterior samples
post_PhyloSor <- as_draws_df(brm_PhyloSor)

# Extract the fixed effect for cYEAR
post_b_cYEAR_PhyloSor <- post_PhyloSor$b_cYEAR
hist(post_b_cYEAR_PhyloSor, breaks = 100)
median(post_b_cYEAR_PhyloSor)
quantile(post_b_cYEAR_PhyloSor, probs = c(0.025, 0.975))

# Proportion of distribution in the same direction (sign) as the median
PhyloSor_pdir <- length(which(sign(post_b_cYEAR_PhyloSor) == sign(median(post_b_cYEAR_PhyloSor)))) / length(post_b_cYEAR_PhyloSor)

## Add to posteriors table
PhyloSor_summ <- summary(brm_PhyloSor)$fixed
post_table <- rbind(post_table, data.frame(metric = "PhyloSor",
                                           b_median = median(post_b_cYEAR_PhyloSor),
                                           Int_median = median(post_PhyloSor$b_Intercept),
                                           CI_low = PhyloSor_summ[2,"l-95% CI"],
                                           CI_high = PhyloSor_summ[2,"u-95% CI"], 
                                           pdir = PhyloSor_pdir, 
                                           Rhat = PhyloSor_summ[2,"Rhat"], 
                                           ESS = PhyloSor_summ[2,"Bulk_ESS"]))

# Plot b_cYEAR
# Convert to array and subset for plotting
post_PhyloSor_array <- as.array(brm_PhyloSor)
thin_factor <- 10 # How much thinning for subset?
post_PhyloSor_array_sub <- post_PhyloSor_array[seq(1, dim(post_PhyloSor_array)[1], thin_factor), , ]

color_scheme_set("darkgray")
plot_b_PhyloSor <- mcmc_areas(post_PhyloSor_array_sub, pars = c("b_cYEAR"), prob = 0.95) +
  geom_vline(aes(xintercept = 0), color = "grey20", linetype = "dashed") + 
  labs(title = "posterior distribution of change in PhyloSor",
       x = "b_cYEAR (Overall slope of PhyloSor ~ time)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  coord_cartesian(ylim = c(0.9, 0.95)); plot_b_PhyloSor

ggsave(filename = '../Graphs/PhyloSor_b_cYEAR.jpeg', plot_b_PhyloSor,
       dpi = 300, width = 6, height = 5, units = "in")

# Extract the random slopes (offsets from b_cYEAR) at the rarefyID level
post_r_assemb_PhyloSor <- post_PhyloSor %>% 
  dplyr::select(starts_with("r_rarefyID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_assemb_PhyloSor[1:10, 1:4]

# Extract the random slopes (offsets from b_cYEAR) at the STUDY_ID level
post_r_study_PhyloSor <- post_PhyloSor %>% 
  dplyr::select(starts_with("r_STUDY_ID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_study_PhyloSor[1:10, 1:4]

# Combine the effects to calculate the posterior slopes for each rarefyID
assemb_id_PhyloSor <- colnames(post_r_assemb_PhyloSor) %>% gsub("r_rarefyID\\[|,cYEAR\\]", "", .)
study_id_PhyloSor <- sub("_.*", "", assemb_id_PhyloSor)
study_cols_PhyloSor <- paste0("r_STUDY_ID[", study_id_PhyloSor, ",cYEAR]")

post_b_assemb_PhyloSor <- post_r_assemb_PhyloSor + post_r_study_PhyloSor[study_cols_PhyloSor] + post_b_cYEAR_PhyloSor

# Summarize the posterior distributions for each rarefyID
slopes_summary_PhyloSor <- apply(post_b_assemb_PhyloSor, 2, function(x) {
  c(median = median(x), sd = sd(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
})

slopes_summary_df_PhyloSor <- as.data.frame(t(slopes_summary_PhyloSor))
slopes_summary_df_PhyloSor$rarefyID <- sub(".*\\[([0-9_]+),.*", "\\1", rownames(slopes_summary_df_PhyloSor))

head(slopes_summary_df_PhyloSor)

## Add to assemb.data
assemb.data$median_betaPhyloSor <- slopes_summary_df_PhyloSor$median[match(assemb.data$ass.name, slopes_summary_df_PhyloSor$rarefyID)]

# Proportion of assemblages with a positive slope
length(which(assemb.data$median_betaPhyloSor > 0))
nrow(assemb.data)
length(which(assemb.data$median_betaPhyloSor > 0)) / nrow(assemb.data)

write.csv(slopes_summary_df_PhyloSor, "./results/slopes_summary_PhyloSor.csv", row.names = FALSE)

rm(list = ls(pattern = "PhyloSor"))
gc()



# BRM of CTI ~ TIME ------------------------------------------------------

# load latest CTI BRM model (by date)
CTI_files <- list.files(pattern = 'brm_CTI', path = "./results/")
load(paste0('./results/', CTI_files[length(CTI_files)]))

brm_CTI$formula
brm_CTI$family

# Extract the posterior samples
post_CTI <- as_draws_df(brm_CTI)

# Extract the fixed effect for cYEAR
post_b_cYEAR_CTI <- post_CTI$b_cYEAR
hist(post_b_cYEAR_CTI, breaks = 100)
median(post_b_cYEAR_CTI)
quantile(post_b_cYEAR_CTI, probs = c(0.025, 0.975))

median(post_CTI$b_Intercept)

# Proportion of distribution in the same direction (sign) as the median
CTI_pdir <- length(which(sign(post_b_cYEAR_CTI) == sign(median(post_b_cYEAR_CTI)))) / length(post_b_cYEAR_CTI)

## Add to posteriors table
CTI_summ <- summary(brm_CTI)$fixed
post_table <- rbind(post_table, data.frame(metric = "CTI",
                                           b_median = median(post_b_cYEAR_CTI),
                                           Int_median = median(post_CTI$b_Intercept),
                                           CI_low = CTI_summ[2,"l-95% CI"],
                                           CI_high = CTI_summ[2,"u-95% CI"], 
                                           pdir = CTI_pdir, 
                                           Rhat = CTI_summ[2,"Rhat"], 
                                           ESS = CTI_summ[2,"Bulk_ESS"]))


# Plot b_cYEAR
# Convert to array and subset for plotting
post_CTI_array <- as.array(brm_CTI)
thin_factor <- 10 # How much thinning for subset? use 100 for testing
post_CTI_array_sub <- post_CTI_array[seq(1, dim(post_CTI_array)[1], thin_factor), , ]

color_scheme_set("darkgray")
plot_b_CTI <- mcmc_areas(post_CTI_array_sub, pars = c("b_cYEAR"), prob = 0.95) +
  geom_vline(aes(xintercept = 0), color = "grey20", linetype = "dashed") + 
  labs(title = "posterior distribution of change in CTI",
       x = "b_cYEAR (Overall slope of CTI ~ time)") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  coord_cartesian(ylim = c(0.9, 0.95)); plot_b_CTI

ggsave(filename = "../Graphs/CTI_b_cYEAR.jpeg", plot_b_CTI,
       dpi = 300, width = 6, height = 5, units = "in")


# Extract the random slopes (offsets from b_cYEAR) at the rarefyID level
post_r_assemb_CTI <- post_CTI %>% 
  dplyr::select(starts_with("r_rarefyID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_assemb_CTI[1:10, 1:4]

# Extract the random slopes (offsets from b_cYEAR) at the STUDY_ID level
post_r_study_CTI <- post_CTI %>% 
  dplyr::select(starts_with("r_STUDY_ID[")) %>%
  dplyr::select(ends_with("cYEAR]"))
post_r_study_CTI[1:10, 1:4]

# Combine the effects to calculate the posterior slopes for each ass.name
assemb_id_CTI <- colnames(post_r_assemb_CTI) %>% gsub("r_rarefyID\\[|,cYEAR\\]", "", .)
study_id_CTI <- sub("_.*", "", assemb_id_CTI)
study_cols_CTI <- paste0("r_STUDY_ID[", study_id_CTI, ",cYEAR]")

post_b_assemb_CTI <- post_r_assemb_CTI + post_b_cYEAR_CTI + post_r_study_CTI[study_cols_CTI] # assemblage offsets + Overall slope + study offsets
post_b_assemb_CTI[1:10, 1:4]

# Summarize the posterior distributions for each ass.name
slopes_summary_CTI <- apply(post_b_assemb_CTI, 2, function(x) {
  c(median = median(x), sd = sd(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
})

slopes_summary_df_CTI <- as.data.frame(t(slopes_summary_CTI))
slopes_summary_df_CTI$ass.name <- sub(".*\\[([0-9_]+),.*", "\\1", rownames(slopes_summary_df_CTI))

head(slopes_summary_df_CTI)

## Add to assemb.data
plot(assemb.data$median_betaCTI, slopes_summary_df_CTI$median[match(assemb.data$ass.name, slopes_summary_df_CTI$ass.name)])

assemb.data$median_betaCTI <- slopes_summary_df_CTI$median[match(assemb.data$ass.name, slopes_summary_df_CTI$ass.name)]

# Proportion of assemblages with a positive slope
length(which(assemb.data$median_betaCTI > 0))
nrow(assemb.data)
length(which(assemb.data$median_betaCTI > 0)) / nrow(assemb.data)

# Visualize
plot_ass_cti <- ggplot(assemb.data, aes(x = median_betaCTI)) +
  geom_density(fill = "steelblue", alpha = 0.7) + # Smooth density curve with transparency
  geom_vline(aes(xintercept = median(median_betaCTI)), 
             color = "red", linetype = "solid", linewidth = 1) + # Median line
  geom_vline(xintercept = 0, 
             color = "grey30", linetype = "dashed", linewidth = 1) + # Median line
  labs(x = "Assemblage-level rate of change in CTI",
       y = "Density") +
  theme_classic(base_size = 14) + # Clean theme with larger base font size
  theme(axis.title = element_text(size = 14), 
        axis.text = element_text(size = 12)); plot_ass_cti

# Save data
write.csv(slopes_summary_df_CTI, "./results/slopes_summary_CTI.csv", row.names = FALSE)


# SAVE ASSEMB.DATA --------------------------------------------------------

write.csv(assemb.data, "./data/data_assemb_final.csv", row.names = FALSE)


# SAVE POSTERIOR TABLE --------------------------------------------------------

write.csv(post_table, "../Graphs/posteriors_table.csv", row.names = FALSE)

