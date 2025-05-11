## GENERALIZED LINEAR MIXED MODELS FOR STI ##

## Author: Philippe Fernandez-Fournier, 2025

#setwd("D:/sfuvault/SFU/MooersLab/Ch3_WinnersLosers/Analysis/")
rm(list=ls())

# Add package versions
library(ggplot2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(lme4)
library(splines)
library(marginaleffects)

## Load data 
comm.data <- read.csv("./data/data_pres_final.csv") 

assemb.data <- read.csv("./data/data_assemb_final.csv")
assemb.data$study <- as.factor(as.character(assemb.data$study))

species.data <- read.csv("./data/data_species_final.csv") 
species.data$study <- as.factor(species.data$study)

## Colors
col_scale_sti <- colorRampPalette(c("royalblue3", "firebrick3"))
col_delta_sti   <- col_scale_sti(11)[6] 

# Data for GLM
comm_data_glm <- comm.data %>%
  dplyr::select(species, ass.name, species_ass, study, presence, yr_obs, sti) %>%
  group_by(ass.name) %>% # Center by assemblage
  mutate(sti_cent = sti - mean(sti)) %>% 
  mutate(yr_obs_cent = yr_obs - mean(yr_obs)) %>% 
  ungroup() %>% 
  mutate(sti_scaled = (sti_cent - mean(sti_cent)) / sd(sti_cent)) %>% # scaled ed for modelling
  mutate(yr_obs_scaled = (yr_obs_cent - mean(yr_obs_cent)) / sd(yr_obs_cent)) # scaled year for modelling

# Add mean temperature
comm_data_glm$SST_mean <- assemb.data$SST_mean[match(comm_data_glm$ass.name, assemb.data$ass.name)]
comm_data_glm$SST_mean_scaled <- as.numeric(scale(comm_data_glm$SST_mean))

# Add change in SST
comm_data_glm$betaSST <- assemb.data$betaSST_lm[match(comm_data_glm$ass.name, assemb.data$ass.name)]
comm_data_glm$betaSST_scaled <- as.numeric(scale(comm_data_glm$betaSST))

# Add longitude and latitude
comm_data_glm$longitude <- assemb.data$longitude[match(comm_data_glm$ass.name, assemb.data$ass.name)]
comm_data_glm$abs_latitude <- assemb.data$abs_latitude[match(comm_data_glm$ass.name, assemb.data$ass.name)]
comm_data_glm$abs_lat_scaled <- as.numeric(scale(comm_data_glm$abs_latitude))

# Remove STI outliers
hist(comm_data_glm$sti_cent)
range(comm_data_glm$sti_cent)
hist(comm_data_glm$sti_scaled)

sti_q1 <- quantile(comm_data_glm$sti_scaled, 0.25, na.rm = TRUE)
sti_q3 <- quantile(comm_data_glm$sti_scaled, 0.75, na.rm = TRUE)
sti_iqr <- sti_q3 - sti_q1
sti_low <- sti_q1 - 3 * sti_iqr
sti_upp <- sti_q3 + 3 * sti_iqr

comm_data_glm <- subset(comm_data_glm, sti_scaled >= sti_low & sti_scaled <= sti_upp)
hist(comm_data_glm$sti_scaled)

# Remove betaSST outliers
hist(comm_data_glm$betaSST_scaled)

bsst_q1 <- quantile(comm_data_glm$betaSST_scaled, 0.25, na.rm = TRUE)
bsst_q3 <- quantile(comm_data_glm$betaSST_scaled, 0.75, na.rm = TRUE)
bsst_iqr <- bsst_q3 - bsst_q1
bsst_low <- bsst_q1 - 3 * bsst_iqr
bsst_upp <- bsst_q3 + 3 * bsst_iqr

comm_data_glm <- subset(comm_data_glm, betaSST_scaled >= bsst_low & betaSST_scaled <= bsst_upp)
hist(comm_data_glm$betaSST_scaled)

# Subset Atlantic regions
# West Atlantic assemblages
#comm_data_glm <- comm_data_glm %>% filter(longitude >= -100 & longitude <= -40)

# East Atlantic assemblages
#comm_data_glm <- comm_data_glm %>% filter(longitude >= -20 & longitude <= 30)

# Quantiles of STI to plot
sti_quants <- quantile(comm_data_glm$sti_scaled, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)

# for back-transformation
sti_mean <- mean(comm_data_glm$sti_cent)
sti_sd <- sd(comm_data_glm$sti_cent)
sti_value <- (sti_quants * sti_sd) + sti_mean

yr_obs_mean <- mean(comm_data_glm$yr_obs_cent)
yr_obs_sd <- sd(comm_data_glm$yr_obs_cent)

temp_mean <- mean(comm_data_glm$SST_mean)
temp_sd <- sd(comm_data_glm$SST_mean)

sst_mean <- mean(comm_data_glm$betaSST)
sst_sd <- sd(comm_data_glm$betaSST)


# STI DENSITY PLOT --------------------------------------------------------

data_quant <- data.frame(xintercept = sti_quants, quant = names(sti_quants))
data_quant$quant <- factor(data_quant$quant, levels = data_quant$quant)

plot_stiq <- ggplot(comm_data_glm, aes(x = sti_scaled)) +
  geom_density(fill = "grey50", color = NA, adjust = 3) +
  geom_vline(data = data_quant, 
             aes(xintercept = xintercept, color = quant), linewidth = 1) +
  scale_x_continuous(breaks = sti_quants, labels = round(sti_value, 0),
                     limits = c(-3.2,4)) + # Back-transformed values
  scale_color_manual(values = col_scale_sti(5)) +
  xlab("STI (centered)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.title.y = element_blank(),   
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.25), "in")); plot_stiq

ggsave(filename = '../Graphs/STI_quantiles.jpeg', plot_stiq,
       dpi = 300, width = 2, height = 1.5, units = "in")


# MODEL: MEAN TEMPERATURE ------------------------------------------------------

# Model
fit_pres_sti_temp_glm <- glmer(presence ~ yr_obs_scaled * sti_scaled * SST_mean_scaled + 
                                 (1 + yr_obs_scaled || ass.name) + (1 + yr_obs_scaled || ass.name:species_ass),
                               data = comm_data_glm, 
                               family = binomial(),
                               control = glmerControl(optimizer = "bobyqa",             # nloptwrap, Nelder_Mead (default), bobyqa
                                                      optCtrl = list(maxfun = 100000))) # Increase iterations

#save(fit_pres_sti_temp_glm, file = "./results/glm_betaPRES_STI_TEMP.RData")
load("./results/glm_betaPRES_STI_TEMP.RData")

summary(fit_pres_sti_temp_glm)

# OTHER MODELS ------------------------------------------------------------

## NONLINEAR 3-way interaction term
fit_pres_sti_temp_glm <- glmer(presence ~ yr_obs_scaled + sti_scaled + SST_mean_scaled + 
                                 yr_obs_scaled:sti_scaled + 
                                 yr_obs_scaled:SST_mean_scaled + 
                                 sti_scaled:SST_mean_scaled + 
                                 yr_obs_scaled:sti_scaled:ns(SST_mean_scaled, df = 3) + 
                                 (1 + yr_obs_scaled || ass.name) + (1 + yr_obs_scaled || ass.name:species_ass),
                               data = comm_data_glm, 
                               family = binomial(),
                               control = glmerControl(optimizer = "bobyqa",
                                                      optCtrl = list(maxfun = 100000)))

save(fit_pres_sti_temp_glm, file = "./results/glm_betaPRES_STI_polyTEMP.RData")
#load("./results/glm_betaPRES_STI_polyTEMP.RData")

summary(fit_pres_sti_temp_glm)


## NONLINEAR year * ED interacton term
fit_pres_sti_temp_glm <- glmer(presence ~ yr_obs_scaled + sti_scaled + SST_mean_scaled + 
                                 yr_obs_scaled:ns(sti_scaled, df = 3) + 
                                 yr_obs_scaled:SST_mean_scaled + 
                                 sti_scaled:SST_mean_scaled + 
                                 yr_obs_scaled:sti_scaled:SST_mean_scaled + 
                                 (1 + yr_obs_scaled || ass.name) + (1 + yr_obs_scaled || ass.name:species_ass),
                               data = comm_data_glm, 
                               family = binomial(),
                               control = glmerControl(optimizer = "bobyqa",
                                                      optCtrl = list(maxfun = 100000)))

save(fit_pres_sti_temp_glm, file = "./results/glm_betaPRES_polySTI_TEMP.RData")
#load("./results/glm_betaPRES_polySTI_TEMP.RData")

summary(fit_pres_sti_temp_glm)





# PREDICTIONS: MEAN TEMPERATURE ----------------------------------------------------------

# Temperature to plot
temp_values <- c(5, 10, 15, 20, 25)
temp_groups <- (temp_values - temp_mean) / temp_sd

# year obs to plot (99% of assemblages have 5-25 observed years)
yr_values <- c(-12.5, -10, -5, 0, 5, 10, 12.5)
yr_breaks <- (yr_values - yr_obs_mean) / yr_obs_sd

# Generate predictions
newdat_glm3 <- data.frame(yr_obs_scaled = rep(seq(min(yr_breaks), max(yr_breaks),
                                                  length.out = 200), (length(sti_quants) * length(temp_groups))),
                          sti_scaled = rep(rep(sti_quants, each = 200), length(temp_groups)),
                          SST_mean_scaled = rep(temp_groups, each = 200*length(sti_quants)))

newdat_glm3$SST_mean <- (newdat_glm3$SST_mean_scaled * temp_sd) + temp_mean

newdat_glm3$predicted <- predict(fit_pres_sti_temp_glm, newdata = newdat_glm3, re.form = NA) # exclude the random effects 
newdat_glm3$predicted_prob <- plogis(newdat_glm3$predicted)

# Temp labels
temp_labels <- setNames(paste0(temp_values, "°C"), temp_groups)

# Temp as factor
newdat_glm3$SST_mean_scaled <- factor(newdat_glm3$SST_mean_scaled, levels = unique(newdat_glm3$SST_mean_scaled))

# Plot one column
plot_sti_temp_quant2 <- ggplot(newdat_glm3, aes(x = yr_obs_scaled, y = predicted_prob, 
                                               group = factor(sti_scaled), color = sti_scaled)) +
  geom_line(linewidth = 1, alpha = 1) + 
  facet_wrap(~SST_mean_scaled, labeller = labeller(SST_mean_scaled = temp_labels), 
             ncol = 1, strip.position = "right") +
  scale_x_continuous(breaks = yr_breaks[2:6], labels = yr_values[2:6]) + # Back-transformed values
  scale_color_gradientn(colors = col_scale_sti(5)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  labs(x = "Year of observation", 
       y = "Predicted probability of presence") +
  theme_test() +
  theme(legend.position = "none",
        strip.text.y = element_text(size = 10, angle = 0),  
        strip.background = element_blank(),
        panel.spacing = unit(0.2, "lines")); plot_sti_temp_quant2

ggsave(filename = '../Graphs/pres_STI_TEMP.jpeg', plot_sti_temp_quant2,
       dpi = 300, width = 2.5, height = 5.5, units = "in")


# PLOT 3W INTERACTION MARGINAL EFFECTS -----------------------------------------------------

# Define grid of SST and STI values
pred_grid <- expand.grid(SST_mean_scaled = seq(min(comm_data_glm$SST_mean_scaled), max(comm_data_glm$SST_mean_scaled), length.out = 200),
  sti_scaled = sti_quants)

# Fix year at 0 (center)
pred_grid$yr_obs_scaled <- 0

# Estimate marginal effects of yr_obs_scaled at each SST × STI combo
mfx <- slopes(model = fit_pres_sti_temp_glm,
              variable = "yr_obs_scaled",
              newdata = pred_grid,
              re.form = NA)  # exclude random effects

temp_values2 <- c(0, 5, 10, 15, 20, 25)
temp_groups2 <- (temp_values2 - temp_mean) / temp_sd
temp_labels2 <- setNames(paste0(temp_values2, "°C"), temp_groups2)

# Plot: slope of year effect on presence by SST values, colored by STI quantiles
plot_marg <- ggplot(mfx, aes(x = SST_mean_scaled, y = estimate, color = factor(sti_scaled))) +
  geom_hline(yintercept = 0, linewidth = 1, color = "grey50", linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = factor(sti_scaled)), alpha = 0.2, color = NA) +
  labs(x = "Assemblage mean SST",
       y = "Marginal effect of year on presence (slope)",
       color = "STI (scaled)",
       fill = "STI (scaled)") +
  ylim(c(-0.053, 0.11)) +
  scale_x_continuous(breaks = temp_groups2, labels = temp_labels2) + # Back-transformed values
  scale_color_manual(values = col_scale_sti(5), 
                     name = "STI quantile", 
                     labels = names(sti_quants)) +
  scale_fill_manual(values = col_scale_sti(5), 
                    name = "STI quantile", 
                    labels = names(sti_quants)) +
  theme_bw(); plot_marg

ggsave(filename = '../Graphs/Supp_marginals_3way_TEMP_STI.jpeg', plot_marg,
       dpi = 300, width = 5.5, height = 3.5, units = "in")


# MODEL: betaSST ------------------------------------------------------

hist(comm_data_glm$betaSST,
     #main = "East Atlantic betaSST",
     xlab = "betaSST")

# Model
fit_pres_sti_betasst_glm <- glmer(presence ~ yr_obs_scaled * sti_scaled * betaSST_scaled + 
                                (1 + yr_obs_scaled || ass.name) + (1 + yr_obs_scaled || ass.name:species_ass),
                              data = comm_data_glm, 
                              family = binomial(),
                              control = glmerControl(optimizer = "bobyqa",            # nloptwrap, Nelder_Mead (default), bobyqa
                                                     optCtrl = list(maxfun = 100000))) # Increase iterations

#save(fit_pres_sti_betasst_glm, file = "./results/glm_betaPRES_STI_betaSST.RData")
load("./results/glm_betaPRES_STI_betaSST.RData")

summary(fit_pres_sti_betasst_glm)


# PREDICTIONS: betaSST ----------------------------------------------------------

range(comm_data_glm$betaSST)

# for back-transformation
sst_mean <- mean(comm_data_glm$betaSST)
sst_sd <- sd(comm_data_glm$betaSST)

# Latitudes to plot
sst_values <- c(-0.5, 0, 0.5, 1, 1.5)
sst_groups <- (sst_values - sst_mean) / sst_sd

# year obs to plot (99% of assemblages have 5-25 observed years)
yr_values <- c(-12.5, -10, -5, 0, 5, 10, 12.5)
yr_breaks <- (yr_values - yr_obs_mean) / yr_obs_sd

# Generate predictions
newdat_glm4 <- data.frame(yr_obs_scaled = rep(seq(min(yr_breaks), max(yr_breaks),
                                                  length.out = 200), (length(sti_quants) * length(sst_groups))),
                          sti_scaled = rep(rep(sti_quants, each = 200), length(sst_groups)),
                          betaSST_scaled = rep(sst_groups, each = 200*length(sti_quants)))

newdat_glm4$betaSST <- (newdat_glm4$betaSST_scaled * sst_sd) + sst_mean

newdat_glm4$predicted <- predict(fit_pres_sti_betasst_glm, newdata = newdat_glm4, re.form = NA) # exclude the random effects 
newdat_glm4$predicted_prob <- plogis(newdat_glm4$predicted)

# Latitude labels
sst_labels <- setNames(paste0(c("", rep("+", 4)), sst_values, "°C"), sst_groups)

# Latitude as factor
newdat_glm4$betaSST_scaled <- factor(newdat_glm4$betaSST_scaled, levels = rev(unique(newdat_glm4$betaSST_scaled)))

# Plot
plot_sti_sst_quant2 <- ggplot(newdat_glm4, aes(x = yr_obs_scaled, y = predicted_prob, 
                                               group = factor(sti_scaled), color = sti_scaled)) +
  geom_line(linewidth = 1, alpha = 1) + 
  facet_wrap(~betaSST_scaled, labeller = labeller(betaSST_scaled = sst_labels), 
             ncol = 1, strip.position = "right") +
  scale_x_continuous(breaks = yr_breaks[2:6], labels = yr_values[2:6]) + # Back-transformed values
  scale_color_gradientn(colors = col_scale_sti(5)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  labs(x = "Year of observation (centered)", 
       y = "Predicted probability of presence") +
  theme_test() +
  theme(legend.position = "none",
        strip.text.y = element_text(size = 10, angle = 0),  # Ensure correct element is targeted
        strip.background = element_blank(),
        panel.spacing = unit(0.2, "lines")); plot_sti_sst_quant2

ggsave(filename = '../Graphs/pres_STI_betaSST.jpeg', plot_sti_sst_quant2,
       dpi = 300, width = 2.5, height = 5.5, units = "in")



# SPECIES betaPRES ---------------------------------------
# Based on 3-way interaction model with latitude
# species_ass coefficients (Sum of the random and fixed effects for each level)
fit_bp_sti_coef <- coef(fit_pres_sti_temp_glm)

coef_sp_ass <- fit_bp_sti_coef$`ass.name:species_ass` %>%
  rownames_to_column("assemb_species_ass") %>%
  separate(assemb_species_ass, into = c("ass.name", "species_ass"), sep = ":") %>%
  rename(intercept = '(Intercept)',
         slope     = yr_obs_scaled) %>%
  dplyr::select(species_ass, ass.name, intercept, slope)

# add betaPRES coefficient to species_ass dataframe
sp_ass_dat <- comm_data_glm %>%
  dplyr::select(species_ass, sti_scaled) %>%
  distinct() %>%
  left_join(coef_sp_ass, by = "species_ass")

# Add species slopes to species.data 
species.data$betaPRES <- coef_sp_ass$slope[match(species.data$species_ass, coef_sp_ass$species_ass)]

hist(species.data$betaPRES)


# SPECIES deltaPRES ---------------------------------------

newdat_glm_dp <- comm_data_glm
newdat_glm_dp$predicted <- predict(fit_pres_sti_temp_glm) 
newdat_glm_dp$predicted_prob <- plogis(newdat_glm_dp$predicted)

# First and last obs
newdat_glm_dp <- newdat_glm_dp %>%
  group_by(species_ass) %>%
  arrange(yr_obs_scaled) %>%
  mutate(pres_first = predicted_prob[which.min(yr_obs_scaled)],  # First observation
         pres_last = predicted_prob[which.max(yr_obs_scaled)]) %>%     # Last observation
  mutate(delta_pres = pres_last - pres_first) %>% 
  ungroup() %>% 
  dplyr::select(species_ass, ass.name, sti_scaled, delta_pres) %>%
  distinct() 

save(newdat_glm_dp, file = "./results/glm_turnover_predict_sti.RData")

# Add to species.data
species.data$deltaPRES <- newdat_glm_dp$delta_pres[match(species.data$species_ass, newdat_glm_dp$species_ass)]

hist(species.data$deltaPRES)

plot(species.data$deltaPRES, species.data$betaPRES)


# SAVE SPECIES.DATA --------------------------------------------

write.csv(species.data, "./data/data_species_final.csv", row.names = FALSE) 




# MODEL: LATITUDE ------------------------------------------------------

# species_ass random intercept and slope nested within assemblage uncorrelated with "||" (corr is -0.14 when they are correlated)
fit_pres_sti_lat_glm <- glmer(presence ~ yr_obs_scaled*sti_scaled*abs_lat_scaled + 
                                (1 + yr_obs_scaled || ass.name) + (1 + yr_obs_scaled || ass.name:species_ass),
                              data = comm_data_glm, 
                              family = binomial(),
                              control = glmerControl(optimizer = "bobyqa",             # nloptwrap, Nelder_Mead (default), bobyqa
                                                     optCtrl = list(maxfun = 100000))) # Increase iterations

#save(fit_pres_sti_lat_glm, file = "./results/glm_betaPRES_STI_LAT.RData")
load("./results/glm_betaPRES_STI_LAT.RData")

summary(fit_pres_sti_lat_glm)


# PREDICTIONS: LATITUDE ----------------------------------------------------------

# for back-transformation
lat_mean <- mean(comm_data_glm$abs_latitude)
lat_sd <- sd(comm_data_glm$abs_latitude)

# Latitudes to plot
lat_values <- c(20, 30, 40, 50, 60)
lat_groups <- (lat_values - lat_mean) / lat_sd

# year obs to plot (99% of assemblages have 5-25 observed years)
yr_values <- c(-12.5, -10, -5, 0, 5, 10, 12.5)
yr_breaks <- (yr_values - yr_obs_mean) / yr_obs_sd

# Generate predictions
newdat_glm3 <- data.frame(yr_obs_scaled = rep(seq(min(yr_breaks), max(yr_breaks),
                                                  length.out = 200), (length(sti_quants) * length(lat_groups))),
                          sti_scaled = rep(rep(sti_quants, each = 200), length(lat_groups)),
                          abs_lat_scaled = rep(lat_groups, each = 200*length(sti_quants)))

newdat_glm3$abs_latitude <- (newdat_glm3$abs_lat_scaled * lat_sd) + lat_mean

newdat_glm3$predicted <- predict(fit_pres_sti_lat_glm, newdata = newdat_glm3, re.form = NA) # exclude the random effects 
newdat_glm3$predicted_prob <- plogis(newdat_glm3$predicted)

# Temp labels
lat_labels <- setNames(paste0(lat_values, "°"), lat_groups)

# Temp as factor
newdat_glm3$abs_lat_scaled <- factor(newdat_glm3$abs_lat_scaled, levels = rev(unique(newdat_glm3$abs_lat_scaled)))


# Plot one column
plot_sti_lat_quant2 <- ggplot(newdat_glm3, aes(x = yr_obs_scaled, y = predicted_prob, 
                                               group = factor(sti_scaled), color = sti_scaled)) +
  geom_line(linewidth = 1, alpha = 1) + 
  facet_wrap(~abs_lat_scaled, labeller = labeller(abs_lat_scaled = lat_labels), 
             ncol = 1, strip.position = "right") +
  scale_x_continuous(breaks = yr_breaks[2:6], labels = yr_values[2:6]) + # Back-transformed values
  scale_color_gradientn(colors = col_scale_sti(5)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  labs(x = "Year of observation (centered)", 
       y = "Predicted probability of presence") +
  theme_test() +
  theme(legend.position = "none",
        strip.text.y = element_text(size = 10, angle = 0),  # Ensure correct element is targeted
        strip.background = element_blank(),
        panel.spacing = unit(0.2, "lines")); plot_sti_lat_quant2

ggsave(filename = '../Graphs/pres_STI_LAT.jpeg', plot_sti_lat_quant2,
       dpi = 300, width = 2.5, height = 5.5, units = "in")






# PLOT: SPECIES TRENDS ~ LAT -----------------------------------------------
colnames(comm_data_glm)
newdat_glm_sp3 <- comm_data_glm
newdat_glm_sp3$predicted <- predict(fit_pres_sti_temp_glm) 
newdat_glm_sp3$predicted_prob <- plogis(newdat_glm_sp3$predicted)

# Species in 0-5% quantiles of STI
sp_05 <- unique(newdat_glm_sp3$species_ass[newdat_glm_sp3$sti_scaled <= sti_quants[1]])

# Species in 95-100% quantiles of STI
sp_95 <- unique(newdat_glm_sp3$species_ass[newdat_glm_sp3$sti_scaled >= sti_quants[5]])

# Latitude bins
latbins <- 10  # How big should the intervals be (in degrees)?
trop_lat <- 15 # join low equator latitudes below this value because of low sample size
custom_breaks <- c(0, trop_lat, seq(trop_lat + latbins, 65, by = latbins))


### Plot 5% ###
plot_data05 <- newdat_glm_sp3[newdat_glm_sp3$species_ass %in% sp_05,]

# Normalize yr_obs_scaled
plot_data05 <- plot_data05 %>% 
  group_by(species_ass) %>% 
  mutate(yr_obs_scaled_norm = (yr_obs_scaled - min(yr_obs_scaled)) / (max(yr_obs_scaled) - min(yr_obs_scaled))) %>% 
  ungroup() %>% 
  group_by(ass.name) %>% 
  mutate(yr_obs_scaled_cent = yr_obs_scaled - mean(yr_obs_scaled))

# Latitude bins
plot_data05$latitude_bin <- cut(plot_data05$abs_latitude, breaks = custom_breaks, include.lowest = TRUE, right = FALSE)
plot_data05$latitude_bin <- factor(plot_data05$latitude_bin, levels = rev(levels(plot_data05$latitude_bin)))

# Plot 
plot_spp05 <- ggplot(plot_data05, aes(x = yr_obs_scaled_norm, y = presence)) +
  geom_point(shape = 108, size = 2) +
  geom_line(aes(x = yr_obs_scaled_norm, y = predicted_prob, group = species_ass), 
            linewidth = 0.5, color = "grey50", alpha = 0.1) +  
  stat_smooth(method = "glm", method.args = list(family = binomial()), 
              se = FALSE, color = "royalblue3", linewidth = 1) +
  facet_wrap(~ latitude_bin, ncol = 1, scales = "free_x", strip.position = "right") +
  labs(x = "Year of observation", 
       y = "Presence") +
  scale_y_continuous(breaks = c(0, 1), limits = c(0,1)) +
  theme_test() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines")); plot_spp05

ggsave(filename = '../Graphs/species_logits_temp.jpeg', plot_spp05,
       dpi = 300, width = 2, height = length(unique(plot_data05$latitude_bin)), units = "in")

### Plot 95% ###
plot_data95 <- newdat_glm_sp3[newdat_glm_sp3$species_ass %in% sp_95,]

# Normalize yr_obs_scaled
plot_data95 <- plot_data95 %>% 
  group_by(species_ass) %>% 
  mutate(yr_obs_scaled_norm = (yr_obs_scaled - min(yr_obs_scaled)) / (max(yr_obs_scaled) - min(yr_obs_scaled))) %>% 
  ungroup() %>% 
  group_by(ass.name) %>% 
  mutate(yr_obs_scaled_cent = yr_obs_scaled - mean(yr_obs_scaled))

# Latitude bins
plot_data95$latitude_bin <- cut(plot_data95$abs_latitude, breaks = custom_breaks, include.lowest = TRUE, right = FALSE)
plot_data95$latitude_bin <- factor(plot_data95$latitude_bin, levels = rev(levels(plot_data95$latitude_bin)))

# Plot 
plot_spp95 <- ggplot(plot_data95, aes(x = yr_obs_scaled_norm, y = presence)) +
  geom_point(shape = 108, size = 2) +
  geom_line(aes(x = yr_obs_scaled_norm, y = predicted_prob, group = species_ass), 
            linewidth = 0.5, color = "grey50", alpha = 0.1) +  
  stat_smooth(method = "glm", method.args = list(family = binomial()),
              se = FALSE, color = "firebrick3", linewidth = 1) +
  facet_wrap(~ latitude_bin, ncol = 1, scales = "free_x", strip.position = "right") +
  labs(x = "Year of observation", 
       y = "Presence") +
  scale_y_continuous(breaks = c(0, 1), limits = c(0,1)) +
  theme_test() +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank(),
        strip.text = element_text(size = 8),
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines")); plot_spp95

ggsave(filename = '../Graphs/species_logits_n100_sti_lat_95.jpeg', plot_spp95,
       dpi = 300, width = 2, height = length(unique(plot_data95$latitude_bin)), units = "in")


