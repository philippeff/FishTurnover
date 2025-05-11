## GENERALIZED LINEAR MIXED MODELS FOR LOCAL ED ##

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
col_scale_bpres <- c("#762A83", "#b179cb", "grey70", "#65c26d", "#1c7837")
col_delta_bpres <- col_scale_bpres[length(col_scale_bpres)]

# Data for GLM
# Center year and ED then scale?
comm_data_glm <- comm.data %>%
  dplyr::select(species, ass.name, species_ass, study, presence, yr_obs, ed.std) %>%
  mutate(ed_trans = sqrt(ed.std)) %>% # sqrt-transformed ed
  group_by(ass.name) %>% # Center by assemblage
  mutate(ed_trans_cent = ed_trans - mean(ed_trans)) %>% 
  mutate(yr_obs_cent = yr_obs - mean(yr_obs)) %>% 
  ungroup() %>% 
  mutate(ed_trans_scaled = (ed_trans_cent - mean(ed_trans_cent)) / sd(ed_trans_cent)) %>% # scaled ed for modelling
  mutate(yr_obs_scaled = (yr_obs_cent - mean(yr_obs_cent)) / sd(yr_obs_cent)) # scaled year for modelling

# Add longitude and latitude
comm_data_glm$longitude <- assemb.data$longitude[match(comm_data_glm$ass.name, assemb.data$ass.name)]
comm_data_glm$abs_latitude <- assemb.data$abs_latitude[match(comm_data_glm$ass.name, assemb.data$ass.name)]
comm_data_glm$abs_lat_scaled <- as.numeric(scale(comm_data_glm$abs_latitude))

# Add mean temperature
comm_data_glm$SST_mean <- assemb.data$SST_mean[match(comm_data_glm$ass.name, assemb.data$ass.name)]
comm_data_glm$SST_mean_scaled <- as.numeric(scale(comm_data_glm$SST_mean))

# Add change in SST
comm_data_glm$betaSST <- assemb.data$betaSST_lm[match(comm_data_glm$ass.name, assemb.data$ass.name)]
comm_data_glm$betaSST_scaled <- as.numeric(scale(comm_data_glm$betaSST))

# SR per assemblage (log-transform and scale)
comm_data_glm <- comm_data_glm %>%
  group_by(ass.name) %>%
  mutate(assemb_SR = n_distinct(species)) %>% 
  ungroup()

comm_data_glm$SR_scaled <- as.numeric(scale(log(comm_data_glm$assemb_SR)))
hist(comm_data_glm$SR_scaled)

# Remove betaSST outliers (others are OK)
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


# Quantiles of ED to plot
ed_quants <- quantile(comm_data_glm$ed_trans_scaled, probs = c(0.05, 0.25, 0.50, 0.75, 0.95), na.rm = TRUE)

# Means and SDs for back-transformation
ed_mean <- mean(comm_data_glm$ed_trans_cent)
ed_sd <- sd(comm_data_glm$ed_trans_cent)
ed_value <- (ed_quants * ed_sd) + ed_mean

yr_obs_mean <- mean(comm_data_glm$yr_obs_cent)
yr_obs_sd <- sd(comm_data_glm$yr_obs_cent)

temp_mean <- mean(comm_data_glm$SST_mean)
temp_sd <- sd(comm_data_glm$SST_mean)

sst_mean <- mean(comm_data_glm$betaSST)
sst_sd <- sd(comm_data_glm$betaSST)


# ED DENSITY PLOT --------------------------------------------------------

data_quant <- data.frame(xintercept = ed_quants, quant = names(ed_quants))
data_quant$quant <- factor(c("5th", "25th", "50th", "75th", "95th"),
                           levels = c("5th", "25th", "50th", "75th", "95th"))

plot_edq <- ggplot(comm_data_glm, aes(x = ed_trans_scaled)) +
  geom_density(fill = "grey50", color = NA, adjust = 3) +
  geom_vline(data = data_quant, 
             aes(xintercept = xintercept, color = quant), linewidth = 1) +
  scale_x_continuous(breaks = ed_quants[c(1,3,5)], labels = round(ed_value[c(1,3,5)], 2)) + # Back-transformed values
  scale_color_manual(values = col_scale_bpres) +
  xlab("Local ED (sqrt, centered)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.line.x.bottom = element_line(colour = "black"),
        axis.title.y = element_blank(),   
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8),
        plot.margin = unit(c(0.1, 0.25, 0.1, 0.1), "in")); plot_edq

ggsave(filename = '../Graphs/ED_quantiles.jpeg', plot_edq, dpi = 300, width = 2, height = 1.5, units = "in")


# MODEL: MEAN TEMPERATURE ------------------------------------------------------

# Model
fit_bp_ed_temp_glm <- glmer(presence ~ yr_obs_scaled * ed_trans_scaled * SST_mean_scaled + 
                              SR_scaled + 
                              (1 + yr_obs_scaled || ass.name) + (1 + yr_obs_scaled||ass.name:species_ass),
                            data = comm_data_glm, 
                            family = binomial(),
                            control = glmerControl(optimizer = "bobyqa",            # nloptwrap, Nelder_Mead (default), bobyqa
                                                   optCtrl = list(maxfun = 100000))) # Increase iterations

#save(fit_bp_ed_temp_glm, file = "./results/glm_betaPRES_ED_TEMP.RData")
load("./results/glm_betaPRES_ED_TEMP.RData")

summary(fit_bp_ed_temp_glm)


# OTHER MODELS ------------------------------------------------------------

## NONLINEAR 3-way interaction term
fit_bp_ed_temp_glm <- glmer(presence ~ yr_obs_scaled + ed_trans_scaled + SST_mean_scaled + 
                              yr_obs_scaled:ed_trans_scaled + 
                              yr_obs_scaled:SST_mean_scaled + 
                              ed_trans_scaled:SST_mean_scaled + 
                              yr_obs_scaled:ed_trans_scaled:ns(SST_mean_scaled, df = 3) + 
                              SR_scaled + 
                              (1 + yr_obs_scaled || ass.name) + (1 + yr_obs_scaled || ass.name:species_ass),
                            data = comm_data_glm, 
                            family = binomial(),
                            control = glmerControl(optimizer = "bobyqa",
                                                   optCtrl = list(maxfun = 100000)))

save(fit_bp_ed_temp_glm, file = "./results/glm_betaPRES_ED_polyTEMP.RData")
#load("./results/glm_betaPRES_ED_polyTEMP.RData")

summary(fit_bp_ed_temp_glm)


## NONLINEAR year * ED interacton term
fit_bp_ed_temp_glm <- glmer(presence ~ yr_obs_scaled + ed_trans_scaled + SST_mean_scaled + 
                              yr_obs_scaled:ns(ed_trans_scaled, df = 3) + 
                              yr_obs_scaled:SST_mean_scaled + 
                              ed_trans_scaled:SST_mean_scaled + 
                              yr_obs_scaled:ed_trans_scaled:SST_mean_scaled + 
                              SR_scaled + 
                              (1 + yr_obs_scaled || ass.name) + (1 + yr_obs_scaled || ass.name:species_ass),
                            data = comm_data_glm, 
                            family = binomial(),
                            control = glmerControl(optimizer = "bobyqa",
                                                   optCtrl = list(maxfun = 100000)))

save(fit_bp_ed_temp_glm, file = "./results/glm_betaPRES_polyED_TEMP.RData")
#load("./results/glm_betaPRES_polyED_TEMP.RData")

summary(fit_bp_ed_temp_glm)


# PREDICTIONS: MEAN TEMPERATURE ----------------------------------------------------------

# Temperature to plot
temp_values <- c(5, 10, 15, 20, 25)
temp_groups <- (temp_values - temp_mean) / temp_sd

# year obs to plot
yr_values <- c(-12.5, -10, -5, 0, 5, 10, 12.5)
yr_breaks <- (yr_values - yr_obs_mean) / yr_obs_sd

# Generate predictions
newdat_glm3 <- data.frame(yr_obs_scaled   = rep(seq(min(yr_breaks), max(yr_breaks), length.out = 200), 
                                                (length(ed_quants) * length(temp_groups))),
                          ed_trans_scaled = rep(rep(ed_quants, each = 200), 
                                                length(temp_groups)),
                          SST_mean_scaled = rep(temp_groups, each = 200*length(ed_quants)))

newdat_glm3$SST_mean <- (newdat_glm3$SST_mean_scaled * temp_sd) + temp_mean
newdat_glm3$SR_scaled <- 0  # Holding SR constant at the mean

newdat_glm3$predicted <- predict(fit_bp_ed_temp_glm, newdata = newdat_glm3, re.form = NA) # exclude the random effects 
newdat_glm3$predicted_prob <- plogis(newdat_glm3$predicted)

# Temperature labels
temp_labels <- setNames(paste0(temp_values, "°C"), temp_groups)

# Temperature as factor
newdat_glm3$SST_mean_scaled <- factor(newdat_glm3$SST_mean_scaled, levels = unique(newdat_glm3$SST_mean_scaled))


# Plot 
plot_ed_temp_quant2 <- ggplot(newdat_glm3, aes(x = yr_obs_scaled, y = predicted_prob, 
                                               group = factor(ed_trans_scaled), color = ed_trans_scaled)) +
  geom_line(linewidth = 1, alpha = 1) + 
  facet_wrap(~SST_mean_scaled, labeller = labeller(SST_mean_scaled = temp_labels), 
             ncol = 1, strip.position = "right") +
  scale_x_continuous(breaks = yr_breaks[2:6], labels = yr_values[2:6]) + 
  scale_color_gradientn(colors = col_scale_bpres) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  labs(x = "Year of observation (centered)", 
       y = "Predicted probability of presence") +
  theme_test() +
  theme(legend.position = "none",
        strip.text.y = element_text(size = 10, angle = 0),
        strip.background = element_blank(),
        panel.spacing = unit(0.2, "lines")); plot_ed_temp_quant2

ggsave(filename = '../Graphs/pres_ED_TEMP.jpeg', plot_ed_temp_quant2,
       dpi = 300, width = 2.5, height = 5.5, units = "in")

# PLOT 3W INTERACTION MARGINAL EFFECTS -----------------------------------------------------

# Define grid of SST and ED values
pred_grid <- expand.grid(SST_mean_scaled = seq(min(comm_data_glm$SST_mean_scaled), max(comm_data_glm$SST_mean_scaled), length.out = 200),
                         ed_trans_scaled = ed_quants)

# Fix year at 0 (center)
pred_grid$yr_obs_scaled <- 0

# Fix SR at 0 (mean)
pred_grid$SR_scaled <- 0

# Estimate marginal effects of yr_obs_scaled at each SST × ED combo
mfx <- slopes(model = fit_bp_ed_temp_glm,
              variable = "yr_obs_scaled",
              newdata = pred_grid,
              re.form = NA)  # exclude random effects

temp_values2 <- c(0, 5, 10, 15, 20, 25)
temp_groups2 <- (temp_values2 - temp_mean) / temp_sd
temp_labels2 <- setNames(paste0(temp_values2, "°C"), temp_groups2)

# Plot: slope of year effect on presence by SST values, colored by ED quantiles
plot_marg <- ggplot(mfx, aes(x = SST_mean_scaled, y = estimate, color = factor(ed_trans_scaled))) +
  geom_hline(yintercept = 0, linewidth = 1, color = "grey50", linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = factor(ed_trans_scaled)), alpha = 0.2, color = NA) +
  labs(x = "Assemblage mean SST",
       y = "Marginal effect of year on presence (slope)") +
  scale_x_continuous(breaks = temp_groups2, labels = temp_labels2) + # Back-transformed values
  ylim(c(-0.052, 0.11)) +
  scale_color_manual(values = col_scale_bpres, 
                     name = "ED quantile ", 
                     labels = names(ed_quants)) +
  scale_fill_manual(values = col_scale_bpres, 
                    name = "ED quantile ", 
                    labels = names(ed_quants)) +
  theme_bw(); plot_marg

ggsave(filename = '../Graphs/Supp_marginals_3way_TEMP_ED.jpeg', plot_marg,
       dpi = 300, width = 5.5, height = 3.5, units = "in")



# MODEL: betaSST ------------------------------------------------------

# Model
fit_bp_ed_betasst_glm <- glmer(presence ~ yr_obs_scaled * ed_trans_scaled * betaSST_scaled + 
                                 SR_scaled + 
                                 (1 + yr_obs_scaled || ass.name) + (1 + yr_obs_scaled||ass.name:species_ass),
                               data = comm_data_glm, 
                               family = binomial(),
                               control = glmerControl(optimizer = "bobyqa",            # nloptwrap, Nelder_Mead (default), bobyqa
                                                      optCtrl = list(maxfun = 100000))) # Increase iterations

#save(fit_bp_ed_betasst_glm, file = "./results/glm_betaPRES_ED_betaSST.RData")
load("./results/glm_betaPRES_ED_betaSST.RData")

summary(fit_bp_ed_betasst_glm)


# PREDICTIONS: betaSST ----------------------------------------------------------

range(comm_data_glm$betaSST)

# Latitudes to plot
sst_values <- c(-0.5, 0, 0.5, 1, 1.5)
sst_groups <- (sst_values - sst_mean) / sst_sd

# year obs to plot (99% of assemblages have 5-25 observed years)
yr_values <- c(-12.5, -10, -5, 0, 5, 10, 12.5)
yr_breaks <- (yr_values - yr_obs_mean) / yr_obs_sd

# Generate predictions
newdat_glm4 <- data.frame(yr_obs_scaled   = rep(seq(min(yr_breaks), max(yr_breaks), length.out = 200), 
                                                (length(ed_quants) * length(sst_groups))),
                          ed_trans_scaled = rep(rep(ed_quants, each = 200), length(sst_groups)),
                          betaSST_scaled  = rep(sst_groups, each = 200*length(ed_quants)))

newdat_glm4$betaSST <- (newdat_glm4$betaSST_scaled * sst_sd) + sst_mean
newdat_glm4$SR_scaled <- 0  # Holding SR constant at the mean

newdat_glm4$predicted <- predict(fit_bp_ed_betasst_glm, newdata = newdat_glm4, re.form = NA) # exclude the random effects 
newdat_glm4$predicted_prob <- plogis(newdat_glm4$predicted)

# Latitude labels
sst_labels <- setNames(paste0(c("", rep("+", 4)), sst_values, "°C"), sst_groups)

# Latitude as factor
newdat_glm4$betaSST_scaled <- factor(newdat_glm4$betaSST_scaled, levels = rev(unique(newdat_glm4$betaSST_scaled)))

# Plot 
plot_ed_sst_quant2 <- ggplot(newdat_glm4, aes(x = yr_obs_scaled, y = predicted_prob, 
                                               group = factor(ed_trans_scaled), color = ed_trans_scaled)) +
  geom_line(linewidth = 1, alpha = 1) + 
  facet_wrap(~betaSST_scaled, labeller = labeller(betaSST_scaled = sst_labels), 
             ncol = 1, strip.position = "right") +
  scale_x_continuous(breaks = yr_breaks[2:6], labels = yr_values[2:6]) + 
  scale_color_gradientn(colors = col_scale_bpres) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  labs(x = "Year of observation", 
       y = "Predicted probability of presence") +
  theme_test() +
  theme(legend.position = "none",
        strip.text.y = element_text(size = 10, angle = 0),  # Ensure correct element is targeted
        strip.background = element_blank(),
        panel.spacing = unit(0.2, "lines")); plot_ed_sst_quant2

ggsave(filename = '../Graphs/pres_ED_betaSST.jpeg', plot_ed_sst_quant2,
       dpi = 300, width = 2.5, height = 5.5, units = "in")


# MODEL: LATITUDE ------------------------------------------------------

# Model
fit_bp_ed_lat_glm <- glmer(presence ~ yr_obs_scaled*ed_trans_scaled*abs_lat_scaled + 
                             SR_scaled + 
                             (1 + yr_obs_scaled || ass.name) + (1 + yr_obs_scaled||ass.name:species_ass),
                           data = comm_data_glm, 
                           family = binomial(),
                           control = glmerControl(optimizer = "bobyqa",            # nloptwrap, Nelder_Mead (default), bobyqa
                                                  optCtrl = list(maxfun = 100000))) # Increase iterations

#save(fit_bp_ed_lat_glm, file = "./results/glm_betaPRES_ED_LAT.RData")
load("./results/glm_betaPRES_ED_LAT.RData")

summary(fit_bp_ed_lat_glm)


# PREDICTIONS: LATITUDE  ----------------------------------------------------------

# for back-transformation
lat_mean <- mean(comm_data_glm$abs_latitude)
lat_sd <- sd(comm_data_glm$abs_latitude)

# Latitudes to plot
lat_values <- c(20, 30, 40, 50, 60)
lat_groups <- (lat_values - lat_mean) / lat_sd

# year obs to plot
yr_values <- c(-12.5, -10, -5, 0, 5, 10, 12.5)
yr_breaks <- (yr_values - yr_obs_mean) / yr_obs_sd

# Generate predictions
newdat_glm3 <- data.frame(yr_obs_scaled = rep(seq(min(yr_breaks), max(yr_breaks), length.out = 200), 
                                              (length(ed_quants) * length(lat_groups))),
                          ed_trans_scaled = rep(rep(ed_quants, each = 200), length(lat_groups)),
                          abs_lat_scaled = rep(lat_groups, each = 200*length(ed_quants)))

newdat_glm3$abs_latitude <- (newdat_glm3$abs_lat_scaled * lat_sd) + lat_mean
newdat_glm3$SR_scaled <- 0  # Holding SR constant at the mean

newdat_glm3$predicted <- predict(fit_bp_ed_lat_glm, newdata = newdat_glm3, re.form = NA) # exclude the random effects 
newdat_glm3$predicted_prob <- plogis(newdat_glm3$predicted)

# Latitude labels
lat_labels <- setNames(paste0(lat_values, "°"), lat_groups)

# Latitude as factor
newdat_glm3$abs_lat_scaled <- factor(newdat_glm3$abs_lat_scaled, levels = rev(unique(newdat_glm3$abs_lat_scaled)))


# Plot
plot_ed_lat_quant2 <- ggplot(newdat_glm3, aes(x = yr_obs_scaled, y = predicted_prob, 
                                              group = factor(ed_trans_scaled), color = ed_trans_scaled)) +
  geom_line(linewidth = 1, alpha = 1) + 
  facet_wrap(~abs_lat_scaled, labeller = labeller(abs_lat_scaled = lat_labels), 
             ncol = 1, strip.position = "right") +
  scale_x_continuous(breaks = yr_breaks[2:6], labels = yr_values[2:6]) + # Back-transformed values
  scale_color_gradientn(colors = col_scale_bpres) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1)) +
  labs(x = "Year of observation (centered)", 
       y = "Predicted probability of presence") +
  theme_test() +
  theme(legend.position = "none",
        strip.text.y = element_text(size = 10, angle = 0),  # Ensure correct element is targeted
        strip.background = element_blank(),
        panel.spacing = unit(0.2, "lines")); plot_ed_lat_quant2

ggsave(filename = '../Graphs/pres_ED_LAT.jpeg', plot_ed_lat_quant2,
       dpi = 300, width = 2.5, height = 5.5, units = "in")





