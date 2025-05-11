#*################################################*#
## CALCULATE EVOLUTIONARI DISTINCTIVENESS METRICS ##
#*################################################*#

## Author: Philippe Fernandez-Fournier, 2025

##  This script is designed to run on a Digital Alliance Canada cluster named Cedar
##  See shell script for job submission to the cluster


#setwd("D:/sfuvault/SFU/MooersLab/Ch3_WinnersLosers/Analysis")
rm(list=ls())

# load libraries
library(tidyr)
library(dplyr)
library(ape)
library(picante)
library(doParallel)

# Set env
ncores <- Sys.getenv("SLURM_CPUS_PER_TASK") 
# ncores <- detectCores()-1
registerDoParallel(cores=ncores)

# load data
comm.data <- read.csv("./data/data_pres.csv")
species.data <- read.csv("./data/data_species.csv")

# load phylogenies
phylo_fish  <- read.tree("./data/phylo_fish_PFF.tre")

## Go over each assemblage and calculate local ED and the contribution to local MPD

ass.unique <- levels(factor(comm.data$assemblage))

# Subset
#ass.unique <- levels(factor(comm.data$assemblage))[1:10]

species_evoldata <- foreach::foreach(r = 1:length(ass.unique), .combine = rbind, .packages = c('dplyr', 'ape', 'picante')) %dopar% {
  # r = 2
  # track 
  cat("\r", paste0(r, "/", length(ass.unique)))
  
  ass <- ass.unique[r]
  ass.name_r <- unique(comm.data$ass.name[which(comm.data$assemblage == ass)])
  species.data_r <- species.data[which(species.data$ass.name == ass.name_r),]
  
  # species names in assemblage r
  spp.ass <- unique(comm.data$sp.name[which(comm.data$assemblage == ass)])
  # Filter out species not in phylogeny (should not remove any, but just in case)
  spp.ass <- spp.ass[spp.ass %in% phylo_fish[[1]]$tip.label]
  
  #comm.tree <- keep.tip(phylo_fish[[1]], spp.ass)
  #plot(comm.tree, cex = 0.6)
  
  # Calculate local ED across all trees of phylo_fish
  ed.df <- data.frame()
  
  for (p in 1:length(phylo_fish)) {
    # p = 1
    # Extract and prune the phylogenetic tree
    phy.x <- keep.tip(phylo_fish[[p]], spp.ass)
    
    # Calculate evolutionary distinctiveness (equal splits)
    ed.x <- evol.distinct(phy.x, type = "equal.splits")
    ed.x$ed.std <- ed.x$w / sum(phy.x$edge.length)
    
    # Calculate MPD of the whole assemblage
    coph_all <- ape::cophenetic.phylo(phy.x)
    diag(coph_all) <- NA
    mpd_all <- mean(coph_all[lower.tri(coph_all)])
    
    # Initialize vector for species contributions to MPD
    sp_contrib <- vector()
    
    for (ii in 1:length(spp.ass)) {
      # MPD of the assemblage without species ii
      phy_min1 <- keep.tip(phy.x, spp.ass[spp.ass != spp.ass[ii]])
      coph_min1 <- ape::cophenetic.phylo(phy_min1)
      diag(coph_min1) <- NA
      mpd_min1 <- mean(coph_min1[lower.tri(coph_min1)])
      
      # Contribution of species ii to MPD
      sp_contrib[ii] <- (mpd_all - mpd_min1) / mpd_all
    }
    
    # create df of contribMPD
    contrib_mpd_df <- data.frame(Species = spp.ass,
                                 mpd_effect = sp_contrib)
    
    # Add contributions to the output
    ed.x <- left_join(ed.x, contrib_mpd_df, by = join_by(Species))
    
    # Combine results across iterations
    ed.df <- rbind(ed.df, ed.x)
  }  
  
  ed.mean <- ed.df %>%
    group_by(Species) %>% 
    dplyr::summarise(ed.mean     = mean(w),
                     ed.std.mean = mean(ed.std),
                     mpd_effect  = mean(mpd_effect)) %>%
    left_join(species.data_r %>% select(sp.name, species_ass), 
              by = c("Species" = "sp.name")) %>%
    ungroup()

  ed.mean  
  
}

## Add local ED and contribMPD to comm.data
comm.data$ed     <- species_evoldata$ed.mean[match(comm.data$species_ass, species_evoldata$species_ass)]
comm.data$ed.std <- species_evoldata$ed.std.mean[match(comm.data$species_ass, species_evoldata$species_ass)]
comm.data$contrib_mpd <- species_evoldata$mpd_effect[match(comm.data$species_ass, species_evoldata$species_ass)]

## Add local ED and contribMPD to species.data
species.data$ed.std <- species_evoldata$ed.std.mean[match(species.data$species_ass, species_evoldata$species_ass)]
species.data$contrib_mpd <- species_evoldata$mpd_effect[match(species.data$species_ass, species_evoldata$species_ass)]


# SAVE --------------------------------------------------------------------
write.csv(comm.data, "./data/data_pres_evol.csv", row.names = FALSE)
write.csv(species.data, "./data/data_species_evol.csv", row.names = FALSE)



