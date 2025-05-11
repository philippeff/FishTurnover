# Warm-affinity and Locally Distinctive Species Dominate Turnover in Marine Fish Assemblages

This repository contains R scripts, shell files and data files for a chapter of my PhD thesis

## Abstract
Marine fish assemblages are rapidly changing as a response to human-mediated changes in ocean temperatures. While poleward shifts in fish species ranges driven by ocean warming are becoming better documented, less is known about how these changes are reshaping the evolutionary composition of local assemblages. Using a global database of time series data (BioTIME), we examined how local fish assemblages are changing in terms of species' thermal affinity and evolutionary relatedness. We found that warm-water species are increasingly colonizing local assemblages, a trend that occurs across all ocean temperatures and is amplified in regions experiencing faster warming. This shift results in widespread increases in average thermal affinity (community temperature index; CTI), particularly in temperate waters, highlighting a strong assemblage-level response to rising sea surface temperatures. We concurrently detect consistent shifts in evolutionary relatedness within assemblages: colonizing species tend to be more evolutionarily distinct locally, contributing to a general increase in phylogenetic mean pairwise distance (MPD), i.e. a decrease in average relatedness. These concurrent patterns in thermal affinity and evolutionary relatedness suggest that the tropicalization of local assemblages as a consequence of warming oceans is promoting the colonization of species that are not only associated to warmer waters but also more distant relatives on average. Our findings provide the first evidence at a global scale that local species turnover under climate change is not only altering species identities but also reshaping the evolutionary composition of marine fish assemblages.


---

## Repository Structure

### Data

- `data/data_assemb.csv` — Assemblage-level data on change in diversity metrics
- `data/data_assemb_final.csv` — Assemblage-level data on change in diversity metrics (final)
- `data/data_species.csv` — Species-level data
- `data/data_species_final.csv` — Species-level data (final)
- `data/data_time_final.csv` — yearly diversity metric data per assemblage (final)

### Data Preparation

- `01a_evoldata_cluster.R` — Prepares species-level evolutionary distinctiveness data.
- `01b_GBIF.R` — Downloads and processes GBIF occurrence records.
- `01c_SST.R` — Extracts sea surface temperature (SST).

### Temporal trends in diversity metrics

These scripts fit Bayesian hierarchical models to different diversity metrics:

- `02a_BRM_MPD_time_cluster.R` — Mean Pairwise Distance (MPD)
- `02b_BRM_sesMPD_time_cluster.R` — Standardized effect size MPD (ses.MPD)
- `02c_BRM_Jaccard_time_cluster.R` — Jaccard similarity index
- `02d_BRM_PhyloSor_time_cluster.R` — Phylogenetic Sørensen index
- `02e_BRM_MNTD_time_cluster.R` — Mean Nearest Taxon Distance (MNTD)
- `02f_BRM_SR_time_cluster.R` — Species Richness
- `02g_BRM_CTI_time_cluster.R` — Community Temperature Index (CTI)

### Post-Processing and Synthesis

- `03_time_posteriors.R` — Summarizes posterior draws of temporal trends models.
- `04a_GLM_turnover_STI.R` — GLM relating turnover to Species Temperature Index, sea surface temperature and latitude.
- `04b_GLM_turnover_ED.R` — GLM relating turnover to local Evolutionary Distinctiveness, sea surface temperature and latitude.
- `05_analysis.R` — Final analysis combining model outputs and plots.

### Job Submission Scripts

Shell scripts for submitting cluster jobs:

- `submit_01a_evoldata_cluster.sh`
- `submit_02a_BRM_MPD_time_cluster.sh`
- `submit_02b_BRM_sesMPD_time_cluster.sh`
- `submit_02c_BRM_Jaccard_time_cluster.sh`
- `submit_02d_BRM_PhyloSor_time_cluster.sh`
- `submit_02e_BRM_MNTD_time_cluster.sh`
- `submit_02f_BRM_SR_time_cluster.sh`
- `submit_02g_CTI_time_cluster.sh`

---

## Data Availability

This repository does **not** include all necessary data due to copyright and access restrictions:

- **Community composition data**: `data_pres_final.csv` is based on BioTIME data and cannot be redistributed here. If you wish to use this dataset, please contact **Philippe Fernandez-Fournier**.
- **Phylogenetic tree**: `phylo_fish_PFF.tre` contains curated phylogenetic data also not publicly archived due to licensing constraints. You may request it from the author.
- **Sea surface temperature data**: `era5_sst.nc` can be downloaded from the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/).

---

