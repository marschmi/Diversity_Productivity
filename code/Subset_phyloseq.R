# Prepare only surface and WholeFree/WholePart fractions for analysis
# Marian L Schmidt
# May 26, 2017



############################# LOAD LIBRARIES ############################### 
library(ggplot2)
library(devtools)
library(phyloseq)
library(tidyr)
library(dplyr)
library(cowplot)
library(picante) # Will also include ape and vegan 
source("code/Muskegon_functions.R")
source("code/set_colors.R")


############################# LOAD DATA ############################### 
# Loads a phyloseq object named otu_merged_musk_pruned)
load("data/otu_merged_musk_pruned.RData")
# The name of the phyloseq object is: 
otu_merged_musk_pruned 


############################# PRUNE DATA ############################### 
# Remove all samples except the ones used in this study! 
surface_PAFL_otu <- subset_samples(otu_merged_musk_pruned, fraction %in% c("WholeFree", "WholePart") & limnion == "Top")

## Add total sequences to metadata frame 
metdf <- sample_data(surface_PAFL_otu) %>%
  dplyr::select(-c(42:67, norep_water_name, dnaconcrep2)) %>%
  # Productivity measurements are reliable only up to 1 decimal
  dplyr::mutate(tot_bacprod = round(tot_bacprod, digits = 1),
                SD_tot_bacprod = round(SD_tot_bacprod, digits = 1),
                frac_bacprod = round(frac_bacprod, digits = 1),
                SD_frac_bacprod = round(SD_frac_bacprod, digits = 1),
                fraction_bac_abund = as.numeric(fraction_bac_abund),
                fracprod_per_cell = frac_bacprod/(1000*fraction_bac_abund),
                fracprod_per_cell_noinf = ifelse(fracprod_per_cell == Inf, NA, fracprod_per_cell))

row.names(metdf) = metdf$norep_filter_name
# Rename the sample data 
sample_data(surface_PAFL_otu) <- metdf

# Remove taxa that are no longer in the phyloseq object
surface_PAFL_otu_pruned_raw <- prune_taxa(taxa_sums(surface_PAFL_otu) > 0, surface_PAFL_otu) 
surface_PAFL_otu_pruned_raw

############################# WRITE OUT DATA IN DATA FOLDER ############################### 
# Create a new file called "surface_PAFL_otu_pruned_RAREFIED_rm2.RData" that has the phyloseq object
save("surface_PAFL_otu_pruned_raw", file=paste0("data/surface_PAFL_otu_pruned_raw.RData")) 
