# Prepare Final tree for Phylogenetic Analysis
# Marian L Schmidt
# May 25, 2017


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

# Remove taxa that are no longer in the phyloseq object
surface_PAFL_otu_pruned <- prune_taxa(taxa_sums(surface_PAFL_otu) > 0, surface_PAFL_otu) 
surface_PAFL_otu_pruned

############################# VERIFY DATA ############################### 
sums_otu <- data.frame(rowSums(otu_table(surface_PAFL_otu_pruned)))
colnames(sums_otu) <- "Sample_TotalSeqs"
sums_otu$names <- row.names(sums_otu)
sums_otu <- arrange(sums_otu, Sample_TotalSeqs) 
sums_otu <- make_metadata_norep(sums_otu)

plot_sample_sums(dataframe = sums_otu, x_total_seqs = "Sample_TotalSeqs", fill_variable = "fraction")


############################# RAREFY DATA ############################### 
# Since we are going to create a tree, we first need to rarefy the samples
min_lib <- min(sample_sums(surface_PAFL_otu_pruned)) - 1
min_lib

# Set the seed for the randomization of rarefy-ing
set.seed(777)

surface_PAFL_otu_pruned_RAREFIED <- rarefy_even_depth(surface_PAFL_otu_pruned, sample.size = min_lib,
                                                       verbose = FALSE, replace = TRUE)

############################# PRUNE OUT DOUBLETONS ############################### 
###  REMOVE ALL OTUS THAT HAVE LESS THAN 2 ACROSS THE ENTIRE DATASET
surface_PAFL_otu_pruned_RAREFIED_rm2 <- prune_taxa(taxa_sums(surface_PAFL_otu_pruned_RAREFIED) > 2, surface_PAFL_otu_pruned_RAREFIED) 
surface_PAFL_otu_pruned_RAREFIED_rm2



############################# WRITE OUT DATA IN PHYLOTREE FOLDER ############################### 
# Create a new file called "surface_PAFL_otu_pruned_RAREFIED_rm2.RData" that has the phyloseq object
save("surface_PAFL_otu_pruned_RAREFIED_rm2", file=paste0("data/PhyloTree/surface_PAFL_otu_pruned_RAREFIED_rm2.RData")) 

# Obtain the OTU names that were retained in the dataset
tax <- data.frame(tax_table(surface_PAFL_otu_pruned_RAREFIED_rm2))
vector_of_otus <- as.vector(tax$OTU)

write(vector_of_otus, file = "data/PhyloTree/OTUnames_rarefied_rm2.txt",
      ncolumns = 1,
      append = FALSE, sep = "\n")
