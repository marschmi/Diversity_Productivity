## Purpose: Create the .RData for the raw OTU phyloseq object
## Author: Marian L. Schmidt 

# Set up file
library(phyloseq)
library(dplyr)
source("Muskegon_functions.R")


############################################################################################################################################
############################################################################################################################################
###################################################################### LOAD IN THE DATA 
#  Load in the taxonomy and shared data
otu_tax <- "../data/mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy"
otu_shared <- "../data/mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared"
otu_data <- import_mothur(mothur_shared_file = otu_shared, mothur_constaxonomy_file = otu_tax)

# Fix the taxonomy names
colnames(tax_table(otu_data)) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

############################################################################################################################################
############################################################################################################################################
###################################################################### ADD THE PROTEOBACTERIA TO THE PHYLA
phy <- data.frame(tax_table(otu_data))
Phylum <- as.character(phy$Phylum)
Class <- as.character(phy$Class)

for  (i in 1:length(Phylum)){ 
  if (Phylum[i] == "Proteobacteria"){
    Phylum[i] <- Class[i]
  } 
}

phy$Phylum <- Phylum # Add the new phylum level data back to phy
phy$OTU <- row.names(tax_table(otu_data)) # Make a column for OTU

t <- tax_table(as.matrix(phy))

tax_table(otu_data) <- t
################################################

# Sample Names
samp_names <- colnames(otu_table(otu_data))

# Create metadata info
df <- data.frame(matrix(NA, ncol = 1, nrow = length(samp_names)))
colnames(df) <- c("names")
df$names <- samp_names


############################################################################################################################################
############################################################################################################################################
###################################################################### LOAD IN METADATA
# Create metadata info
meta_df <- make_muskegon_metadata(df)

# Create a norep_water_name column
meta_df$norep_water_name <- paste(substr(meta_df$names,1,4),substr(meta_df$names,7,9), sep = "")

# Load in the extra metadata for the Muskegon Lake Project
musk_data <- read.csv("../data/metadata/processed_muskegon_metadata.csv", header = TRUE) # Load in the extra metada
musk_data_subset <- select(musk_data, -c(lakesite, limnion, month, year, project, season))
complete_meta_df <- left_join(meta_df, musk_data_subset, by = "norep_water_name")
row.names(complete_meta_df) <- complete_meta_df$names

complete_meta_df$water_name <- paste(substr(complete_meta_df$names,1,5),substr(complete_meta_df$names,7,9), sep = "")
complete_meta_df$norep_filter_name <- paste(substr(complete_meta_df$names,1,4),substr(complete_meta_df$names,6,9), sep = "")


############################################################################################################################################
############################################################################################################################################
###################################################################### MERGE REPLICATE SAMPLES 
# Add the metadata to our phyloseq object! 
##  Add metadata to OTUs
sample_data(otu_data) <- complete_meta_df
otu_data <- prune_taxa(taxa_sums(otu_data) > 0, otu_data)


# Merged metadata
merged_complete_meta_df<- 
  select(complete_meta_df, -c(names, replicate, nuc_acid_type, water_name)) %>% 
  distinct() %>%
  arrange(norep_filter_name)

## Add Production data from GVSU AWRI provided by the lab of Bopi Biddanda
bopi_data <- read.csv("../data/metadata/production_data.csv") %>%
  dplyr::rename(norep_filter_name = names)  %>% # rename "names" column to "norep_filter_name"
  select(-c(X, limnion, fraction, month, year, season))

## Merge two metadata files together!
df1 <- full_join(merged_complete_meta_df, bopi_data) %>%
  select(-c(Depth, month)) 

# provide row.names to match sample
row.names(df1) <- df1$norep_filter_name



############################################################################################################################################
############################################################################################################################################
###################################################################### LOAD IN FLOW CYTOMETRY DATA

flow_cy <- read.csv2("../data/metadata/flow_cytometry.csv", header = TRUE) %>%
  select(-c(X, Season, Month, Year, Fraction, Site, Depth)) %>%
  filter(Lake == "Muskegon") %>%
  # Rename columns so they are more representative of the data within them
  rename(volume_uL = volume, 
         raw_counts = counts,
         HNA_counts = HNA,
         LNA_counts = LNA) %>%
  mutate(cells_per_uL = (raw_counts*dilution)/volume_uL,
         HNA_per_uL = (HNA_counts*dilution)/volume_uL,
         HNA_percent = (HNA_counts)/raw_counts*100,
         LNA_per_uL = (LNA_counts*dilution)/volume_uL,
         LNA_percent = (LNA_counts)/raw_counts,
         Nuc_acid_ratio = HNA_per_uL/LNA_per_uL) %>%
  select(-Lake)

# Create a new column so we can merge with other data frames
flow_cy$norep_water_name <- paste(substring(flow_cy$Sample_16S,1,4), substring(flow_cy$Sample_16S,7,9), sep = "")

# Join original metadata and flow cytometry metadata
df2 <- left_join(df1, flow_cy, by = "norep_water_name")
row.names(df2) <- row.names(df1)

# Fix the factor levels for nice plotting
df2$lakesite <- factor(df2$lakesite,  levels = c("MOT", "MDP", "MBR", "MIN"))
df2$station <- factor(df2$station,  levels = c("Channel", "Deep", "Bear", "River"))
df2$season <- factor(df2$season,  levels = c("Spring", "Summer", "Fall"))
df2$year <- factor(df2$year,  levels = c("2014", "2015"))


# merge samples for OTUs
otu_merged <- merge_samples(otu_data, group = "norep_filter_name", fun = "sum")

# Add nice sample information
sample_data(otu_merged) <- df2





############################################################################################################################################
############################################################################################################################################
###################################################################### SUBSET MUSKEGON LAKE SAMPLES 
# Subset Mukegon Lake samples out of the OTU phyloseq object
otu_merged_musk <- subset_samples(physeq = otu_merged, project == "Muskegon_Lake")
otu_merged_musk_pruned <- prune_taxa(taxa_sums(otu_merged_musk) > 0, otu_merged_musk) 

# Remove the flow cytometry data from the particle samples!
df3 <- sample_data(otu_merged_musk_pruned)


# # Subset the particle samples 
PA_df3 <- filter(df3, fraction %in% c("WholePart", "Particle"))
flow_cy_columns <- colnames(flow_cy)

##  Make all of the clow cytometry columns NA
for(i in flow_cy_columns){
  PA_df3[,i] <- NA
}

# Subset the FL samples
FL_df3 <- filter(df3, fraction %in% c("WholeFree", "Free", "Sediment"))

## Recombine the dataframe back together with flow cytometry information as NA for Particle Samples
final_meta_dataframe <- bind_rows(FL_df3, PA_df3)
row.names(final_meta_dataframe) <- final_meta_dataframe$norep_filter_name

sample_data(otu_merged_musk_pruned) <- final_meta_dataframe

# Create a new file called "Phyloseq.RData" that has the phyloseq object
save(list="otu_merged_musk_pruned", file=paste0("../data/otu_merged_musk_pruned.RData")) 

