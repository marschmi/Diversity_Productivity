## Purpose: Create the .RData for the raw OTU phyloseq object
## Author: Marian L. Schmidt 

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
musk_data <- read.csv("data/metadata/processed_muskegon_metadata.csv", header = TRUE) # Load in the extra metada
musk_data_subset <- select(musk_data, -c(lakesite, limnion, month, year, project, season))
complete_meta_df <- left_join(meta_df, musk_data_subset, by = "norep_water_name")
row.names(complete_meta_df) <- complete_meta_df$names

complete_meta_df$water_name <- paste(substr(complete_meta_df$names,1,5),substr(complete_meta_df$names,7,9), sep = "")
complete_meta_df$norep_filter_name <- paste(substr(complete_meta_df$names,1,4),substr(complete_meta_df$names,6,9), sep = "")



