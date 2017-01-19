# Why this script?  

##  OUTPUT:  data.frame that has calculated Hill Diversity metrics with Phenoflow package
### 1. iNEXT_oligo.RData = iNEXT and estimateD data object with oligotyping table

## INPUT: Phyloseq object with otu_table that has rows as samples and columns as taxa

# WHERE:  Run on flux in the shell
# Run the script in the same working directory where the data files are available 

# TO EXECUTE:
# R make_phenoflow_diversity.R path/to/phyloseq_object  phyloseq_name & > log.txt
      ### 2nd argument (phyloseq_name) must be what the phyloseq name was while it was being saved (variable name)

################################################################################## SETUP
##################################################################################
## Load the R-packages necessary for the script below to run
library(Phenoflow) # iNEXT package for alpha diversity measurements 

###  Set the user preferences  
userprefs <- commandArgs(trailingOnly = TRUE)
phyloseq_object <- userprefs[1]
phyloseq_name <- userprefs[2]

# Write an update in the commandline that the script is working...
cat(paste0(date(),"\tCalculating Phenoflow Hill Diversity Estimates with phyloseq data\n"))


################################################################################## iNEXT ON OLIGOTYPES 
### Make Oligotyping iNEXT object
# Read in the oligotyping file to R as a data frame, assumes a tab delimited file 
load(phyloseq_object) 

# Calculate Oligotype richness, shannon and simpson Hill Diversity estimates with iNEXT function 
phenoflow_diversity <- Diversity_16S(phyloseq_name, R=100, brea = FALSE) # 100 bootstraps and do not do breakaway analysis


################################################################################## FINALIZE
# Finalize script
cat(paste0(date(),"\tExporting Phenoflow Hill Diversity Estimates to Phenoflow_diversity.RData\n")) # Update user that .Rdata is being exported

# Create a new file called "iNEXT_oligo.RData" that has the oligotyping inext object
save("phenoflow_diversity", file=paste0("Phenoflow_diversity.RData")) 

# Script has finished!
cat(paste0(date(),"\tDone! Time to work on your alpha diversity analysis!\n")) 


# SOURCES:
# Props R, Monsieurs P, Mysara M, Clement L and Boon N (2016). _Phenoflow_. R package version 1.0, <URL:https://github.com/rprops/Phenoflow_package>.

# fin