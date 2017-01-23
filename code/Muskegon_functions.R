###  Functions for Mukegon Analysis

# 1. make_metadata = create metadata from sample ID information in $names column
# 2. matround = a "real" rounding function
# 3. scale_reads = scale reads to equal sequencing depth  
# 4. make_metadata_norep = create metadata from sample ID information in $names column without a replicate
# 5. plot_sample_sums =  Draws a 3 paneled plot including: 1. histogram, 2. density, and 3. boxplot of the total number of sequences per sample.
# 6. muskegon_lakesite = fixes the name of the Muskegon lakesites (e.g. "MIN" = "Inlet")
# 7. show_otu = show the top 10 rows and columns of an otu table



#################################################################################### 1
#################################################################################### 1

## Written for this data by Marian Schmidt on July 21st, 2016
## This function adds the categorical metadata to a dataframe based on the sample name
# IMPORTANT!!!  The dataframe MUST have column named "names"
make_muskegon_metadata <- function(dataframe){ 
  
  # Create new columnes based on information in the sample name
  dataframe$lakesite <- substr(dataframe$names,1,3) # 1st 3 letters of string = lake site info
  dataframe$limnion <- substr(dataframe$names, 4, 4) # 4th letter = where in water column sample was taken
  dataframe$replicate <- substr(dataframe$names, 5, 5) # 5th letter = Which sample replicate 
  dataframe$fraction <- substr(dataframe$names, 6,6) # 6th letter = filter fraction (e.g. particle, whole, free)
  dataframe$month <- substr(dataframe$names, 7,7) # 7th letter = month sampled
  dataframe$year <- substr(dataframe$names, 8,9) # 8-9th letters = year sampled
  dataframe$nuc_acid_type <- substr(dataframe$names, 10,10) # 10th letter = nucleic acid sequenced (e.g. DNA,RNA,cDNA)
  
  # Let's make these letters meaningful to a human
  dataframe$lakesite <- as.factor(dataframe$lakesite)
  dataframe$lakesite <- factor(dataframe$lakesite,levels = c("MLB", "M15", "M45", "110", 
                                                             "MIN", "MBR", "MDP", "MOT", "EXT", "PBS", "MIL"))
  
  
  # Depth in water column
  dataframe$limnion <- ifelse(dataframe$limnion == "E", "Top", 
                              ifelse(dataframe$limnion == "S", "Top", 
                                     ifelse(dataframe$limnion == "D","Bottom",
                                            ifelse(dataframe$limnion == "H","Bottom",
                                                   ifelse(dataframe$limnion == "M","Middle",
                                                          ifelse(dataframe$limnion == "B", "Benthic",
                                                                 "MilliQ"))))))
  dataframe$limnion <- as.factor(dataframe$limnion)
  dataframe$limnion <- factor(dataframe$limnion,levels = c("Top", "Middle", "Bottom", "Benthic","MilliQ"))
  
  
  # Replicate
  dataframe$replicate <- ifelse(dataframe$replicate == "4", "2", 
                                ifelse(dataframe$replicate == "A", "1", 
                                       ifelse(dataframe$replicate == "Q","1",
                                              ifelse(dataframe$replicate == "I","1",
                                                     "1"))))  
  dataframe$replicate <- as.factor(dataframe$replicate)
  dataframe$replicate <- factor(dataframe$replicate,levels = c("1", "2"))
  
  
  # fraction Fraction
  dataframe$fraction <- ifelse(dataframe$fraction == "F", "Free", 
                               ifelse(dataframe$fraction == "P", "Particle", 
                                      ifelse(dataframe$fraction == "W","Whole",
                                             ifelse(dataframe$fraction == "C","Extraction_control",
                                                    ifelse(dataframe$fraction == "J","WholePart",
                                                           ifelse(dataframe$fraction == "K","WholeFree",
                                                                  ifelse(dataframe$fraction == "N","Nano",
                                                                         ifelse(dataframe$fraction == "S","Sediment",
                                                                                "PBS_Rinse"))))))))
  dataframe$fraction <- as.factor(dataframe$fraction)
  dataframe$fraction <- factor(dataframe$fraction,levels = c("Particle", "Free", "Whole", "Sediment","WholePart", "WholeFree", "Nano", "Extraction_control", "PBS_Rinse"))
  
  
  # Month
  dataframe$month <- ifelse(dataframe$month == "4", "April", 
                            ifelse(dataframe$month == "5", "May", 
                                   ifelse(dataframe$month == "6","June",
                                          ifelse(dataframe$month == "7","July",
                                                 ifelse(dataframe$month == "8","August",
                                                        ifelse(dataframe$month == "9","September",
                                                               ifelse(dataframe$month == "1","October",
                                                                      "NA")))))))
  dataframe$month <- as.factor(dataframe$month)
  dataframe$month <- factor(dataframe$month,levels = c("April", "May", "June", "July", "August", "September", "October"))
  
  
  #Year 
  dataframe$year <- ifelse(dataframe$year == "14", "2014", 
                           ifelse(dataframe$year == "15", "2015", 
                                  "NA")) 
  dataframe$year <- as.factor(dataframe$year)
  dataframe$year <- factor(dataframe$year,levels = c("2014", "2015"))
  
  # Nucleic Acid Type
  dataframe$nuc_acid_type <- ifelse(dataframe$nuc_acid_type == "D", "DNA",
                                    ifelse(dataframe$nuc_acid_type == "R", "RNA",
                                           ifelse(dataframe$nuc_acid_type == "C", "cDNA",
                                                  "DNA")))
  dataframe$nuc_acid_type <- as.factor(dataframe$nuc_acid_type)
  dataframe$nuc_acid_type <- factor(dataframe$nuc_acid_type,levels = c("DNA", "RNA", "cDNA"))
  
  
  
  #### What Project are the samples from??  -> Start with a new column named project from lakesite
  dataframe$project <- dataframe$lakesite
  # Change the information in the column to reflect which project
  dataframe$project <- ifelse(dataframe$project %in% c("MLB","M15","M45","110"), "Lake_Michigan", 
                              ifelse(dataframe$project %in% c("MBR","MDP","MIN","MOT"), "Muskegon_Lake", 
                                     ifelse(dataframe$project == "EXT","Extraction_Control",
                                            ifelse(dataframe$project == "MIL","MilliQ",
                                                   "PBS_Control"))))
  dataframe$project <- as.factor(dataframe$project)
  
  ## Add the Season  
  dataframe$season <- dataframe$month
  dataframe$season <- ifelse(dataframe$season == "April", "Spring", 
                             ifelse(dataframe$season == "May", "Spring", 
                                    ifelse(dataframe$season == "June", "Summer",
                                           ifelse(dataframe$season == "July", "Summer",
                                                  ifelse(dataframe$season == "August", "Summer",
                                                         ifelse(dataframe$season == "September", "Fall",
                                                                ifelse(dataframe$season == "October", "Fall",
                                                                       "NA")))))))
  dataframe$season <- as.factor(dataframe$season)
  dataframe$season <- factor(dataframe$season,levels = c("Spring", "Summer", "Fall"))
  
  # Add rownames as the names
  row.names(dataframe) <- dataframe$names
  
  # Give us the dataframe, please!
  return(dataframe)
  
}  
#################################################################################### 1
#################################################################################### 1




#################################################################################### 2
#################################################################################### 2
# This function was written by Michelle Berry , available at http://deneflab.github.io/MicrobeMiseq/
# Better rounding function than R's base round
matround <- function(x){trunc(x+0.5)}
#################################################################################### 2
#################################################################################### 2





#################################################################################### 3
#################################################################################### 3

# Modified from code written by Michelle Berry, available at http://deneflab.github.io/MicrobeMiseq/ 
# Scales reads by 
# 1) taking proportions
# 2) multiplying by a given library size of n
# 3) rounding 
# Default for n is the minimum sample size in your library
# Default for round is floor
scale_reads <- function(physeq, n = min(sample_sums(physeq)), round = "round") {
  
  # transform counts to n
  physeq.scale <- transform_sample_counts(physeq, 
                                          function(x) {(n * x/sum(x))}
  )
  
  # Pick the rounding functions
  if (round == "floor"){
    otu_table(physeq.scale) <- floor(otu_table(physeq.scale))
  } else if (round == "round"){
    otu_table(physeq.scale) <- round(otu_table(physeq.scale))
  } else if (round == "matround"){
    otu_table(physeq.scale) <- matround(otu_table(physeq.scale))
  }
  
  # Prune taxa and return new phyloseq object
  physeq.scale <- prune_taxa(taxa_sums(physeq.scale) > 0, physeq.scale)
  return(physeq.scale)
}

#################################################################################### 3
#################################################################################### 3




#################################################################################### 4
#################################################################################### 4

## Written for this data by Marian Schmidt on July 21st, 2016
## This function adds the categorical metadata to a dataframe based on the sample name
# IMPORTANT!!!  The dataframe MUST have column named "names" that DOES NOT have the REPLICATE
make_metadata_norep <- function(dataframe){ 
  
  # Create new columnes based on information in the sample name
  dataframe$lakesite <- substr(dataframe$names,1,3) # 1st 3 letters of string = lake site info
  dataframe$limnion <- substr(dataframe$names, 4, 4) # 4th letter = where in water column sample was taken
  ###  NO REPLICATE COLUMN
  dataframe$fraction <- substr(dataframe$names, 5,5) # 6th letter = filter fraction (e.g. particle, whole, free)
  dataframe$month <- substr(dataframe$names, 6,6) # 7th letter = month sampled
  dataframe$year <- substr(dataframe$names, 7,8) # 8-9th letters = year sampled
  dataframe$nuc_acid_type <- substr(dataframe$names, 10,10) # 10th letter = nucleic acid sequenced (e.g. DNA,RNA,cDNA)
  
  # Let's make these letters meaningful to a human
  dataframe$lakesite <- as.factor(dataframe$lakesite)
  dataframe$lakesite <- factor(dataframe$lakesite,levels = c("MLB", "M15", "M45", "110", 
                                                             "MIN", "MBR", "MDP", "MOT", "EXT", "PBS", "MIL"))
  
  
  # Depth in water column
  dataframe$limnion <- ifelse(dataframe$limnion == "E", "Top", 
                              ifelse(dataframe$limnion == "S", "Top", 
                                     ifelse(dataframe$limnion == "D","Bottom",
                                            ifelse(dataframe$limnion == "H","Bottom",
                                                   ifelse(dataframe$limnion == "M","Middle",
                                                          ifelse(dataframe$limnion == "B", "Benthic",
                                                                 "MilliQ"))))))
  dataframe$limnion <- as.factor(dataframe$limnion)
  dataframe$limnion <- factor(dataframe$limnion,levels = c("Top", "Middle", "Bottom", "Benthic","MilliQ"))
  
  
  # fraction Fraction
  dataframe$fraction <- ifelse(dataframe$fraction == "F", "Free", 
                               ifelse(dataframe$fraction == "P", "Particle", 
                                      ifelse(dataframe$fraction == "W","Whole",
                                             ifelse(dataframe$fraction == "C","Extraction_control",
                                                    ifelse(dataframe$fraction == "J","WholePart",
                                                           ifelse(dataframe$fraction == "K","WholeFree",
                                                                  ifelse(dataframe$fraction == "N","Nano",
                                                                         ifelse(dataframe$fraction == "S","Sediment",
                                                                                "PBS_Rinse"))))))))
  dataframe$fraction <- as.factor(dataframe$fraction)
  dataframe$fraction <- factor(dataframe$fraction,levels = c("Particle", "Free", "Whole", "Sediment","WholePart", "WholeFree", "Nano", "Extraction_control", "PBS_Rinse"))
  
  
  # Month
  dataframe$month <- ifelse(dataframe$month == "4", "April", 
                            ifelse(dataframe$month == "5", "May", 
                                   ifelse(dataframe$month == "6","June",
                                          ifelse(dataframe$month == "7","July",
                                                 ifelse(dataframe$month == "8","August",
                                                        ifelse(dataframe$month == "9","September",
                                                               ifelse(dataframe$month == "1","October",
                                                                      "NA")))))))
  dataframe$month <- as.factor(dataframe$month)
  dataframe$month <- factor(dataframe$month,levels = c("April", "May", "June", "July", "August", "September", "October"))
  
  
  #Year 
  dataframe$year <- ifelse(dataframe$year == "14", "2014", 
                           ifelse(dataframe$year == "15", "2015", 
                                  "NA")) 
  dataframe$year <- as.factor(dataframe$year)
  dataframe$year <- factor(dataframe$year,levels = c("2014", "2015"))
  
  # Nucleic Acid Type
  dataframe$nuc_acid_type <- ifelse(dataframe$nuc_acid_type == "D", "DNA",
                                    ifelse(dataframe$nuc_acid_type == "R", "RNA",
                                           ifelse(dataframe$nuc_acid_type == "C", "cDNA",
                                                  "DNA")))
  dataframe$nuc_acid_type <- as.factor(dataframe$nuc_acid_type)
  dataframe$nuc_acid_type <- factor(dataframe$nuc_acid_type,levels = c("DNA", "RNA", "cDNA"))
  
  
  
  #### What Project are the samples from??  -> Start with a new column named project from lakesite
  dataframe$project <- dataframe$lakesite
  # Change the information in the column to reflect which project
  dataframe$project <- ifelse(dataframe$project %in% c("MLB","M15","M45","110"), "Lake_Michigan", 
                              ifelse(dataframe$project %in% c("MBR","MDP","MIN","MOT"), "Muskegon_Lake", 
                                     ifelse(dataframe$project == "EXT","Extraction_Control",
                                            ifelse(dataframe$project == "MIL","MilliQ",
                                                   "PBS_Control"))))
  dataframe$project <- as.factor(dataframe$project)
  
  ## Add the Season  
  dataframe$season <- dataframe$month
  dataframe$season <- ifelse(dataframe$season == "April", "Spring", 
                             ifelse(dataframe$season == "May", "Spring", 
                                    ifelse(dataframe$season == "June", "Summer",
                                           ifelse(dataframe$season == "July", "Summer",
                                                  ifelse(dataframe$season == "August", "Summer",
                                                         ifelse(dataframe$season == "September", "Fall",
                                                                ifelse(dataframe$season == "October", "Fall",
                                                                       "NA")))))))
  dataframe$season <- as.factor(dataframe$season)
  dataframe$season <- factor(dataframe$season,levels = c("Spring", "Summer", "Fall"))
  
  # Add rownames as the names
  row.names(dataframe) <- dataframe$names
  
  # Give us the dataframe, please!
  return(dataframe)
  
}  

#################################################################################### 4
#################################################################################### 4



#################################################################################### 5
#################################################################################### 5

# This function will make a 3 paneled plot including: 
# A. histogram, B. density, and C. boxplot of the total number of sequences per sample.

# Arguments:
# dataframe = a dataframe with information
# x_total_seqs = a column in dataframe with the total number of sequences per sample
# fill variable = a categorical variable to fill in the boxplot, hostogram and density plots by
# brewer_pal = a color brewer pallette to color the plot by


plot_sample_sums <- function(dataframe, x_total_seqs, fill_variable, brewer_pal = "Set1"){
  
  # Histogram
  year_hist <- ggplot(dataframe, aes_string(x = x_total_seqs, fill = fill_variable)) +
    xlab("# of Seqs/Sample") + ylab("Count") +
    geom_histogram(position = "dodge",binwidth = 2000) +
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
    scale_fill_brewer(palette=brewer_pal) + 
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  
  
  # Density Plot
  year_dens <- ggplot(dataframe, aes_string(x = x_total_seqs, fill = fill_variable)) +
    xlab("# of Seqs/Sample") + ylab("Density") +
    geom_density(alpha = 0.5) +
    scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
    scale_fill_brewer(palette=brewer_pal) + 
    theme(legend.position = "none", 
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  
  # Boxplot
  year_box <- ggplot(dataframe, aes_string(y = x_total_seqs, x = fill_variable, fill = fill_variable)) + 
    xlab(fill_variable) + ylab("# of Seqs/Sample") +
    geom_boxplot(color = "black") + 
    scale_fill_brewer(palette=brewer_pal) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  
  
  return(ggdraw() +
            draw_plot(year_hist, 0, 0, 0.3, 1) +
            draw_plot(year_dens, 0.28, 0, 0.33, 1) + # 1st = where to start drawing, 3rd = width, 4th = height
            draw_plot(year_box, 0.6, 0, 0.4, 1) + # 3rd = width, 4th = height
            draw_plot_label(c("A", "B", "C"), 
                            c(0, 0.3, 0.58), # Where along the x-axis would you like the labels?
                            c(1, 1, 1),  # Put the label at the top of the plotting space (y-axis)
                            size = 15) # Size of label
  )
}

#################################################################################### 5
#################################################################################### 5





#################################################################################### 6
#################################################################################### 6
# Make a new column called "station" with the human name of the Muskegon sampling stations
muskegon_lakesite <- function(dataframe){
  
  ## Meaningful lakesites to a human 
  dataframe$station <- dataframe$lakesite
  dataframe$station <- ifelse(dataframe$lakesite == "MIN", "River", 
                              ifelse(dataframe$lakesite == "MBR", "Bear", 
                                     ifelse(dataframe$lakesite == "MDP", "Deep",
                                            ifelse(dataframe$lakesite == "MOT", "Channel",
                                                   ifelse(dataframe$lakesite == "PBS", "Control",
                                                          ifelse(dataframe$lakesite == "MIL", "Control",
                                                                 "NA"))))))
  # Give us the dataframe, please!
  return(dataframe)
}

#################################################################################### 6
#################################################################################### 6 


#################################################################################### 7
#################################################################################### 7
# This function was written by Marian Schmidt
# Show the top of an OTU table
show_otu <- function(dataframe, nrow = 10, ncol = 10){
  dataframe[1:nrow,1:ncol]
}
#################################################################################### 7
#################################################################################### 7

