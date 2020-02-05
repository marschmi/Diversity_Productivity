###  Functions for Mukegon Analysis

# 1. make_metadata = create metadata from sample ID information in $names column
# 2. matround = a "real" rounding function
# 3. scale_reads = scale reads to equal sequencing depth  
# 4. make_metadata_norep = create metadata from sample ID information in $names column without a replicate
# 5. plot_sample_sums =  Draws a 3 paneled plot including: 1. histogram, 2. density, and 3. boxplot of the total number of sequences per sample.
# 6. muskegon_lakesite = fixes the name of the Muskegon lakesites (e.g. "MIN" = "Inlet")
# 7. show_otu = show the top 10 rows and columns of an otu table  
# 8. grid_arrange_shared_legend = Create one shared legend in a multiplot figure



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


#################################################################################### 8
#################################################################################### 8
# This function was originally written by Ruben Props and further modified by Marian

# This function takes in a linear model output and a dataframe to do an analysis of 
# the residual and predicted values.
      # INPUT
      # 1. Linear model object from lm() function
      # 2. vector of observed y values (y-variable)
      # 3. Title in quotes for figure

      # OUTPUT
      # 3 plots:
      # Plot 1: QQPlot: Theoretical quantiles (x-axis) vs studentized residuals (y-axis) 
      # Plot 2:  Fitted values (x-axis) vs studentized residuals (y-axis) 
      # Plot 3:  Fitted values (x-axis) vs observed values (y-axis) 

plot_residuals <- function(lm_model, lm_observed_y, main_title){
  
  # set seed
  set.seed(777)
  
  # Make a 3 paneled plot
  par(mfrow=c(1,3), oma=c(0,0,0,0)) # c(bottom, left, top, right) 
  
  # 1st plot:  A qqplot; Theoretical quantiles (x-axis) vs studentized residuals (y-axis) 
  qqPlot(lm_model, col="blue", reps=10000, 
         ylab="Studentized residuals", 
         xlab="Theoretical quantiles (t-distribution)",
         cex=1.5, las=1)
  
  # 2nd plot:  Fitted values (x-axis) vs studentized residuals (y-axis) 
  plot(y=studres(lm_model), x=predict(lm_model), col="blue",
       las=1,ylab="Studentized residuals", xlab="Fitted Values", cex=1.5)
  # Draw a line at y=0
  lines(x=c(-10,60), y=c(0,0), lty=2)
  
  # If you studentized residuals is greater than 3, state that there's an outlier
  ifelse(sum(studres(lm_model) > 3) >= 1, 
         # Do this if the length of the studentized residuals is equal to or greater than 1
         print(paste("WARNING:You have", 
               length(studres(lm_model)[studres(lm_model) > 3]), 
               "high-leverage point(s)!")), 
         # Do this if the length of the studentized residuals is equal to or greater than 1
         print("There are no high leverage points in this model.")) 
  
  # Add a title to all the plots
  mtext(text = main_title, side = 3, line = 2, outer = FALSE)
  #  side:  1=bottom, 2=left, 3=top, 4=right
  
  # 3rd plot: Fitted values (x-axis) vs observed values (y-axis) 
  plot(y=lm_observed_y, x=predict(lm_model), col="blue",
       ylab="Observed values",
       xlab="Fitted Values", cex=1.5,
       las=1)
  # draw a 1:1 line
  abline(a = 0, b = 1, lty=2)
  

 }


#################################################################################### 8
#################################################################################### 8



#################################################################################### 9
#################################################################################### 9
# This function is from: http://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots/28594060#28594060
# Downloaded on March 3rd, 2017

grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}


#################################################################################### 9
#################################################################################### 9


#################################################################################### 10
#################################################################################### 10
# Calculate Alpha Diversity Values and combine into dataframe 


calc_mean_alphadiv <- function(physeq, richness_df, simpson_df, exp_shannon_df){
  
  # How many samples?
  nsamp <- nsamples(physeq)
  
  # Check that all sample names are the same
  stopifnot(row.names(richness_df) == row.names(simpson_df))
  stopifnot(row.names(richness_df) == row.names(exp_shannon_df))
  
  # Create sample names for the dataframes
  norep_filter_name <- row.names(richness_df)
  
  # Create a new dataframe to hold the means and standard deviations of richness estimates
  mean <- apply(richness_df, 1, mean)
  sd <- apply(richness_df, 1, sd)
  measure <- rep("Richness", nsamp)
  otu_rich_stats <- data.frame(norep_filter_name, mean, sd, measure)
  
  # Create a new dataframe to hold the means and standard deviations of evenness estimates
  mean <- apply(simpson_df, 1, mean)
  sd <- apply(simpson_df, 1, sd)
  measure <- rep("Inverse_Simpson", nsamp)
  otu_simps_stats <- data.frame(norep_filter_name, mean, sd, measure)
  
  # Create a new dataframe to hold the means and standard deviations of shannon entropy estimates
  mean <- apply(exp_shannon_df, 1, mean)
  sd <- apply(exp_shannon_df, 1, sd)
  measure <- rep("Exponential_Shannon", nsamp)
  otu_shan_stats <- data.frame(norep_filter_name, mean, sd, measure)
  
  # Combine alpha diversity into one dataframe 
  otu_alpha <- rbind(otu_rich_stats, otu_simps_stats, otu_shan_stats)
  s <- data.frame(sample_data(physeq))
  otu_alphadiv <- merge(otu_alpha, s, by = "norep_filter_name")
  
  # Give us the dataframe, please!
  return(otu_alphadiv)
}
#################################################################################### 10
#################################################################################### 10




#################################################################################### 11
#################################################################################### 11
# Take a phyloseq object and:
    # 1. Rarefy at the minimum sequencing depth -1
    # 2. Calculate (1) the observed richness, (2) inverse simpson, (3) Shannon entropy.
    # 3. Output a list of 3 matrices with data from #2

# Modified from original code written by Michelle Berry.  

# INPUT:
    # A phyloseq object

# OUTPUT:
    # A named list of three matrices where:
          # Matrix 1  = Richness
          # Matrix 2  = Inverse Simpson
          # Matrix 3  = Exp(Shannon Entropy)

calc_alpha_diversity <- function(physeq, iters = 100){

 nsamp <- nsamples(physeq)
 min_seqs <- min(sample_sums(physeq)) - 1
 
 richness <- matrix(nrow = nsamp, ncol = iters)
 row.names(richness) <- sample_names(physeq)
 
 simpson <- matrix(nrow = nsamp, ncol = iters)
 row.names(simpson) <- sample_names(physeq)
 
 shannon <- matrix(nrow = nsamp, ncol = iters)
 row.names(shannon) <- sample_names(physeq)

# Set the seed to have reproducible results
set.seed(777)

for (i in 1:iters) {
    
    # Subsample
    r <- rarefy_even_depth(physeq, sample.size = min_seqs, verbose = FALSE, replace = TRUE)
   
    # Calculate richness; hill = 0
    calc_rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
    richness[ ,i] <- calc_rich
    
    # Calculate the Shannon Exponential; hill = 1 
    calc_shan <- as.numeric(as.matrix(exp(estimate_richness(r, measures = "Shannon"))))
    shannon[ ,i] <- calc_shan
    
    # Calculate the Inverse Simpson; hill = 2
    calc_simps <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
    simpson[ ,i] <- calc_simps

  }

  # Return a list of 3 named matrices 
  return(list(Richness = richness, Inverse_Simpson = simpson, Shannon = shannon))

}
#################################################################################### 11
#################################################################################### 11

# The input data frame must have:
        # A column named fraction
        # A column named Removed
        # Already be subsetted for a specific diversity measure s 

lm_fraction_output <- function(dataframe) {
  
  # Check if factor level order is correct!
  removed_factors <- unique(dataframe$Removed)
  removed_factor_order <- c("1-tons", "5-tons", "10-tons", "20-tons", "30-tons", "60-tons", "90-tons", "150-tons", "225-tons", "300-tons")
  stopifnot(removed_factor_order == removed_factors) # Stop if it is incorrect!
  
  lm_01_part_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe,Removed == "1-tons" & fraction == "Particle"))
  lm_05_part_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "5-tons" & fraction == "Particle"))
  lm_10_part_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "10-tons" & fraction == "Particle"))
  lm_20_part_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "20-tons" & fraction == "Particle"))
  lm_30_part_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "30-tons" & fraction == "Particle"))
  lm_60_part_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "60-tons" & fraction == "Particle"))
  lm_90_part_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "90-tons" & fraction == "Particle"))
  lm_150_part_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "150-tons" & fraction == "Particle"))
  lm_225_part_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "225-tons" & fraction == "Particle"))
  lm_300_part_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "300-tons" & fraction == "Particle"))
  
  part_lm_adj_r2_comm <- c(
    round(summary(lm_01_part_comm)$adj.r.squared, digits = 2),
    round(summary(lm_05_part_comm)$adj.r.squared, digits = 2),
    round(summary(lm_10_part_comm)$adj.r.squared, digits = 2),
    round(summary(lm_20_part_comm)$adj.r.squared, digits = 2),
    round(summary(lm_30_part_comm)$adj.r.squared, digits = 2),
    round(summary(lm_60_part_comm)$adj.r.squared, digits = 2),
    round(summary(lm_90_part_comm)$adj.r.squared, digits = 2),
    round(summary(lm_150_part_comm)$adj.r.squared, digits = 2),
    round(summary(lm_225_part_comm)$adj.r.squared, digits = 2),
    round(summary(lm_300_part_comm)$adj.r.squared, digits = 2))
  
  part_lm_pval_comm <- c(
    round(unname(summary(lm_01_part_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_05_part_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_10_part_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_20_part_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_30_part_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_60_part_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_90_part_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_150_part_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_225_part_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_300_part_comm)$coefficients[,4][2]), digits = 4))
  
  part_lm_df_comm <- data.frame(removed_factors, part_lm_adj_r2_comm, part_lm_pval_comm) %>%
    rename(Removed = removed_factors, Adj_R2 = part_lm_adj_r2_comm, pval = part_lm_pval_comm) %>%
    mutate(fraction = "Particle", test = "Community-Wide Production")
  
  
  ########## PER CAPITA PRODUCTION
  lm_01_part_percap <- lm(log10(fracprod_per_cell_noinf) ~ mean, data = dplyr::filter(dataframe,Removed == "1-tons" & fraction == "Particle"))
  lm_05_part_percap <- lm(log10(fracprod_per_cell_noinf)  ~ mean, data = dplyr::filter(dataframe, Removed == "5-tons" & fraction == "Particle"))
  lm_10_part_percap <- lm(log10(fracprod_per_cell_noinf)  ~ mean, data = dplyr::filter(dataframe, Removed == "10-tons" & fraction == "Particle"))
  lm_20_part_percap <- lm(log10(fracprod_per_cell_noinf)  ~ mean, data = dplyr::filter(dataframe, Removed == "20-tons" & fraction == "Particle"))
  lm_30_part_percap <- lm(log10(fracprod_per_cell_noinf)  ~ mean, data = dplyr::filter(dataframe, Removed == "30-tons" & fraction == "Particle"))
  lm_60_part_percap <- lm(log10(fracprod_per_cell_noinf)  ~ mean, data = dplyr::filter(dataframe, Removed == "60-tons" & fraction == "Particle"))
  lm_90_part_percap <- lm(log10(fracprod_per_cell_noinf)  ~ mean, data = dplyr::filter(dataframe, Removed == "90-tons" & fraction == "Particle"))
  lm_150_part_percap <- lm(log10(fracprod_per_cell_noinf)  ~ mean, data = dplyr::filter(dataframe, Removed == "150-tons" & fraction == "Particle"))
  lm_225_part_percap <- lm(log10(fracprod_per_cell_noinf)  ~ mean, data = dplyr::filter(dataframe, Removed == "225-tons" & fraction == "Particle"))
  lm_300_part_percap <- lm(log10(fracprod_per_cell_noinf)  ~ mean, data = dplyr::filter(dataframe, Removed == "300-tons" & fraction == "Particle"))
  
  part_lm_adj_r2_percap <- c(
    round(summary(lm_01_part_percap)$adj.r.squared, digits = 2),
    round(summary(lm_05_part_percap)$adj.r.squared, digits = 2),
    round(summary(lm_10_part_percap)$adj.r.squared, digits = 2),
    round(summary(lm_20_part_percap)$adj.r.squared, digits = 2),
    round(summary(lm_30_part_percap)$adj.r.squared, digits = 2),
    round(summary(lm_60_part_percap)$adj.r.squared, digits = 2),
    round(summary(lm_90_part_percap)$adj.r.squared, digits = 2),
    round(summary(lm_150_part_percap)$adj.r.squared, digits = 2),
    round(summary(lm_225_part_percap)$adj.r.squared, digits = 2),
    round(summary(lm_300_part_percap)$adj.r.squared, digits = 2))
  
  part_lm_pval_percap <- c(
    round(unname(summary(lm_01_part_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_05_part_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_10_part_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_20_part_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_30_part_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_60_part_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_90_part_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_150_part_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_225_part_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_300_part_percap)$coefficients[,4][2]), digits = 4))
  
  part_lm_df_percap <- data.frame(removed_factors, part_lm_adj_r2_percap, part_lm_pval_percap) %>%
    rename(Removed = removed_factors, Adj_R2 = part_lm_adj_r2_percap, pval = part_lm_pval_percap) %>%
    mutate(fraction = "Particle", test = "Per-Capita Production")
  
  
  ########### FREE LIVING ANALYSIS
  ##### Community-wide production
  lm_01_free_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe,  Removed == "1-tons" & fraction == "Free"))
  lm_05_free_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "5-tons" & fraction == "Free"))
  lm_10_free_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "10-tons" & fraction == "Free"))
  lm_20_free_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "20-tons" & fraction == "Free"))
  lm_30_free_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe,  Removed == "30-tons" & fraction == "Free"))
  lm_60_free_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "60-tons" & fraction == "Free"))
  lm_90_free_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "90-tons" & fraction == "Free"))
  lm_150_free_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe,  Removed == "150-tons" & fraction == "Free"))
  lm_225_free_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "225-tons" & fraction == "Free"))
  lm_300_free_comm <- lm(frac_bacprod ~ mean, data = dplyr::filter(dataframe, Removed == "300-tons" & fraction == "Free"))
  
  free_lm_adj_r2_comm <- c(
    round(summary(lm_01_free_comm)$adj.r.squared, digits = 2),
    round(summary(lm_05_free_comm)$adj.r.squared, digits = 2),
    round(summary(lm_10_free_comm)$adj.r.squared, digits = 2),
    round(summary(lm_20_free_comm)$adj.r.squared, digits = 2),
    round(summary(lm_30_free_comm)$adj.r.squared, digits = 2),
    round(summary(lm_60_free_comm)$adj.r.squared, digits = 2),
    round(summary(lm_90_free_comm)$adj.r.squared, digits = 2),
    round(summary(lm_150_free_comm)$adj.r.squared, digits = 2),
    round(summary(lm_225_free_comm)$adj.r.squared, digits = 2),
    round(summary(lm_300_free_comm)$adj.r.squared, digits = 2))
  
  free_lm_pval_comm <- c(
    round(unname(summary(lm_01_free_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_05_free_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_10_free_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_20_free_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_30_free_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_60_free_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_90_free_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_150_free_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_225_free_comm)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_300_free_comm)$coefficients[,4][2]), digits = 4))
  
  free_lm_df_comm <- data.frame(removed_factors, free_lm_adj_r2_comm, free_lm_pval_comm) %>%
    rename(Removed = removed_factors, Adj_R2 = free_lm_adj_r2_comm, pval = free_lm_pval_comm) %>%
    mutate(fraction = "Free", test = "Community-Wide Production")
  
  
  ##### Per-capita 
  lm_01_free_percap <- lm(log10(fracprod_per_cell_noinf) ~ mean, data = dplyr::filter(dataframe,  Removed == "1-tons" & fraction == "Free"))
  lm_05_free_percap <- lm(log10(fracprod_per_cell_noinf) ~ mean, data = dplyr::filter(dataframe, Removed == "5-tons" & fraction == "Free"))
  lm_10_free_percap <- lm(log10(fracprod_per_cell_noinf) ~ mean, data = dplyr::filter(dataframe, Removed == "10-tons" & fraction == "Free"))
  lm_20_free_percap <- lm(log10(fracprod_per_cell_noinf) ~ mean, data = dplyr::filter(dataframe, Removed == "20-tons" & fraction == "Free"))
  lm_30_free_percap <- lm(log10(fracprod_per_cell_noinf) ~ mean, data = dplyr::filter(dataframe,  Removed == "30-tons" & fraction == "Free"))
  lm_60_free_percap <- lm(log10(fracprod_per_cell_noinf) ~ mean, data = dplyr::filter(dataframe, Removed == "60-tons" & fraction == "Free"))
  lm_90_free_percap <- lm(log10(fracprod_per_cell_noinf) ~ mean, data = dplyr::filter(dataframe, Removed == "90-tons" & fraction == "Free"))
  lm_150_free_percap <- lm(log10(fracprod_per_cell_noinf) ~ mean, data = dplyr::filter(dataframe,  Removed == "150-tons" & fraction == "Free"))
  lm_225_free_percap <- lm(log10(fracprod_per_cell_noinf) ~ mean, data = dplyr::filter(dataframe, Removed == "225-tons" & fraction == "Free"))
  lm_300_free_percap <- lm(log10(fracprod_per_cell_noinf) ~ mean, data = dplyr::filter(dataframe, Removed == "300-tons" & fraction == "Free"))
  
  free_lm_adj_r2_percap <- c(
    round(summary(lm_01_free_percap)$adj.r.squared, digits = 2),
    round(summary(lm_05_free_percap)$adj.r.squared, digits = 2),
    round(summary(lm_10_free_percap)$adj.r.squared, digits = 2),
    round(summary(lm_20_free_percap)$adj.r.squared, digits = 2),
    round(summary(lm_30_free_percap)$adj.r.squared, digits = 2),
    round(summary(lm_60_free_percap)$adj.r.squared, digits = 2),
    round(summary(lm_90_free_percap)$adj.r.squared, digits = 2),
    round(summary(lm_150_free_percap)$adj.r.squared, digits = 2),
    round(summary(lm_225_free_percap)$adj.r.squared, digits = 2),
    round(summary(lm_300_free_percap)$adj.r.squared, digits = 2))
  
  free_lm_pval_percap <- c(
    round(unname(summary(lm_01_free_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_05_free_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_10_free_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_20_free_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_30_free_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_60_free_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_90_free_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_150_free_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_225_free_percap)$coefficients[,4][2]), digits = 4),
    round(unname(summary(lm_300_free_percap)$coefficients[,4][2]), digits = 4))
  
  free_lm_df_percap <- data.frame(removed_factors, free_lm_adj_r2_percap, free_lm_pval_percap) %>%
    rename(Removed = removed_factors, Adj_R2 = free_lm_adj_r2_percap, pval = free_lm_pval_percap) %>%
    mutate(fraction = "Free", test = "Per-Capita Production")  
  
  return(list(free_comm = free_lm_df_comm, free_percap = free_lm_df_percap, 
              part_comm = part_lm_df_comm, part_percap = part_lm_df_percap))
}

 

