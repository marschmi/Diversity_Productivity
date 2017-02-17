# Thursday, Feb. 16th, 2017


library(ggplot2)
library(devtools)
library(phyloseq)
library(tidyr)
library(dplyr)
library(cowplot)
library(picante) # Will also include ape and vegan 
source("Muskegon_functions.R")
source("set_colors.R")


# Load in community relative abundance matrix 
comm_RAREFIED_rm10 <- read.csv("comm_RAREFIED_rm10.csv", header = TRUE)
row.names(comm_RAREFIED_rm10) <- comm_RAREFIED_rm10$X
comm_RAREFIED_rm10$X = NULL
comm_RAREFIED_rm10 <- as.matrix(comm_RAREFIED_rm10)

# Load in Tree
phy_RAREFIED_rm10 <- read.nexus("phy_RAREFIED_rm10.csv")
# Calculate the phylogenetic distances
phy_dist_RAREFIED_rm10 <- cophenetic(phy_RAREFIED_rm10)


####################################################################################################
######################################## TAXA.LABELS ###############################################
####################################################################################################

######################  MEAN PAIRWISE DISTANCE 
# calculate standardized effect size mean pairwise distance (ses.mpd)
unweighted_sesmpd_taxalab <- ses.mpd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "taxa.labels", 
                                     abundance.weighted = FALSE, runs = 999)
# Write the file
write.csv(unweighted_sesmpd_taxalab, file = "mpd_mntd/unweighted_sesmpd_taxalab.csv", row.names = TRUE)


# CALCULATE THE WEIGHTED  MPD
weighted_sesmpd_taxalab <- ses.mpd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "taxa.labels", 
                              abundance.weighted = TRUE, runs = 999)
# Write the file 
write.csv(weighted_sesmpd_taxalab, file = "mpd_mntd/weighted_sesmpd_taxalab.csv", row.names = TRUE)



######################  MEAN NEAREST TAXON DISTANCE
# calculate UNWEIGHTED ses.mntd
unweighted_sesMNTD_taxalab <- ses.mntd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "taxa.labels", 
                                       abundance.weighted = FALSE, runs = 999)
# Write the file 
write.csv(unweighted_sesMNTD_taxalab, file = "mpd_mntd/unweighted_sesMNTD_taxalab.csv", row.names = TRUE)



# calculate WEIGHTED ses.mntd
weighted_sesMNTD_taxalab <- ses.mntd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "taxa.labels", 
                                   abundance.weighted = TRUE, runs = 999)
# Write the file 
write.csv(weighted_sesMNTD_taxalab, file = "mpd_mntd/weighted_sesMNTD_taxalab.csv", row.names = TRUE)


####################################################################################################
###################################### INDEPENDENT SWAP ############################################
####################################################################################################

# calculate standardized effect size mean pairwise distance (ses.mpd)
unweighted_sesMPD_indepswap <- ses.mpd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "independentswap", 
                                     abundance.weighted = FALSE, runs = 999)
# Write the file
write.csv(unweighted_sesMPD_indepswap, file = "mpd_mntd/unweighted_sesMPD_indepswap.csv", row.names = TRUE)


# CALCULATE THE WEIGHTED  MPD
weighted_sesMPD_indepswap <- ses.mpd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "independentswap", 
                                   abundance.weighted = TRUE, runs = 999)
# Write the file 
write.csv(weighted_sesMPD_indepswap, file = "mpd_mntd/weighted_sesMPD_indepswap.csv", row.names = TRUE)



######################  MEAN NEAREST TAXON DISTANCE
# calculate UNWEIGHTED ses.mntd
unweighted_sesMNTD_indepswap <- ses.mntd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "independentswap", 
                                       abundance.weighted = FALSE, runs = 999)
# Write the file 
write.csv(unweighted_sesMNTD_indepswap, file = "mpd_mntd/unweighted_sesMNTD_indepswap.csv", row.names = TRUE)


# calculate WEIGHTED ses.mntd
weighted_sesMNTD_indepswap <- ses.mntd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "independentswap", 
                                     abundance.weighted = TRUE, runs = 999)
# Write the file 
write.csv(weighted_sesMNTD_indepswap, file = "mpd_mntd/weighted_sesMNTD_indepswap.csv", row.names = TRUE)


####################################################################################################
########################################## RICHNESS ################################################
####################################################################################################

# calculate standardized effect size mean pairwise distance (ses.mpd)
unweighted_sesMPD_rich <- ses.mpd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "richness", 
                                       abundance.weighted = FALSE, runs = 999)
# Write the file
write.csv(unweighted_sesMPD_rich, file = "mpd_mntd/unweighted_sesMPD_rich.csv", row.names = TRUE)


# CALCULATE THE WEIGHTED  MPD
weighted_sesMPD_rich <- ses.mpd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "richness", 
                                     abundance.weighted = TRUE, runs = 999)
# Write the file 
write.csv(weighted_sesMPD_rich, file = "mpd_mntd/weighted_sesMPD_rich.csv", row.names = TRUE)



######################  MEAN NEAREST TAXON DISTANCE
# calculate UNWEIGHTED ses.mntd
unweighted_sesMNTD_rich <- ses.mntd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "richness", 
                                         abundance.weighted = FALSE, runs = 999)
# Write the file 
write.csv(unweighted_sesMNTD_rich, file = "mpd_mntd/unweighted_sesMNTD_rich.csv", row.names = TRUE)


# calculate WEIGHTED ses.mntd
weighted_sesMNTD_rich <- ses.mntd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "richness", 
                                       abundance.weighted = TRUE, runs = 999)
# Write the file 
write.csv(weighted_sesMNTD_rich, file = "mpd_mntd/weighted_sesMNTD_rich.csv", row.names = TRUE)



####################################################################################################
######################################### FREQUENCY ################################################
####################################################################################################

# calculate standardized effect size mean pairwise distance (ses.mpd)
unweighted_sesMPD_freq <- ses.mpd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "frequency", 
                                  abundance.weighted = FALSE, runs = 999)
# Write the file
write.csv(unweighted_sesMPD_freq, file = "mpd_mntd/unweighted_sesMPD_freq.csv", row.names = TRUE)


# CALCULATE THE WEIGHTED  MPD
weighted_sesMPD_freq <- ses.mpd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "frequency", 
                                abundance.weighted = TRUE, runs = 999)
# Write the file 
write.csv(weighted_sesMPD_freq, file = "mpd_mntd/weighted_sesMPD_freq.csv", row.names = TRUE)



######################  MEAN NEAREST TAXON DISTANCE
# calculate UNWEIGHTED ses.mntd
unweighted_sesMNTD_freq <- ses.mntd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "frequency", 
                                    abundance.weighted = FALSE, runs = 999)
# Write the file 
write.csv(unweighted_sesMNTD_freq, file = "mpd_mntd/unweighted_sesMNTD_freq.csv", row.names = TRUE)


# calculate WEIGHTED ses.mntd
weighted_sesMNTD_freq <- ses.mntd(comm_RAREFIED_rm10, phy_dist_RAREFIED_rm10, null.model = "frequency", 
                                  abundance.weighted = TRUE, runs = 999)
# Write the file 
write.csv(weighted_sesMNTD_freq, file = "mpd_mntd/weighted_sesMNTD_freq.csv", row.names = TRUE)
