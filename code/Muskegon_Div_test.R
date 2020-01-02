# Test for diversity differences between 2015 muskegon and rest of muskegon samples 
# Marian L Schmidt
# May 14, 2018



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


############################# PRUNE DATA for Muskegon Lake ############################### 
# Remove all samples except the ones used in this study! 
musk_lake_physeq <- subset_samples(otu_merged_musk_pruned, project == "Muskegon_Lake" & fraction != "Sediment")
musk_lake_physeq

# Remove taxa that are no longer in the phyloseq object
musk_lake_physeq_rm2 <- prune_taxa(taxa_sums(musk_lake_physeq) > 2, musk_lake_physeq) 
musk_lake_physeq_rm2


############################# REMOVE SAMPLES WITH FEW READS #############################
# What's the minimum sequencing size?
min(sample_sums(musk_lake_physeq_rm2))

# Calculate and plot the sample sequencing sums
sums_otu <- data.frame(rowSums(otu_table(musk_lake_physeq_rm2))) 
colnames(sums_otu) <- "Sample_TotalSeqs"
sums_otu$norep_filter_name <- row.names(sums_otu)
sums_otu <- arrange(sums_otu, Sample_TotalSeqs)

ggplot(filter(sums_otu, Sample_TotalSeqs < 6489), aes(x=reorder(norep_filter_name, Sample_TotalSeqs), y = Sample_TotalSeqs)) + 
  ylab("# of Sequences/Sample") +
  geom_bar(stat = "identity", colour="black",fill="slateblue1")  + xlab("Sample Name") + 
  ggtitle("OTU: Samples less than 10,000 reads") + 
  theme(axis.text.x = element_text(colour = "black", size=10, angle=45, hjust = 1, vjust = 1))

# Remove samples that are less than 6,488 (to match study)
musk_lake_physeq_rm2_pruned <- prune_samples(sample_sums(musk_lake_physeq_rm2) > 6488, musk_lake_physeq_rm2)
min(sample_sums(musk_lake_physeq_rm2_pruned))

############################# CALC DIVERSITY #############################


# Set the seed for reproducibility
set.seed(777)

# Calculate the alpha diversity with 100 repsampling events
alpha_calc <- calc_alpha_diversity(physeq = musk_lake_physeq_rm2_pruned)

# What was the minimum sample size? 
min(sample_sums(musk_lake_physeq_rm2_pruned)) - 1

# Put it altogether in a dataframe 
otu_alphadiv <- calc_mean_alphadiv(physeq = musk_lake_physeq_rm2_pruned,
                                   richness_df = alpha_calc$Richness, 
                                   evenness_df = alpha_calc$Inverse_Simpson, 
                                   shannon_df = alpha_calc$Shannon) %>%
  mutate(fraction = factor(fraction, levels = c("WholePart","Particle", "WholeFree","Free")),
         lakesite = factor(lakesite,  levels = c("Outlet", "Deep", "Bear", "River")),
         measure = factor(measure, levels = c("Richness",  "Shannon_Entropy", "Inverse_Simpson", "Simpsons_Evenness")),
         norep_water_name = paste(substr(norep_filter_name, 1, 4), substr(norep_filter_name, 6, 8), sep = ""))

otu_alphadiv <- otu_alphadiv %>%
  left_join(pca_scores_df, by = "norep_water_name")

# Make subsetted dataframes for each diversity metric
richness <- filter(otu_alphadiv, measure == "Richness")
shannon <- filter(otu_alphadiv, measure == "Shannon_Entropy")
invsimps <- filter(otu_alphadiv, measure == "Inverse_Simpson")
simpseven <- filter(otu_alphadiv, measure == "Simpsons_Evenness")

ggplot(richness, aes(y = mean, x = fraction, color = fraction, fill = fraction)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + geom_jitter() +
  scale_color_manual(values = fraction_colors) +
  scale_fill_manual(values = fraction_colors) 

rich_top <- filter(richness, limnion == "Top")
rich_top_PA <- filter(rich_top, fraction %in% c("WholePart", "Particle"))
rich_top_FL <- filter(rich_top, fraction %in% c("WholeFree", "Free"))
fit_PA <- aov(mean ~ fraction, data = rich_top_PA)
summary(fit_PA)
t.test(filter(rich_top_PA, fraction == "WholePart")$mean, filter(rich_top_PA, fraction == "Particle")$mean)
fit_FL <- aov(mean ~ fraction, data = rich_top_FL)
summary(fit_FL)
t.test(filter(rich_top_FL, fraction == "WholeFree")$mean, filter(rich_top_FL, fraction == "Free")$mean)


colorzz <- c(
  Particle = "firebrick3",
  Free = "cornflowerblue",
  WholePart =  "#FF6600", 
  WholeFree = "skyblue")

ggplot(rich_top, aes(y = mean, x = fraction, color = fraction, fill = fraction)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + geom_jitter(size = 2.5) +
  scale_color_manual(values = colorzz) + ylab("Observed Richness") +
  scale_fill_manual(values = colorzz) + xlab("Samples") +
  scale_y_continuous(limits = c(0, 1250), breaks = seq(0, 1250, by = 250), expand = c(0,0)) +
  annotate("text", x=1.5, y=1150, size = 4, color = "gray40", label= paste("NS \n |-- p = 0.709 --|")) +
  annotate("text", x=3.5, y=750, size = 4, color = "gray40", label= paste("NS \n|-- p = 0.281 --|")) +
  scale_x_discrete(labels = c('Particle-Associated \n This Study','Particle-Associated \n Previous Data',
                              'Free-Living \n This Study','Free-Living \n Previous Data')) +
  theme(legend.position = "none", legend.title = element_blank(),
        axis.title = element_text(face = "bold"))

ggsave("/Users/marschmi/git_repos/Diversity_Productivity/code/Div_comparison.png", dpi = 300)
