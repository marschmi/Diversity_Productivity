# This file is where I am putting code from old analyses:
# Created: January 24th, 2016


# Vegan Alpha Diversity Analysis
## Remove samples with too few reads
#```{r rarefy-and-scale, eval = FALSE}
##  OTU: SUBSAMPLE AT 6600
# Prune samples out that have too few reads
otu_merged_pruned <- prune_samples(sample_sums(otu_merged) > 6600, otu_merged)
otu_merged_pruned <- prune_taxa(taxa_sums(otu_merged_pruned) > 0, otu_merged_pruned)
min(sample_sums(otu_merged_pruned))

# Scale the samples
scaled_otu_merged_pruned <- otu_merged_pruned %>%
  scale_reads(round = "round") 
#```



## Rarefy Read Depth Analysis ran on December 16th, 2016
#```{r calculate-alphadiv, eval = FALSE}
# # Initialize matrices to store richness and evenness estimates
# nsamp <- nsamples(otu_merged_pruned)
# min_lib <- min(sample_sums(otu_merged_pruned)) - 1
# 
# otu_richness <- matrix(nrow = nsamp, ncol = 100)
# row.names(otu_richness) <- sample_names(otu_merged_pruned)
# 
# otu_evenness <- matrix(nrow = nsamp, ncol = 100)
# row.names(otu_evenness) <- sample_names(otu_merged_pruned)
# 
# otu_shannon <- matrix(nrow = nsamp, ncol = 100)
# row.names(otu_shannon) <- sample_names(otu_merged_pruned)
# 
# # It is always important to set a seed when you subsample so your result is replicable
# set.seed(777)
# 
# for (i in 1:100) {
#  # Subsample
#  r <- rarefy_even_depth(otu_merged_pruned, sample.size = min_lib, verbose = FALSE, replace = TRUE)
# 
#  # Calculate richness
#  rich <- as.numeric(as.matrix(estimate_richness(r, measures = "Observed")))
#  otu_richness[ ,i] <- rich
# 
#  # Calculate evenness
#  even <- as.numeric(as.matrix(estimate_richness(r, measures = "InvSimpson")))
#  otu_evenness[ ,i] <- even
# 
#  # Calculate evenness
#  shannon <- as.numeric(as.matrix(estimate_richness(r, measures = "Shannon")))
#  otu_shannon[ ,i] <- shannon
# 
# }

#  write.table(otu_evenness, "data/metadata/otu_evenness100_rarefy6664", row.names = TRUE)
#  write.table(otu_richness, "data/metadata/otu_richness100_rarefy6664", row.names = TRUE)
#  write.table(otu_shannon, "data/metadata/otu_shannon100_rarefy6664", row.names = TRUE)
```

## Vegan Alpha Diversity Analysis
### OTU Analysis 
#```{r vegan-otu-alphadiv, eval = FALSE}
# Load values
nsamp <- nsamples(otu_merged_pruned)
min_lib <- min(sample_sums(otu_merged_pruned)) - 1

# Read in the files 
otu_richness <- read.table("metadata/otu_richness100_rarefy6664",  header = TRUE)
otu_evenness <- read.table("metadata/otu_evenness100_rarefy6664", header = TRUE)
otu_shannon <- read.table("metadata/otu_shannon100_rarefy6664", header = TRUE)

# Create a new dataframe to hold the means and standard deviations of richness estimates
norep_filter_name <- row.names(otu_richness)
mean <- apply(otu_richness, 1, mean)
sd <- apply(otu_richness, 1, sd)
measure <- rep("Richness", nsamp)
otu_rich_stats <- data.frame(norep_filter_name, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
norep_filter_name <- row.names(otu_evenness)
mean <- apply(otu_evenness, 1, mean)
sd <- apply(otu_evenness, 1, sd)
measure <- rep("Inverse_Simpson", nsamp)
otu_even_stats <- data.frame(norep_filter_name, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of shannon entropy estimates
norep_filter_name <- row.names(otu_shannon)
mean <- apply(otu_shannon, 1, mean)
sd <- apply(otu_shannon, 1, sd)
measure <- rep("Shannon_Entropy", nsamp)
otu_shan_stats <- data.frame(norep_filter_name, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of simpsons evenness estimates
norep_filter_name <- row.names(otu_evenness)
mean <- apply(otu_evenness, 1, mean)
sd <- apply(otu_evenness, 1, sd)
measure <- rep("Inverse_Simpson", nsamp)
otu_simpseven_stats <- data.frame(norep_filter_name, mean, sd, measure)

# Calculate Simpson's Evenness into new df called "simps_evenness"
otu_simps_evenness <- inner_join(otu_rich_stats, otu_even_stats, by = "norep_filter_name") %>%
  mutate(mean = mean.y/mean.x,
         sd = sd(mean),
         measure = "Simpsons_Evenness") %>%
  select(norep_filter_name, mean, sd, measure)

# Combine alpha diversity into one dataframe 
otu_alpha <- rbind(otu_rich_stats, otu_even_stats, otu_simps_evenness, otu_shan_stats)
s <- data.frame(sample_data(otu_merged_pruned))
otu_alphadiv <- merge(otu_alpha, s, by = "norep_filter_name") %>%
  filter(project == "Muskegon_Lake")

ggplot(filter(otu_alphadiv, measure == "Richness"), aes(x = norep_filter_name, y = mean, color = lakesite)) +
  geom_point() + facet_grid(project~fraction, scales = "free") + 
  xlab("Sample Name") + ylab("Richness") +
  theme(axis.text.x = element_text(angle = 30))  #Set the x-axis labels)

ggplot(filter(otu_alphadiv, measure == "Inverse_Simpson"), aes(x = norep_filter_name, y = mean, color = lakesite)) +
  geom_point() + facet_grid(project~fraction, scales = "free") + 
  xlab("Sample Name") + ylab("Inverse Simpson") +
  theme(axis.text.x = element_text(angle = 30))  #Set the x-axis labels)

ggplot(filter(otu_alphadiv, measure == "Shannon_Entropy"), aes(x = norep_filter_name, y = mean, color = lakesite)) +
  geom_point() + facet_grid(project~fraction, scales = "free") + 
  xlab("Sample Name") + ylab("Shannon Entropy") +
  theme(axis.text.x = element_text(angle = 30))  #Set the x-axis labels)


###############################
ML_otu_rich_stats <- filter(otu_alphadiv, measure == "Richness" & project == "Muskegon_Lake" & fraction %in% c("WholePart", "Free") & year == "2015")

# Can fractional production be predicted by richness? 
free_ML_otu_rich_stats <- filter(ML_otu_rich_stats, fraction == "Free")
freeprod_ML_otu_rich <- lm(frac_bacprod ~ mean, data = free_ML_otu_rich_stats)
freeprod_ML_otu_rich
summary(freeprod_ML_otu_rich)

part_ML_abs_rich_stats <- filter(ML_otu_rich_stats, fraction == "WholePart")
partprod_MLotu_rich <- lm(frac_bacprod ~ mean, data = part_ML_abs_rich_stats)
partprod_MLotu_rich
summary(partprod_MLotu_rich)

# Plot 
otu_rich_vegan <- ggplot(ML_otu_rich_stats, aes(x=mean, y=frac_bacprod, color = fraction)) + 
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), width = 0.2) + 
  ggtitle("OTU: Vegan") +
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "Free")) +
  scale_x_continuous(expand = c(0,0), limits = c(0,1250)) + 
  scale_y_continuous(limits = c(0,70),expand = c(0,0)) + 
  ylab("Production (μgC/L/hr)") + xlab("Observed Richness (D0)") +
  #geom_smooth(data=subset(ML_otu_rich_stats, fraction == "Free"), method='lm') + 
  geom_smooth(data=subset(ML_otu_rich_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.15,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 800, y=35, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_ML_otu_rich)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(freeprod_ML_otu_rich)$coefficients[,4][2]), digits = 4))) + 
  annotate("text", x = 900, y=5, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_MLotu_rich)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(partprod_MLotu_rich)$coefficients[,4][2]), digits = 4)));



######################################################### SHANNON ENTROPY
#  Subset a dataframe with the key information
ML_otu_shannon_stats <- filter(otu_alphadiv, 
                               measure == "Shannon_Entropy" & 
                                 project == "Muskegon_Lake" & 
                                 fraction %in% c("WholePart", "Free") & 
                                 year == "2015")

# Can fractional production be predicted by simpsevenness? 
free_ML_otu_shannon_stats <- filter(ML_otu_shannon_stats, fraction == "Free")
freeprod_ML_otu_shannon <- lm(frac_bacprod ~ mean, data = free_ML_otu_shannon_stats)
freeprod_ML_otu_shannon
summary(freeprod_ML_otu_shannon)

part_ML_abs_shannon_stats <- filter(ML_otu_shannon_stats, fraction == "WholePart")
partprod_MLotu_shannon <- lm(frac_bacprod ~ mean, data = part_ML_abs_shannon_stats)
partprod_MLotu_shannon
summary(partprod_MLotu_shannon)

# Plot 
otu_shannon_vegan <- ggplot(ML_otu_shannon_stats, aes(x=mean, y=frac_bacprod, color = fraction)) + 
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), width = 0.2) + 
  ggtitle("OTU: Vegan") +
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "Free")) +
  scale_x_continuous(expand = c(0,0), limits=c(3, 6)) + 
  scale_y_continuous(limits = c(0,70),expand = c(0,0)) + 
  ylab("Production (μgC/L/hr)") + xlab("Shannon Entropy (D1)") +
  #geom_smooth(data=subset(ML_otu_shannon_stats, fraction == "Free"), method='lm') + 
  geom_smooth(data=subset(ML_otu_shannon_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.15,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 3.5, y=35, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_ML_otu_shannon)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(freeprod_ML_otu_shannon)$coefficients[,4][2]), digits = 4))) + 
  annotate("text", x = 5.35, y=5, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_MLotu_shannon)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(partprod_MLotu_shannon)$coefficients[,4][2]), digits = 4))); 





######################################################### INVERSE SIMPSON
#  Subset a dataframe with the key information
ML_otu_invsimps_stats <- filter(otu_alphadiv, 
                                measure == "Inverse_Simpson" & 
                                  project == "Muskegon_Lake" & 
                                  fraction %in% c("WholePart", "Free") & 
                                  year == "2015")

# Can fractional production be predicted by invsimpsness? 
free_ML_otu_invsimps_stats <- filter(ML_otu_invsimps_stats, fraction == "Free")
freeprod_ML_otu_invsimps <- lm(frac_bacprod ~ mean, data = free_ML_otu_invsimps_stats)
freeprod_ML_otu_invsimps
summary(freeprod_ML_otu_invsimps)

part_ML_abs_invsimps_stats <- filter(ML_otu_invsimps_stats, fraction == "WholePart")
partprod_MLotu_invsimps <- lm(frac_bacprod ~ mean, data = part_ML_abs_invsimps_stats)
partprod_MLotu_invsimps
summary(partprod_MLotu_invsimps)

# Plot Simpson's Evenness
otu_invsimps_vegan <- ggplot(ML_otu_invsimps_stats, aes(x=mean, y=frac_bacprod, color = fraction)) + 
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), width = 0.2) + 
  ggtitle("OTU: Vegan") +
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "Free")) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,70),expand = c(0,0)) + 
  ylab("Production (μgC/L/hr)") + xlab("Inverse Simpson (D2)") +
  #geom_smooth(data=subset(ML_otu_invsimps_stats, fraction == "Free"), method='lm') + 
  geom_smooth(data=subset(ML_otu_invsimps_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.85,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 58, y=35, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_ML_otu_invsimps)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(freeprod_ML_otu_invsimps)$coefficients[,4][2]), digits = 4))) + 
  annotate("text", x = 63, y=5, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_MLotu_invsimps)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(partprod_MLotu_invsimps)$coefficients[,4][2]), digits = 4))); 




######################################################### SIMPSON'S EVENNESS
#  Subset a dataframe with the key information
ML_otu_simpseven_stats <- filter(otu_alphadiv, 
                                 measure == "Simpsons_Evenness" & 
                                   project == "Muskegon_Lake" & 
                                   fraction %in% c("WholePart", "Free") & 
                                   year == "2015")

# Can fractional production be predicted by simpsevenness? 
free_ML_otu_simpseven_stats <- filter(ML_otu_simpseven_stats, fraction == "Free")
freeprod_ML_otu_simpseven <- lm(frac_bacprod ~ mean, data = free_ML_otu_simpseven_stats)
freeprod_ML_otu_simpseven
summary(freeprod_ML_otu_simpseven)

part_ML_abs_simpseven_stats <- filter(ML_otu_simpseven_stats, fraction == "WholePart")
partprod_MLotu_simpseven <- lm(frac_bacprod ~ mean, data = part_ML_abs_simpseven_stats)
partprod_MLotu_simpseven
summary(partprod_MLotu_simpseven)

# Plot 
otu_simpseven_vegan <- ggplot(ML_otu_simpseven_stats, aes(x=mean, y=frac_bacprod, color = fraction)) + 
  geom_point(size = 3.5) +
  ggtitle("OTU: Vegan") +
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "Free")) +
  scale_x_continuous(expand = c(0,0), limits=c(0, 0.09)) + 
  scale_y_continuous(limits = c(0,70),expand = c(0,0)) + 
  ylab("Production (μgC/L/hr)") + xlab("Simpson's Evenness") +
  #geom_smooth(data=subset(ML_otu_simpseven_stats, fraction == "Free"), method='lm') + 
  geom_smooth(data=subset(ML_otu_simpseven_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.15,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 0.03, y=35, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_ML_otu_simpseven)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(freeprod_ML_otu_simpseven)$coefficients[,4][2]), digits = 4))) + 
  annotate("text", x = 0.015, y=7, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_MLotu_simpseven)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(partprod_MLotu_simpseven)$coefficients[,4][2]), digits = 4))); 

otu_vegan <- plot_grid(otu_rich_vegan, otu_simpseven_vegan, otu_shannon_vegan, otu_invsimps_vegan,
                       labels = c("A", "B", "C", "D"), 
                       align = "h", nrow = 2, ncol = 2)
otu_vegan

#ggsave("Figures/vegan_otu_alpha_vs_prod.png", otu_vegan, dpi = 600, units = "in", width = 10, height = 8)
#```
