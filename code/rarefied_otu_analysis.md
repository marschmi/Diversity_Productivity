# Rarefied Diversity Analysis
Marian L. Schmidt  
February 14, 2017  
<style>
pre code, pre, code {
  white-space: pre !important;
  overflow-x: scroll !important;
  word-break: keep-all !important;
  word-wrap: initial !important;
}
</style>





```r
library(ggplot2)
library(devtools)
library(phyloseq)
library(tidyr)
library(dplyr)
library(cowplot)
library(picante) # Will also include ape and vegan 
library(car) # For residual analysis
source("Muskegon_functions.R")
source("set_colors.R")
```


# Rarefied Alpha Diversity Analysis
## Remove samples with too few reads

```r
# Loads a phyloseq object named otu_merged_musk_pruned)
load("../data/otu_merged_musk_pruned.RData")
# The name of the phyloseq object is: 
otu_merged_musk_pruned 
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 52980 taxa and 163 samples ]
## sample_data() Sample Data:       [ 163 samples by 70 sample variables ]
## tax_table()   Taxonomy Table:    [ 52980 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 52980 tips and 52978 internal nodes ]
```

```r
# Remove tree because it's too computationally intensive
otu_merged_musk_pruned <- merge_phyloseq(tax_table(otu_merged_musk_pruned), sample_data(otu_merged_musk_pruned), otu_table(otu_merged_musk_pruned))
otu_merged_musk_pruned
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 52980 taxa and 163 samples ]
## sample_data() Sample Data:       [ 163 samples by 70 sample variables ]
## tax_table()   Taxonomy Table:    [ 52980 taxa by 8 taxonomic ranks ]
```

```r
# Remove MOTHJ715 and MBRHP715
otu_merged_musk_pruned_noMOTHJ715_MBRHP715 <- subset_samples(otu_merged_musk_pruned, 
                                                             norep_filter_name != "MOTHJ715" & 
                                                               norep_filter_name != "MBRHP715")
otu_merged_musk_pruned_noMOTHJ715_MBRHP715
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 52980 taxa and 161 samples ]
## sample_data() Sample Data:       [ 161 samples by 70 sample variables ]
## tax_table()   Taxonomy Table:    [ 52980 taxa by 8 taxonomic ranks ]
```



```r
# Check the sequencing depth of each sample by 
sums_otu <- data.frame(rowSums(otu_table(otu_merged_musk_pruned_noMOTHJ715_MBRHP715))) 
colnames(sums_otu) <- "Sample_TotalSeqs"
sums_otu$names <- row.names(sums_otu)
sums_otu <- arrange(sums_otu, Sample_TotalSeqs) 
sums_otu <- make_metadata_norep(sums_otu)

## Add total sequences to metadata frame 
metdf <- sample_data(otu_merged_musk_pruned_noMOTHJ715_MBRHP715)
sums_otu$norep_filter_name <- sums_otu$names
sums_otu2 <- select(sums_otu, norep_filter_name, Sample_TotalSeqs)
metdf_num2 <- left_join(metdf, sums_otu2, by = "norep_filter_name") %>%
  dplyr::select(-one_of("D0", "D0_chao", "D1", "D2", "D0_SD", "D1_sd", "D0_chao_sd"))
row.names(metdf_num2) <- metdf_num2$norep_filter_name
# Rename the sample data 
sample_data(otu_merged_musk_pruned_noMOTHJ715_MBRHP715) <- metdf_num2
```

## Metadata Frames for Subsets

```r
### PREPARE DATA FRAMES FOR PHENOFLOW ALPHA DIVERSITY ANALYSIS
free_meta_data <- filter(metdf_num2, fraction %in% c("Free", "WholeFree") & norep_filter_name != "MOTHJ715" & limnion == "Top")
part_meta_data <- filter(metdf_num2, fraction %in% c("Particle", "WholePart") & norep_filter_name != "MOTHJ715" & limnion == "Top")
nosed_meta_data <- filter(metdf_num2, fraction != "Sediment" & limnion == "Top")

# Only "True Free Living and Particle" From 20um prefiltered samples 
free_only <- filter(free_meta_data, fraction == "Free")
part_only <- filter(part_meta_data, fraction == "Particle")

## 2015 specific samples that have NOT been prefiltered (whole water)
wholefree_only <- filter(free_meta_data, fraction == "WholeFree")
wholepart_only <- filter(part_meta_data, fraction == "WholePart")
```


#### Rarefy Read Depth Analysis ran on February 14th, 2016


# Fraction Diversity-Production Analysis 

```r
# Load values
nsamp <- nsamples(otu_merged_musk_pruned_noMOTHJ715_MBRHP715)
min_lib <- min(sample_sums(otu_merged_musk_pruned_noMOTHJ715_MBRHP715)) - 1
min_lib
```

```
## [1] 2895
```

```r
# Read in the files 
otu_richness <- read.table("../data/metadata/otu_richness100_rarefy2895",  header = TRUE)
otu_evenness <- read.table("../data/metadata/otu_evenness100_rarefy2895", header = TRUE)
otu_shannon <- read.table("../data/metadata/otu_shannon100_rarefy2895", header = TRUE)

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
s <- data.frame(sample_data(otu_merged_musk_pruned_noMOTHJ715_MBRHP715))
otu_alphadiv <- merge(otu_alpha, s, by = "norep_filter_name") %>%
  filter(project == "Muskegon_Lake" & limnion == "Top" & fraction != "Sediment") %>%
  mutate(fraction = factor(fraction, levels = c("WholePart", "Particle", "WholeFree", "Free")),
         lakesite = factor(lakesite,  levels = c("MOT", "MDP", "MBR", "MIN")),
         measure = factor(measure, levels = c("Richness", "Simpsons_Evenness", "Shannon_Entropy", "Inverse_Simpson"))) %>%
  # Calculate the per cell production rates 
  mutate(fraction_bac_abund = as.numeric(fraction_bac_abund),
         fracprod_per_cell = frac_bacprod/(1000*fraction_bac_abund),
         fracprod_per_cell_noinf = ifelse(fracprod_per_cell == Inf, NA, fracprod_per_cell))
```


## Subset Diversity Data 

```r
######################################################### RICHNESS
# Subset only richness data 
ML_otu_rich_stats <- filter(otu_alphadiv, measure == "Richness" & 
                              project == "Muskegon_Lake" & 
                              fraction %in% c("WholePart", "WholeFree") & year == "2015")

######################################################### SHANNON ENTROPY
# Subset only Shannon_Entropy data 
ML_otu_shannon_stats <- filter(otu_alphadiv, 
                               measure == "Shannon_Entropy" & 
                                 project == "Muskegon_Lake" & 
                                 fraction %in% c("WholePart", "WholeFree") & 
                                 year == "2015")

######################################################### INVERSE SIMPSON
# Subset only Inverse_Simpson data 
ML_otu_invsimps_stats <- filter(otu_alphadiv, 
                                measure == "Inverse_Simpson" & 
                                  project == "Muskegon_Lake" & 
                                  fraction %in% c("WholePart", "WholeFree") & 
                                  year == "2015")

######################################################### SIMPSON'S EVENNESS
# Subset only Simpsons_Evenness data 
ML_otu_simpseven_stats <- filter(otu_alphadiv, 
                                 measure == "Simpsons_Evenness" & 
                                   project == "Muskegon_Lake" & 
                                   fraction %in% c("WholePart", "WholeFree") & 
                                   year == "2015")
```



# Cell Count and Production Rates 

```r
######################################################### Fraction ABUNDANCe 
frac_abund_wilcox <- wilcox.test(log10(as.numeric(fraction_bac_abund)) ~ fraction, 
             data = filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & measure == "Richness"))
frac_abund_wilcox
```

```
## 
## 	Wilcoxon rank sum test
## 
## data:  log10(as.numeric(fraction_bac_abund)) by fraction
## W = 0, p-value = 1.479e-06
## alternative hypothesis: true location shift is not equal to 0
```

```r
filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & measure == "Richness") %>%
  group_by(fraction) %>%
  summarize(mean(as.numeric(fraction_bac_abund)))
```

```
## # A tibble: 2 × 2
##    fraction `mean(as.numeric(fraction_bac_abund))`
##      <fctr>                                  <dbl>
## 1 WholePart                               41168.88
## 2 WholeFree                              734522.25
```

```r
# Make a data frame to draw significance line in boxplot (visually calculated)
dat1 <- data.frame(a = c(1.15,1.15,1.85,1.85), b = c(6.45,6.5,6.5,6.45)) # WholePart vs WholeFree

poster_a <- ggplot(filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & 
                            measure == "Richness" & norep_filter_name != "MOTEJ515"), 
       aes(y = log10(as.numeric(fraction_bac_abund)), x = fraction)) +
  scale_color_manual(values = fraction_colors) + 
  scale_fill_manual(values = fraction_colors) +
  geom_jitter(size = 3, aes(color = fraction, fill = fraction), width = 0.2) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  ylab("Log10(Bacterial Cells/mL)") +
  scale_x_discrete(breaks=c("WholePart", "WholeFree"),
                      labels=c("Particle-\nAssociated", "Free-\nLiving")) + 
  #####  WHOLE PARTICLE VS WHOLE FREE CELL ABUNDANCES
  geom_path(data = dat1, aes(x = a, y = b), linetype = 1, color = "gray40") +
  annotate("text", x=1.5, y=6.5, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("***\np =", round(frac_abund_wilcox$p.value, digits = 6))) +
  theme(legend.position = "none",
        axis.title.x = element_blank())



######################################################### TOTAL PRODUCTION 
totprod_wilcox <- wilcox.test(frac_bacprod ~ fraction, 
             data = filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & measure == "Richness"))
totprod_wilcox
```

```
## 
## 	Wilcoxon rank sum test
## 
## data:  frac_bacprod by fraction
## W = 33, p-value = 0.02418
## alternative hypothesis: true location shift is not equal to 0
```

```r
filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & measure == "Richness") %>%
  group_by(fraction) %>%
  summarize(mean(frac_bacprod))
```

```
## # A tibble: 2 × 2
##    fraction `mean(frac_bacprod)`
##      <fctr>                <dbl>
## 1 WholePart              9.95449
## 2 WholeFree             24.07018
```

```r
# Make a data frame to draw significance line in boxplot (visually calculated)
dat2 <- data.frame(a = c(1.15,1.15,1.85,1.85), b = c(67,68,68,67)) # WholePart vs WholeFree


poster_b <- ggplot(filter(otu_alphadiv, 
                          fraction %in% c("WholePart", "WholeFree") & measure == "Richness"), 
       aes(y = frac_bacprod, x = fraction)) + 
  scale_color_manual(values = fraction_colors) + 
  scale_fill_manual(values = fraction_colors) +
  geom_jitter(size = 3, aes(color = fraction, fill = fraction), width = 0.2) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  ylab("Secondary Production (μgC/L/hr)") +
  scale_x_discrete(breaks=c("WholePart", "WholeFree"),
                    labels=c("Particle-\nAssociated", "Free-\nLiving")) + 
  #####  WHOLE PARTICLE VS WHOLE FREE TOTAL PRODUCTION 
  geom_path(data = dat2, aes(x = a, y = b), linetype = 1, color = "gray40") +
  annotate("text", x=1.5, y=68, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("***\np =", round(totprod_wilcox$p.value, digits = 3))) +
  theme(legend.position = "none",
        axis.title.x = element_blank())



######################################################### TOTAL PRODUCTION 
percellprod_wilcox <- wilcox.test(log10(fracprod_per_cell) ~ fraction, 
             data = filter(otu_alphadiv, 
                          fraction %in% c("WholePart", "WholeFree") & 
                            measure == "Richness" &
                            norep_filter_name != "MOTEJ515"))
percellprod_wilcox
```

```
## 
## 	Wilcoxon rank sum test
## 
## data:  log10(fracprod_per_cell) by fraction
## W = 125, p-value = 6.656e-05
## alternative hypothesis: true location shift is not equal to 0
```

```r
filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & measure == "Richness" &
       norep_filter_name != "MOTEJ515" & norep_filter_name != "MOTEP515") %>%
  group_by(fraction) %>%
  summarize(mean(fracprod_per_cell))
```

```
## # A tibble: 2 × 2
##    fraction `mean(fracprod_per_cell)`
##      <fctr>                     <dbl>
## 1 WholePart              4.818158e-07
## 2 WholeFree              3.869446e-08
```

```r
# Make a data frame to draw significance line in boxplot (visually calculated)
dat3 <- data.frame(a = c(1.15,1.15,1.85,1.85), b = c(-5.05,-5,-5,-5.05)) # WholePart vs WholeFree


poster_c <- ggplot(filter(otu_alphadiv, 
                          fraction %in% c("WholePart", "WholeFree") & 
                            measure == "Richness" &
                            norep_filter_name != "MOTEJ515" & norep_filter_name != "MOTEP515"), 
       aes(y = log10(fracprod_per_cell), x = fraction)) +
  scale_color_manual(values = fraction_colors) + 
  scale_fill_manual(values = fraction_colors) +
  geom_jitter(size = 3, aes(color = fraction, fill = fraction), width = 0.2) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  ylim(c(-8.5, -5)) + 
  ylab("log10(Production/cell) (μgC/cell/hr)") +
  scale_x_discrete(breaks=c("WholePart", "WholeFree"),
                    labels=c("Particle-\nAssociated", "Free-\nLiving")) + 
  #####  WHOLE PARTICLE VS WHOLE FREE PER CELL PRODUCTION 
  geom_path(data = dat3, aes(x = a, y = b), linetype = 1, color = "gray40") +
  annotate("text", x=1.5, y=-5, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("***\np =", round(percellprod_wilcox$p.value, digits = 5))) +
  theme(legend.position = "none",
        axis.title.x = element_blank())

plot_grid(poster_a, poster_b, poster_c,
          labels = c("A", "B", "C"),
          ncol = 3)
```

<img src="Rarefied_Figures/boxplot-cellcount-prod-1.png" style="display: block; margin: auto;" />



# Diversity Comparison

```r
######################################################### RICHNESS 
rich_wilcox <- wilcox.test(mean ~ fraction, 
             data = filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & measure == "Richness"))

filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & measure == "Richness") %>%
  group_by(fraction) %>%
  summarize(mean(mean),  median(mean))
```

```
## # A tibble: 2 × 3
##    fraction `mean(mean)` `median(mean)`
##      <fctr>        <dbl>          <dbl>
## 1 WholePart     454.6417        413.935
## 2 WholeFree     281.4133        257.595
```

```r
# Make a data frame to draw significance line in boxplot (visually calculated)
nums1 <- data.frame(a = c(1.15,1.15,1.85,1.85), b = c(790,800,800,790)) # WholePart vs WholeFree

poster_rich1 <-  ggplot(filter(otu_alphadiv, 
                          fraction %in% c("WholePart", "WholeFree") & measure == "Richness"), 
       aes(y = mean, x = fraction)) +
  ylab("Observed Richness") + 
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = fraction), alpha = 0.7) + # X-axis errorbars
  geom_jitter(size = 3, aes(color = fraction, fill = fraction), width = 0.2) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  scale_fill_manual(values = fraction_colors, 
                    breaks=c("WholeFree", "WholePart"), 
                    labels=c("Free-Living", "Particle-Associated")) +
  scale_color_manual(values = fraction_colors,
                 breaks=c("WholeFree", "WholePart"), 
                 labels=c("Free-Living", "Particle-Associated")) + 
  geom_path(data = nums1, aes(x = a, y = b), linetype = 1, color = "gray40") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(180, 810), breaks = c(200, 400, 600, 800)) 
     
     
poster_rich <- poster_rich1 + 
  annotate("text", x=1.5, y=800, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("***\np =", round(rich_wilcox$p.value, digits = 3))) +
  theme(legend.position = c(0.7, 0.75), 
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.spacing.y = unit(-0.5, "cm"),
        legend.spacing.x = unit(-0.3, "cm"),
        legend.box = "horizontal") +
  scale_fill_manual(values = fraction_colors, 
                    breaks=c("WholeFree", "WholePart"), 
                    labels=c("Free", "Particle")) +
  scale_color_manual(values = fraction_colors,
                 breaks=c("WholeFree", "WholePart"), 
                 labels=c("Free", "Particle")) +
  geom_jitter(size = 3, aes(color = fraction, fill = fraction, shape = lakesite), width = 0.2) + 
  scale_shape_discrete(guide = guide_legend(reverse=TRUE))



######################################################### SHANNON ENTROPY  
shannon_wilcox <- wilcox.test(mean ~ fraction, 
             data = filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & measure == "Shannon_Entropy"))

filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & measure == "Shannon_Entropy") %>%
  group_by(fraction) %>%
  summarize(mean(mean),  median(mean))
```

```
## # A tibble: 2 × 3
##    fraction `mean(mean)` `median(mean)`
##      <fctr>        <dbl>          <dbl>
## 1 WholePart     4.555905       4.464302
## 2 WholeFree     4.009610       4.001229
```

```r
# Make a data frame to draw significance line in boxplot (visually calculated)
nums2 <- data.frame(a = c(1.15,1.15,1.85,1.85), b = c(5.75,5.8,5.8,5.75)) # WholePart vs WholeFree

poster_shannon1 <- ggplot(filter(otu_alphadiv, 
                          fraction %in% c("WholePart", "WholeFree") & measure == "Shannon_Entropy"), 
       aes(y = mean, x = fraction)) +  
  ylab("Shannon Entropy") +   
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.7) + # X-axis errorbars
  geom_jitter(size = 3, aes(color = fraction, fill = fraction), width = 0.2) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  scale_fill_manual(values = fraction_colors, 
                    breaks=c("WholeFree", "WholePart"), 
                    labels=c("Free-Living", "Particle-Associated")) +
  scale_color_manual(values = fraction_colors,
                 breaks=c("WholeFree", "WholePart"), 
                 labels=c("Free-Living", "Particle-Associated")) + 
  geom_path(data = nums2, aes(x = a, y = b), linetype = 1, color = "gray40") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(3.4, 5.85), breaks = c(3.5, 4, 4.5, 5, 5.5))  

  
  
poster_shannon <- poster_shannon1 + 
  annotate("text", x=1.5, y=5.8, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("***\np =", round(shannon_wilcox$p.value, digits = 3))) +
  theme(legend.position = "none", axis.title.x = element_blank())


#########################################################  INVERSE SIMPSON 
simpson_wilcox <- wilcox.test(mean ~ fraction, 
             data = filter(otu_alphadiv, 
                          fraction %in% c("WholePart", "WholeFree") & 
                            measure == "Inverse_Simpson" &
                            norep_filter_name != "MOTEJ515" & norep_filter_name != "MOTEP515"))
simpson_wilcox
```

```
## 
## 	Wilcoxon rank sum test
## 
## data:  mean by fraction
## W = 70, p-value = 0.8328
## alternative hypothesis: true location shift is not equal to 0
```

```r
filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & measure == "Inverse_Simpson" &
       norep_filter_name != "MOTEJ515" & norep_filter_name != "MOTEP515") %>%
  group_by(fraction) %>%
  summarize(mean(mean), median(mean))
```

```
## # A tibble: 2 × 3
##    fraction `mean(mean)` `median(mean)`
##      <fctr>        <dbl>          <dbl>
## 1 WholePart     37.27907       19.45268
## 2 WholeFree     24.85521       23.88347
```

```r
nums3 <- data.frame(a = c(1.15,1.15,1.85,1.85), b = c(88,89, 89, 88)) # WholePart vs WholeFree


poster_invsimps1 <-   ggplot(filter(otu_alphadiv, 
                          fraction %in% c("WholePart", "WholeFree") & 
                            measure == "Inverse_Simpson" &
                            norep_filter_name != "MOTEJ515" & norep_filter_name != "MOTEP515"), 
       aes(y = mean, x = fraction)) +
  ylab("Inverse Simpson") +   
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.7) + # X-axis errorbars
  geom_jitter(size = 3, aes(color = fraction, fill = fraction), width = 0.2) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  scale_fill_manual(values = fraction_colors, 
                    breaks=c("WholeFree", "WholePart"), 
                    labels=c("Free-Living", "Particle-Associated")) +
  scale_color_manual(values = fraction_colors,
                 breaks=c("WholeFree", "WholePart"), 
                 labels=c("Free-Living", "Particle-Associated")) + 
  geom_path(data = nums3, aes(x = a, y = b), linetype = "dotted", color = "gray40") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0,95), breaks = c(0, 20, 40, 60, 80), expand = c(0,0))  


poster_invsimps <- poster_invsimps1 +
  annotate("text", x=1.5, y=89, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("NS\np =", round(simpson_wilcox$p.value, digits = 2))) +
  theme(legend.position = "none", axis.title.x = element_blank())

#########################################################  SIMPSON'S EVENNESS
simpseven_wilcox <- wilcox.test(mean ~ fraction, 
             data = filter(otu_alphadiv, 
                          fraction %in% c("WholePart", "WholeFree") & 
                            measure == "Simpsons_Evenness" &
                            norep_filter_name != "MOTEJ515" & norep_filter_name != "MOTEP515"))
simpseven_wilcox
```

```
## 
## 	Wilcoxon rank sum test
## 
## data:  mean by fraction
## W = 48, p-value = 0.2875
## alternative hypothesis: true location shift is not equal to 0
```

```r
filter(otu_alphadiv, fraction %in% c("WholePart", "WholeFree") & measure == "Simpsons_Evenness" &
       norep_filter_name != "MOTEJ515" & norep_filter_name != "MOTEP515") %>%
  group_by(fraction) %>%
  summarize(mean(mean), median(mean))
```

```
## # A tibble: 2 × 3
##    fraction `mean(mean)` `median(mean)`
##      <fctr>        <dbl>          <dbl>
## 1 WholePart   0.07295660     0.05540021
## 2 WholeFree   0.08859526     0.08346213
```

```r
nums3 <- data.frame(a = c(1.15,1.15,1.85,1.85), b = c(0.145, 0.15, 0.15, 0.145)) # WholePart vs WholeFree


poster_simpseven1 <-   ggplot(filter(otu_alphadiv, 
                          fraction %in% c("WholePart", "WholeFree") & 
                            measure == "Simpsons_Evenness" &
                            norep_filter_name != "MOTEJ515" & norep_filter_name != "MOTEP515"), 
       aes(y = mean, x = fraction)) +
  ylab("Simpson's Evenness") +   
  #geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.7) + # X-axis errorbars
  geom_jitter(size = 3, aes(color = fraction, fill = fraction), width = 0.2) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  scale_fill_manual(values = fraction_colors, 
                    breaks=c("WholeFree", "WholePart"), 
                    labels=c("Free-Living", "Particle-Associated")) +
  scale_color_manual(values = fraction_colors,
                 breaks=c("WholeFree", "WholePart"), 
                 labels=c("Free-Living", "Particle-Associated")) + 
  geom_path(data = nums3, aes(x = a, y = b), linetype = "dotted", color = "gray40") +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0.03,0.151), breaks = c(0.05, 0.1, 0.15))  

  
poster_simpseven <- poster_simpseven1 +
  annotate("text", x=1.5, y=0.15, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("NS\np =", round(simpseven_wilcox$p.value, digits = 2))) +
  theme(legend.position = "none", axis.title.x = element_blank())

plot_grid(poster_rich, poster_shannon, poster_invsimps, poster_simpseven,
          labels = c("A", "B", "C", "D"),
          ncol = 4)
```

<img src="Rarefied_Figures/boxplot-diversity-comparison-1.png" style="display: block; margin: auto;" />


# Diversity vs Fraction Production 

```r
######################################################### RICHNESS
# Free-Living Richness vs fractional production 
freeprod_ML_otu_rich <- lm(frac_bacprod ~ mean, data = filter(ML_otu_rich_stats, fraction == "WholeFree"))
summary(freeprod_ML_otu_rich) 
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean, data = filter(ML_otu_rich_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -21.741 -12.261  -1.552   8.540  29.066 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)  0.51721   18.32443   0.028    0.978
## mean         0.08370    0.06275   1.334    0.212
## 
## Residual standard error: 16.94 on 10 degrees of freedom
## Multiple R-squared:  0.151,	Adjusted R-squared:  0.06612 
## F-statistic: 1.779 on 1 and 10 DF,  p-value: 0.2119
```

```r
# Particle-Associated Richness vs fractional production 
partprod_MLotu_rich <- lm(frac_bacprod ~ mean, data = filter(ML_otu_rich_stats, fraction == "WholePart"))
summary(partprod_MLotu_rich) 
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean, data = filter(ML_otu_rich_stats, 
##     fraction == "WholePart"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -6.9920 -3.6940 -0.5982  3.2391 11.3706 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept) -9.04480    5.47348  -1.652  0.12944   
## mean         0.04179    0.01149   3.638  0.00455 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.673 on 10 degrees of freedom
## Multiple R-squared:  0.5696,	Adjusted R-squared:  0.5265 
## F-statistic: 13.23 on 1 and 10 DF,  p-value: 0.004553
```

```r
# Plot 
otu_rich_vegan <-  ggplot(ML_otu_rich_stats, aes(x=mean, y=frac_bacprod, color = fraction)) + 
  geom_point(size = 3.5) + 
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), alpha = 0.7) + # X-axis errorbars
  # Y-axis errorbars
  geom_errorbar(aes(ymin = frac_bacprod - SD_frac_bacprod, ymax = frac_bacprod + SD_frac_bacprod),  alpha = 0.5) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), 
                     limits = c("WholePart", "WholeFree"),
                     breaks=c("WholePart", "WholeFree"),
                     labels=c("Particle", "Free")) + 
  ylab("Secondary Production (μgC/L/hr)") + xlab("Observed Richness") +
  geom_smooth(data=subset(ML_otu_rich_stats, fraction == "WholePart"), method='lm') + 
  scale_x_continuous(limits = c(180, 810), breaks = c(200, 400, 600, 800)) + 
  theme(legend.position=c(0.15,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 250, y=55, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_ML_otu_rich)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(freeprod_ML_otu_rich)$coefficients[,4][2]), digits = 2))) + 
  annotate("text", x = 650, y=3, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_MLotu_rich)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(partprod_MLotu_rich)$coefficients[,4][2]), digits = 4)));


######################################################### SHANNON ENTROPY
# Free-Living Shannon vs fractional production 
freeprod_ML_otu_shannon <- lm(frac_bacprod ~ mean, data = filter(ML_otu_shannon_stats, fraction == "WholeFree"))
summary(freeprod_ML_otu_shannon)
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean, data = filter(ML_otu_shannon_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -22.642  -8.587  -3.915   7.135  34.659 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)   -31.96      67.18  -0.476    0.645
## mean           13.97      16.71   0.836    0.422
## 
## Residual standard error: 17.78 on 10 degrees of freedom
## Multiple R-squared:  0.06538,	Adjusted R-squared:  -0.02808 
## F-statistic: 0.6995 on 1 and 10 DF,  p-value: 0.4225
```

```r
# Particle-Associated Shannon vs fractional production 
partprod_MLotu_shannon <- lm(frac_bacprod ~ mean, data = filter(ML_otu_shannon_stats, fraction == "WholePart"))
summary(partprod_MLotu_shannon)
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean, data = filter(ML_otu_shannon_stats, 
##     fraction == "WholePart"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -8.0980 -3.2043 -0.3301  1.3486 12.1080 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)  -38.540     13.486  -2.858  0.01702 * 
## mean          10.644      2.938   3.623  0.00467 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.686 on 10 degrees of freedom
## Multiple R-squared:  0.5676,	Adjusted R-squared:  0.5243 
## F-statistic: 13.13 on 1 and 10 DF,  p-value: 0.004667
```

```r
# Plot 
otu_shannon_vegan <- ggplot(ML_otu_shannon_stats, aes(x=mean, y=frac_bacprod, color = fraction)) + 
  geom_point(size = 3.5) + 
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), alpha = 0.7) + # X-axis errorbars
  # Y-axis errorbars
  geom_errorbar(aes(ymin = frac_bacprod - SD_frac_bacprod, ymax = frac_bacprod + SD_frac_bacprod),  alpha = 0.5) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), 
                     limits = c("WholePart", "WholeFree"),
                     breaks=c("WholeFree", "WholePart"),
                     labels=c("Free-Living", "Particle-Associated")) + 
  ylab("Secondary Production (μgC/L/hr)") + xlab("Shannon Entropy") +
  scale_x_continuous(limits = c(3.4, 5.85), breaks = c(3.5, 4, 4.5, 5, 5.5)) + 
  geom_smooth(data=subset(ML_otu_shannon_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.15,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 3.75, y=55, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_ML_otu_shannon)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(freeprod_ML_otu_shannon)$coefficients[,4][2]), digits = 2))) + 
  annotate("text", x = 5.35, y=3, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_MLotu_shannon)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(partprod_MLotu_shannon)$coefficients[,4][2]), digits = 4))); 


######################################################### INVERSE SIMPSON
# Free-Living Inverse Simpson vs fractional production 
freeprod_ML_otu_invsimps <- lm(frac_bacprod ~ mean, data = filter(ML_otu_invsimps_stats, fraction == "WholeFree"))
summary(freeprod_ML_otu_invsimps)
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean, data = filter(ML_otu_invsimps_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -20.569 -10.734  -4.013   6.290  34.864 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)   9.8575    15.6009   0.632    0.542
## mean          0.5718     0.5935   0.963    0.358
## 
## Residual standard error: 17.59 on 10 degrees of freedom
## Multiple R-squared:  0.08494,	Adjusted R-squared:  -0.006562 
## F-statistic: 0.9283 on 1 and 10 DF,  p-value: 0.358
```

```r
# Particle-Associated Inverse Simpson vs fractional production 
partprod_MLotu_invsimps <- lm(frac_bacprod ~ mean, data = filter(ML_otu_invsimps_stats, fraction == "WholePart"))
summary(partprod_MLotu_invsimps)
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean, data = filter(ML_otu_invsimps_stats, 
##     fraction == "WholePart"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -7.4169 -2.1453 -0.2015  0.9740  7.8468 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -0.12549    2.37217  -0.053 0.958852    
## mean         0.26871    0.05264   5.105 0.000461 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 4.554 on 10 degrees of freedom
## Multiple R-squared:  0.7227,	Adjusted R-squared:  0.6949 
## F-statistic: 26.06 on 1 and 10 DF,  p-value: 0.0004608
```

```r
# Plot Simpson's Evenness
otu_invsimps_vegan <- ggplot(ML_otu_invsimps_stats, aes(x=mean, y=frac_bacprod, color = fraction)) + 
  geom_point(size = 3.5) +  
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), alpha = 0.7) + # X-axis errorbars
  # Y-axis errorbars
  geom_errorbar(aes(ymin = frac_bacprod - SD_frac_bacprod, ymax = frac_bacprod + SD_frac_bacprod),  alpha = 0.5) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), 
                     limits = c("WholePart", "WholeFree"),
                     breaks=c("WholeFree", "WholePart"),
                     labels=c("Free-Living", "Particle-Associated")) + 
  scale_x_continuous(limits = c(0,95), breaks = c(0, 20, 40, 60, 80), expand = c(0,0)) + 
  ylab("Secondary Production (μgC/L/hr)") + xlab("Inverse Simpson") +
  geom_smooth(data=subset(ML_otu_invsimps_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.85,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 15, y=55, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_ML_otu_invsimps)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(freeprod_ML_otu_invsimps)$coefficients[,4][2]), digits = 2))) + 
  annotate("text", x = 63, y=3, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_MLotu_invsimps)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(partprod_MLotu_invsimps)$coefficients[,4][2]), digits = 4))); 



######################################################### SIMPSON'S EVENNESS
# Free-Living Simpson's Evenness vs fractional production 
freeprod_ML_otu_simpseven <- lm(frac_bacprod ~ mean, data = filter(ML_otu_simpseven_stats, fraction == "WholeFree"))
summary(freeprod_ML_otu_simpseven)
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean, data = filter(ML_otu_simpseven_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -19.709 -12.801  -2.287   5.341  39.838 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)    20.41      24.17   0.845    0.418
## mean           41.29     266.14   0.155    0.880
## 
## Residual standard error: 18.37 on 10 degrees of freedom
## Multiple R-squared:  0.002401,	Adjusted R-squared:  -0.09736 
## F-statistic: 0.02406 on 1 and 10 DF,  p-value: 0.8798
```

```r
# Particle-Associated Simpson's Evenness vs fractional production 
partprod_MLotu_simpseven <- lm(frac_bacprod ~ mean, data = filter(ML_otu_simpseven_stats, fraction == "WholePart"))
summary(partprod_MLotu_simpseven)
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean, data = filter(ML_otu_simpseven_stats, 
##     fraction == "WholePart"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -10.3301  -1.9616  -0.8216   1.2507  10.5267 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)   -5.080      4.277  -1.188  0.26236   
## mean         199.911     52.742   3.790  0.00354 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.54 on 10 degrees of freedom
## Multiple R-squared:  0.5896,	Adjusted R-squared:  0.5486 
## F-statistic: 14.37 on 1 and 10 DF,  p-value: 0.003541
```

```r
# Plot 
otu_simpseven_vegan <- ggplot(ML_otu_simpseven_stats, aes(x=mean, y=frac_bacprod, color = fraction)) + 
  geom_point(size = 3.5) +  
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), alpha = 0.7) + # X-axis errorbars
  # Y-axis errorbars
  geom_errorbar(aes(ymin = frac_bacprod - SD_frac_bacprod, ymax = frac_bacprod + SD_frac_bacprod),  alpha = 0.5) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), 
                     limits = c("WholePart", "WholeFree"),
                     breaks=c("WholeFree", "WholePart"),
                     labels=c("Free-Living", "Particle-Associated")) + 
  scale_x_continuous(limits = c(0.03,0.151), breaks = c(0.05, 0.1, 0.15))  +
  ylab("Secondary Production (μgC/L/hr)") + xlab("Simpson's Evenness") +
  geom_smooth(data=subset(ML_otu_simpseven_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.15,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 0.05, y=55, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_ML_otu_simpseven)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(freeprod_ML_otu_simpseven)$coefficients[,4][2]), digits = 2))) + 
  annotate("text", x = 0.125, y=3, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_MLotu_simpseven)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(partprod_MLotu_simpseven)$coefficients[,4][2]), digits = 4))); 

#otu_vegan <- plot_grid(otu_rich_vegan, otu_simpseven_vegan, otu_shannon_vegan, otu_invsimps_vegan,
#                       labels = c("A", "B", "C", "D"), 
#                       align = "h", nrow = 2, ncol = 2)
#otu_vegan
```


## Residual Analysis 

```r
##########################################################################
#############################   RESIDUALS   ##############################
##########################################################################

######################################################### RICHNESS
# Residual analysis of the RICHNESS Models
plot_residuals(lm_model = partprod_MLotu_rich, 
               lm_observed_y = filter(ML_otu_rich_stats, fraction == "WholePart")$frac_bacprod)
```

<img src="Rarefied_Figures/check-lm-residuals-1.png" style="display: block; margin: auto;" />

```r
######################################################### SHANNON ENTROPY
# Residual analysis of the SHANNON ENTROPY Models
plot_residuals(lm_model = partprod_MLotu_shannon, 
               lm_observed_y = filter(ML_otu_shannon_stats, fraction == "WholePart")$frac_bacprod)
```

<img src="Rarefied_Figures/check-lm-residuals-2.png" style="display: block; margin: auto;" />

```r
######################################################### INVERSE SIMPSON
# Residual analysis of the INVERSE SIMPSON Models
plot_residuals(lm_model = partprod_MLotu_invsimps, 
               lm_observed_y = filter(ML_otu_invsimps_stats, fraction == "WholePart")$frac_bacprod)
```

<img src="Rarefied_Figures/check-lm-residuals-3.png" style="display: block; margin: auto;" />

```r
######################################################### SIMPSONS EVENNESS
# Residual analysis of the INVERSE SIMPSON Models
plot_residuals(lm_model = partprod_MLotu_simpseven, 
               lm_observed_y = filter(ML_otu_simpseven_stats, fraction == "WholePart")$frac_bacprod)
```

<img src="Rarefied_Figures/check-lm-residuals-4.png" style="display: block; margin: auto;" />


## Correlations 


```r
##########################################################################
###########################   CORRELATIONS   #############################
##########################################################################
# RICHNESS vs SHANNON
cor(filter(ML_otu_rich_stats, fraction == "WholePart")$mean,
    filter(ML_otu_shannon_stats, fraction == "WholePart")$mean) # YES
```

```
## [1] 0.961703
```

```r
# SHANNON VS INVERSE SIMPSON
cor(filter(ML_otu_shannon_stats, fraction == "WholePart")$mean,
    filter(ML_otu_invsimps_stats, fraction == "WholePart")$mean) # YES
```

```
## [1] 0.9682435
```

```r
# INVERSE SIMPSON VS SIMPSONS EVENNESS
cor(filter(ML_otu_invsimps_stats, fraction == "WholePart")$mean,
    filter(ML_otu_simpseven_stats, fraction == "WholePart")$mean) # YES
```

```
## [1] 0.9284642
```

```r
# SIMPSONS EVENNESS VS RICHNESS
cor(filter(ML_otu_simpseven_stats, fraction == "WholePart")$mean,
    filter(ML_otu_rich_stats, fraction == "WholePart")$mean) # YES
```

```
## [1] 0.7683792
```


## Post-hoc analysis 

```r
# Are the fractions different from each other in predicting fraction production?
lm_by_divmeasure <- lm(frac_bacprod ~ mean/measure, data=otu_alphadiv)
summary(lm_by_divmeasure)
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean/measure, data = otu_alphadiv)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -15.504 -11.818  -5.381   9.078  46.809 
## 
## Coefficients:
##                                Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                   14.786690   3.548610   4.167 4.71e-05 ***
## mean                           0.004952   0.010170   0.487    0.627    
## mean:measureSimpsons_Evenness 32.272480  45.724725   0.706    0.481    
## mean:measureShannon_Entropy    0.470810   0.940750   0.500    0.617    
## mean:measureInverse_Simpson    0.063010   0.099439   0.634    0.527    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 15.07 on 187 degrees of freedom
##   (88 observations deleted due to missingness)
## Multiple R-squared:  0.00297,	Adjusted R-squared:  -0.01836 
## F-statistic: 0.1393 on 4 and 187 DF,  p-value: 0.9675
```

```r
# Run a post-hoc test
library(multcomp)
post_hoc_measure <- glht(lm_by_divmeasure, linfct = mcp(measure = "Tukey", interaction_average=TRUE),
                vcov=vcovHC(lm_by_divmeasure, type = "HC0"))
```

```
## Error in glht.mcp(lm_by_divmeasure, linfct = mcp(measure = "Tukey", interaction_average = TRUE), : could not find function "vcovHC"
```

```r
summary(post_hoc_measure)
```

```
## Error in summary(post_hoc_measure): object 'post_hoc_measure' not found
```

```r
detach("package:multcomp", unload=TRUE) # This package masks the dplyr select function = :(
```






```r
plot_grid(poster_rich1 + xlab("Fraction \n") + ylab("Observed Richness") + 
                theme(legend.position = "none", axis.text.y = element_blank()) +coord_flip() +
                annotate("text", x=1.5, y=700, fontface = "bold",  size = 4, color = "gray40",
                          label= paste("***\np =", round(rich_wilcox$p.value, digits = 3))), 
          poster_shannon1 + xlab("Fraction  \n") + ylab("Shannon Entropy") + 
                theme(legend.position = "none", axis.text.y = element_blank()) +
                coord_flip() +
                annotate("text", x=1.5, y=5.45, fontface = "bold",  size = 4, color = "gray40",
                          label= paste("**\np =", round(shannon_wilcox$p.value, digits = 3))),  
          poster_invsimps1 + xlab("Fraction \n") + ylab("Inverse Simpson") + 
                theme(legend.position = c(0.78, 0.90), axis.text.y = element_blank(), legend.title = element_blank()) +
                coord_flip() +
                annotate("text", x=1.5, y=78, fontface = "bold",  size = 4, color = "gray40",
                          label= paste("NS\np =", round(simpson_wilcox$p.value, digits = 2))), 
          poster_simpseven1 + xlab("Fraction \n") + ylab("Simpson's Evenness") + 
                theme(legend.position = "none", axis.text.y = element_blank()) +
                coord_flip() +
                annotate("text", x=1.5, y=0.135, fontface = "bold",  size = 4, color = "gray40",
                          label= paste("NS\np =", round(simpseven_wilcox$p.value, digits = 2))), 
          otu_rich_vegan + theme(legend.position = "none"), 
          otu_shannon_vegan + theme(legend.position = "none"), 
          otu_invsimps_vegan + theme(legend.position = c(0.7, 0.80)), 
          otu_simpseven_vegan + theme(legend.position = "none"), 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 4, nrow = 2)
```

<img src="Rarefied_Figures/richness-div-boxplot-BEF-1.png" style="display: block; margin: auto;" />




```r
plot_grid(poster_rich1 + xlab("\n Fraction \n") + ylab("Observed Richness") + 
                theme(legend.position = "none", axis.text.y = element_blank()) +coord_flip() +
                annotate("text", x=1.5, y=700, fontface = "bold",  size = 4, color = "gray40",
                          label= paste("***\np =", round(rich_wilcox$p.value, digits = 3))), 
          poster_invsimps1 + xlab("\n Fraction \n") + ylab("Inverse Simpson") + 
                theme(legend.position = c(0.78, 0.90), 
                      axis.text.y = element_blank(), 
                      legend.title = element_blank()) +
                coord_flip() +
                annotate("text", x=1.5, y=78, fontface = "bold",  size = 4, color = "gray40",
                          label= paste("NS\np =", round(simpson_wilcox$p.value, digits = 2))), 
          otu_rich_vegan +  ylab("Heterotrophic Production \n(μgC/L/hr)") + 
                theme(legend.position = "none"), 
          otu_invsimps_vegan + ylab("Heterotrophic Production \n(μgC/L/hr)") +
                theme(legend.position = c(0.7, 0.80)), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
```

<img src="Rarefied_Figures/rich-and-invsimps-1.png" style="display: block; margin: auto;" />


```r
plot_grid(poster_shannon1 + xlab("\n Fraction  \n") + ylab("Shannon Entropy") + 
                theme(legend.position = c(0.78, 0.90), axis.text.y = element_blank()) +
                coord_flip() +
                annotate("text", x=1.5, y=5.45, fontface = "bold",  size = 4, color = "gray40",
                          label= paste("**\np =", round(shannon_wilcox$p.value, digits = 3))),  
          poster_simpseven1 + xlab("\n Fraction \n") + ylab("Simpson's Evenness") + 
                theme(legend.position = "none", axis.text.y = element_blank()) +
                coord_flip() +
                annotate("text", x=1.5, y=0.135, fontface = "bold",  size = 4, color = "gray40",
                          label= paste("NS\np =", round(simpseven_wilcox$p.value, digits = 2))), 
          otu_shannon_vegan + ylab("Heterotrophic Production \n(μgC/L/hr)") + 
                theme(legend.position = c(0.7, 0.75)), 
          otu_simpseven_vegan + ylab("Heterotrophic Production \n(μgC/L/hr)") + 
                theme(legend.position = "none"), 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
```

<img src="Rarefied_Figures/shannon-and-simpseven-1.png" style="display: block; margin: auto;" />



# Per cell Fraction Production vs Diversity

```r
#########################################################  RICHNESS 
# Free-Living Richness vs fractional production per cell 
freeprod_percell_ML_otu_rich <- lm(log10(fracprod_per_cell) ~ mean, data = filter(ML_otu_rich_stats, fraction == "WholeFree"))
summary(freeprod_percell_ML_otu_rich)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(ML_otu_rich_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.7056 -0.1326  0.0992  0.2334  0.4650 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -8.231139   0.392770 -20.957 1.36e-09 ***
## mean         0.002332   0.001345   1.734    0.114    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3632 on 10 degrees of freedom
## Multiple R-squared:  0.2311,	Adjusted R-squared:  0.1543 
## F-statistic: 3.006 on 1 and 10 DF,  p-value: 0.1136
```

```r
# Particle-Associated Richness vs fractional production per cell 
partprod_percell_MLotu_rich <- lm(log10(fracprod_per_cell) ~ mean, data = filter(filter(ML_otu_rich_stats, fraction == "WholePart" & fracprod_per_cell != Inf)))
summary(partprod_percell_MLotu_rich)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(filter(ML_otu_rich_stats, 
##     fraction == "WholePart" & fracprod_per_cell != Inf)))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.4778 -0.2188 -0.0434  0.1136  0.6824 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -7.9295357  0.3569740 -22.213 3.59e-09 ***
## mean         0.0026197  0.0007395   3.543  0.00629 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3628 on 9 degrees of freedom
## Multiple R-squared:  0.5824,	Adjusted R-squared:  0.536 
## F-statistic: 12.55 on 1 and 9 DF,  p-value: 0.006288
```

```r
# Plot 
rich_vs_fracprod_percell <- ggplot(filter(ML_otu_rich_stats, fracprod_per_cell != Inf),
       aes(x=mean, y=log10(fracprod_per_cell), color = fraction)) + 
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  scale_x_continuous(limits = c(180, 810), breaks = c(200, 400, 600, 800)) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") + 
  xlab("Observed Richness") +
  geom_smooth(data=subset(ML_otu_rich_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 500, y=-8, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_percell_ML_otu_rich)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(freeprod_percell_ML_otu_rich)$coefficients[,4][2]), digits = 2))) + 
  annotate("text", x = 650, y=-7, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_percell_MLotu_rich)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(partprod_percell_MLotu_rich)$coefficients[,4][2]), digits = 3)));


#########################################################  SHANNON
# Free-Living Shannon vs fractional production per cell 
freeprod_percell_ML_otu_shannon <- lm(log10(fracprod_per_cell) ~ mean, data = filter(ML_otu_shannon_stats, fraction == "WholeFree"))
summary(freeprod_percell_ML_otu_shannon)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(ML_otu_shannon_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.73908 -0.12432  0.06301  0.16396  0.60637 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  -9.3026     1.4662  -6.345 8.41e-05 ***
## mean          0.4309     0.3646   1.182    0.265    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.388 on 10 degrees of freedom
## Multiple R-squared:  0.1226,	Adjusted R-squared:  0.03482 
## F-statistic: 1.397 on 1 and 10 DF,  p-value: 0.2646
```

```r
# Particle-Associated Shannon vs fractional production per cell 
partprod_percell_MLotu_shannon <- lm(log10(fracprod_per_cell) ~ mean, 
                                     data = filter(filter(ML_otu_shannon_stats, fraction == "WholePart" & fracprod_per_cell != Inf)))
summary(partprod_percell_MLotu_shannon)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(filter(ML_otu_shannon_stats, 
##     fraction == "WholePart" & fracprod_per_cell != Inf)))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.34887 -0.25847 -0.00133  0.03518  0.71166 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  -9.7582     0.8571 -11.386  1.2e-06 ***
## mean          0.6666     0.1869   3.567  0.00605 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3614 on 9 degrees of freedom
## Multiple R-squared:  0.5857,	Adjusted R-squared:  0.5397 
## F-statistic: 12.72 on 1 and 9 DF,  p-value: 0.006052
```

```r
# Plot 
shannon_vs_fracprod_percell <- ggplot(filter(ML_otu_shannon_stats, fracprod_per_cell != Inf),
       aes(x=mean, y=log10(fracprod_per_cell), color = fraction)) + 
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") + 
  xlab("Shannon Entropy") +
  scale_x_continuous(limits = c(3.4, 5.85), breaks = c(3.5, 4, 4.5, 5, 5.5)) + 
  geom_smooth(data=subset(ML_otu_shannon_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 4.75, y=-8, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_percell_ML_otu_shannon)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(freeprod_percell_ML_otu_shannon)$coefficients[,4][2]), digits = 2))) + 
  annotate("text", x = 5.5, y=-7, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_percell_MLotu_shannon)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(partprod_percell_MLotu_shannon)$coefficients[,4][2]), digits = 3)));


#########################################################  INVERSE SIMPSON 
# Free-Living Inverse Simpson vs fractional production per cell 
freeprod_percell_ML_otu_invsimps <- lm(log10(fracprod_per_cell) ~ mean, 
                                       data = filter(ML_otu_invsimps_stats, fraction == "WholeFree"))
summary(freeprod_percell_ML_otu_invsimps)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(ML_otu_invsimps_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.67727 -0.13613  0.01755  0.18500  0.59968 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -8.05281    0.33073 -24.349 3.11e-10 ***
## mean         0.01923    0.01258   1.528    0.157    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3729 on 10 degrees of freedom
## Multiple R-squared:  0.1894,	Adjusted R-squared:  0.1083 
## F-statistic: 2.336 on 1 and 10 DF,  p-value: 0.1574
```

```r
# Particle-Associated Inverse Simpson vs fractional production per cell 
partprod_percell_MLotu_invsimps <- lm(log10(fracprod_per_cell)  ~ mean, 
                              data = filter(ML_otu_invsimps_stats, fraction == "WholePart" & fracprod_per_cell != Inf))
summary(partprod_percell_MLotu_invsimps)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(ML_otu_invsimps_stats, 
##     fraction == "WholePart" & fracprod_per_cell != Inf))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.29800 -0.18081 -0.11663  0.08102  0.56016 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -7.340126   0.157481 -46.610 4.82e-12 ***
## mean         0.016480   0.003462   4.761  0.00103 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.2993 on 9 degrees of freedom
## Multiple R-squared:  0.7158,	Adjusted R-squared:  0.6842 
## F-statistic: 22.66 on 1 and 9 DF,  p-value: 0.001029
```

```r
# Plot Simpson's Evenness
invsimps_vs_fracprod_percell <- ggplot(filter(ML_otu_invsimps_stats, fracprod_per_cell != Inf),
       aes(x=mean, y=log10(fracprod_per_cell) , color = fraction)) + 
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  scale_x_continuous(limits = c(0,95), breaks = c(0, 20, 40, 60, 80), expand = c(0,0)) + 
  ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") + 
  xlab("Inverse Simpson") +
  geom_smooth(data=subset(ML_otu_invsimps_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 50, y=-8, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_percell_ML_otu_invsimps)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(freeprod_percell_ML_otu_invsimps)$coefficients[,4][2]), digits = 2))) + 
  annotate("text", x = 75, y=-7, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_percell_MLotu_invsimps)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(partprod_percell_MLotu_invsimps)$coefficients[,4][2]), digits = 3))); 



#########################################################  SIMPSON'S EVENNESS
# Free-Living Inverse Simpson vs fractional production per cell 
freeprod_percell_ML_otu_simpseven <- lm(log10(fracprod_per_cell) ~ mean, 
                                       data = filter(ML_otu_simpseven_stats, fraction == "WholeFree"))
summary(freeprod_percell_ML_otu_simpseven)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(ML_otu_simpseven_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.6446 -0.2578  0.0211  0.1669  0.7773 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   -7.823      0.539 -14.513  4.8e-08 ***
## mean           2.801      5.936   0.472    0.647    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.4096 on 10 degrees of freedom
## Multiple R-squared:  0.02178,	Adjusted R-squared:  -0.07604 
## F-statistic: 0.2227 on 1 and 10 DF,  p-value: 0.6471
```

```r
# Particle-Associated Inverse Simpson vs fractional production per cell 
partprod_percell_MLotu_simpseven <- lm(log10(fracprod_per_cell)  ~ mean, 
                              data = filter(ML_otu_simpseven_stats, fraction == "WholePart" & fracprod_per_cell != Inf))
summary(partprod_percell_MLotu_simpseven)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(ML_otu_simpseven_stats, 
##     fraction == "WholePart" & fracprod_per_cell != Inf))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.45099 -0.20370 -0.08701  0.11225  0.58419 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  -7.6709     0.2715 -28.253 4.24e-10 ***
## mean         12.9548     3.4302   3.777  0.00437 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3492 on 9 degrees of freedom
## Multiple R-squared:  0.6131,	Adjusted R-squared:  0.5701 
## F-statistic: 14.26 on 1 and 9 DF,  p-value: 0.004371
```

```r
# Plot Simpson's Evenness
simpseven_vs_fracprod_percell <- ggplot(filter(ML_otu_simpseven_stats, fracprod_per_cell != Inf),
       aes(x=mean, y=log10(fracprod_per_cell) , color = fraction)) + 
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), 
                     limits = c("WholePart", "WholeFree"),
                     breaks=c("WholePart", "WholeFree"),
                     labels=c("Particle-Associated", "Free-Living")) + 
  scale_x_continuous(limits = c(0.03,0.151), breaks = c(0.05, 0.1, 0.15))  +
  ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") + 
  xlab("Simpson's Evenness") +
  geom_smooth(data=subset(ML_otu_simpseven_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 0.11, y=-8, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_percell_ML_otu_simpseven)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(freeprod_percell_ML_otu_simpseven)$coefficients[,4][2]), digits = 2))) + 
  annotate("text", x =0.13, y=-7, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_percell_MLotu_simpseven)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(partprod_percell_MLotu_simpseven)$coefficients[,4][2]), digits = 3))); 




plot_grid(rich_vs_fracprod_percell + theme(legend.position= "none"), 
          shannon_vs_fracprod_percell + theme(legend.position= "none"),  
          invsimps_vs_fracprod_percell + theme(legend.position= "none"),
          simpseven_vs_fracprod_percell + theme(legend.position= c(0.35,0.9)),
          labels = c("A", "B", "C", "D"), 
          ncol = 4)
```

<img src="Rarefied_Figures/fracprod_percell-vs-div-1.png" style="display: block; margin: auto;" />



# Fraction Production Altogether

```r
#########################################################  RICHNESS 
# Free-Living Richness vs fractional production per cell 
total_prodpercell_rich <- lm(log10(fracprod_per_cell) ~ mean, data = filter(filter(ML_otu_rich_stats, fracprod_per_cell != Inf)))
summary(total_prodpercell_rich)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(filter(ML_otu_rich_stats, 
##     fracprod_per_cell != Inf)))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -0.8487 -0.2245  0.0577  0.2207  0.7277 
## 
## Coefficients:
##               Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -8.4027594  0.2145671 -39.161  < 2e-16 ***
## mean         0.0033662  0.0005434   6.195 3.81e-06 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3826 on 21 degrees of freedom
## Multiple R-squared:  0.6463,	Adjusted R-squared:  0.6295 
## F-statistic: 38.38 on 1 and 21 DF,  p-value: 3.808e-06
```

```r
# Plot 
combined_richness <- ggplot(filter(ML_otu_rich_stats, fracprod_per_cell != Inf), aes(x=mean, y=log10(fracprod_per_cell))) + 
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") + 
  xlab("Observed Richness") +
  geom_smooth(data= ML_otu_rich_stats, method='lm', color = "black") + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 500, y=-8, color = "black", fontface = "bold",
           label = paste("R2 =", round(summary(total_prodpercell_rich)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(total_prodpercell_rich)$coefficients[,4][2]), digits = 6)))



#########################################################  SHANNON
total_prodpercell_shannon <- lm(log10(fracprod_per_cell) ~ mean, data = filter(filter(ML_otu_shannon_stats,fracprod_per_cell != Inf)))
summary(total_prodpercell_shannon)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(filter(ML_otu_shannon_stats, 
##     fracprod_per_cell != Inf)))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.01023 -0.15351 -0.02587  0.23804  0.84872 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -10.8281     0.7229 -14.979 1.11e-12 ***
## mean          0.8575     0.1681   5.101 4.72e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.4299 on 21 degrees of freedom
## Multiple R-squared:  0.5534,	Adjusted R-squared:  0.5322 
## F-statistic: 26.02 on 1 and 21 DF,  p-value: 4.719e-05
```

```r
# Plot 
combined_shannon <- ggplot(filter(ML_otu_shannon_stats, fracprod_per_cell != Inf),
       aes(x=mean, y=log10(fracprod_per_cell))) + 
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") + 
  xlab("Shannon Entropy") +
  geom_smooth(data= ML_otu_shannon_stats, method='lm', color = "black") + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 4.75, y=-8, color = "black", fontface = "bold",
           label = paste("R2 =", round(summary(total_prodpercell_shannon)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(total_prodpercell_shannon)$coefficients[,4][2]), digits = 6)));


#########################################################  INVERSE SIMPSON 
total_prodpercell_simpson <- lm(log10(fracprod_per_cell)  ~ mean, data = filter(ML_otu_invsimps_stats, fracprod_per_cell != Inf))
summary(total_prodpercell_simpson)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(ML_otu_invsimps_stats, 
##     fracprod_per_cell != Inf))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.95760 -0.26625 -0.04747  0.30807  0.98598 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -7.837452   0.173634 -45.138  < 2e-16 ***
## mean         0.021713   0.004725   4.595 0.000157 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.4542 on 21 degrees of freedom
## Multiple R-squared:  0.5014,	Adjusted R-squared:  0.4776 
## F-statistic: 21.11 on 1 and 21 DF,  p-value: 0.0001567
```

```r
# Plot Simpson's Evenness
combined_invsimps <- ggplot(filter(ML_otu_invsimps_stats, fracprod_per_cell != Inf),
       aes(x=mean, y=log10(fracprod_per_cell))) + 
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) + 
  ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") + 
  xlab("Inverse Simpson") +
  geom_smooth(data = ML_otu_invsimps_stats, method='lm', color = "black") + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 50, y=-8, color = "black", fontface = "bold",
           label = paste("R2 =", round(summary(total_prodpercell_simpson)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(total_prodpercell_simpson)$coefficients[,4][2]), digits = 5))) 


#########################################################  SIMPSON'S EVENNESS
total_prodpercell_simpseven <- lm(log10(fracprod_per_cell) ~ mean, 
                                       data = filter(ML_otu_simpseven_stats, fracprod_per_cell != Inf))
summary(total_prodpercell_simpseven)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ mean, data = filter(ML_otu_simpseven_stats, 
##     fracprod_per_cell != Inf))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -1.07904 -0.43066  0.07355  0.35519  1.50111 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  -7.5192     0.4198 -17.912 3.36e-14 ***
## mean          4.3198     4.9137   0.879    0.389    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.6318 on 21 degrees of freedom
## Multiple R-squared:  0.0355,	Adjusted R-squared:  -0.01043 
## F-statistic: 0.7729 on 1 and 21 DF,  p-value: 0.3893
```

```r
# Plot Simpson's Evenness
combined_simpseven <- ggplot(filter(ML_otu_simpseven_stats, fracprod_per_cell != Inf),
       aes(x=mean, y=log10(fracprod_per_cell))) + 
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  scale_x_continuous(limits = c(0.03,0.151), breaks = c(0.05, 0.1, 0.15))  +
  ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") + 
  xlab("Simpson's Evenness") +
  geom_smooth(data = ML_otu_simpseven_stats, method='lm', color = "black") + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 0.11, y=-8, color = "black", fontface = "bold",
           label = paste("R2 =", round(summary(total_prodpercell_simpseven)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(total_prodpercell_simpseven)$coefficients[,4][2]), digits = 2)))



plot_grid(rich_vs_fracprod_percell + ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") + theme(legend.position = "none"), 
          shannon_vs_fracprod_percell + ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)")+ theme(legend.position = "none"), 
          invsimps_vs_fracprod_percell + ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)")+ theme(legend.position = "none"), 
          simpseven_vs_fracprod_percell + ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") +
                    theme(legend.position = c(0.35,0.9)), 
          combined_richness, combined_shannon, combined_invsimps,combined_simpseven, 
          labels = c("A", "B", "C", "D","E", "F", "G", "H"), 
          ncol = 4, nrow = 2)
```

<img src="Rarefied_Figures/fracprod-together-vs-div-1.png" style="display: block; margin: auto;" />



# Exponential Shannon

```r
#########################################################  EXPONENTIALSHANNON
# Free-Living EXPONENTIAL Shannon vs fractional production per cell 
freeprod_expshannon <- lm(frac_bacprod ~ exp(mean), data = filter(ML_otu_shannon_stats, fraction == "WholeFree"))
summary(freeprod_expshannon)
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ exp(mean), data = filter(ML_otu_shannon_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -21.831 -10.032  -4.348   7.285  35.288 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)  11.9560    16.6699   0.717    0.490
## exp(mean)     0.2094     0.2740   0.764    0.462
## 
## Residual standard error: 17.87 on 10 degrees of freedom
## Multiple R-squared:  0.05518,	Adjusted R-squared:  -0.0393 
## F-statistic: 0.5841 on 1 and 10 DF,  p-value: 0.4624
```

```r
# Particle-Associated EXPONENTIAL Shannon vs fractional production per cell 
partprod_expshannon <- lm(frac_bacprod ~ exp(mean),
                                        data = filter(filter(ML_otu_shannon_stats, fraction == "WholePart" & fracprod_per_cell != Inf)))
summary(partprod_expshannon)
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ exp(mean), data = filter(filter(ML_otu_shannon_stats, 
##     fraction == "WholePart" & fracprod_per_cell != Inf)))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -7.6328 -2.7344 -0.9818  1.9563 12.3093 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)   
## (Intercept)  0.16204    3.29910   0.049  0.96190   
## exp(mean)    0.08538    0.02449   3.486  0.00687 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.936 on 9 degrees of freedom
## Multiple R-squared:  0.5746,	Adjusted R-squared:  0.5273 
## F-statistic: 12.16 on 1 and 9 DF,  p-value: 0.006869
```

```r
# Plot 
expshannon_fracprod <- ggplot(filter(ML_otu_shannon_stats, fracprod_per_cell != Inf),
       aes(x= exp(mean), y=frac_bacprod, color = fraction)) + 
  geom_point(size = 3.5) + #geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  ylab("Fraction Production (μgC/L/hr)") + 
  xlab("Exponential Shannon") +
  #scale_x_continuous(limits = c(3.4, 5.85), breaks = c(3.5, 4, 4.5, 5, 5.5)) + 
  geom_smooth(data=subset(ML_otu_shannon_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 150, y=35, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_expshannon)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(freeprod_expshannon)$coefficients[,4][2]), digits = 2))) + 
  annotate("text", x = 200, y=7, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_expshannon)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(partprod_expshannon)$coefficients[,4][2]), digits = 3)));




#########################################################  EXPONENTIALSHANNON
# Free-Living EXPONENTIAL Shannon vs fractional production per cell 
freeprod_percell_ML_otu_expshannon <- lm(log10(fracprod_per_cell) ~ exp(mean), data = filter(ML_otu_shannon_stats, fraction == "WholeFree"))
summary(freeprod_percell_ML_otu_expshannon)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ exp(mean), data = filter(ML_otu_shannon_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.71802 -0.15455  0.07673  0.16850  0.61743 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -7.972211   0.362962 -21.964 8.57e-10 ***
## exp(mean)    0.006868   0.005965   1.151    0.276    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.3892 on 10 degrees of freedom
## Multiple R-squared:  0.117,	Adjusted R-squared:  0.02875 
## F-statistic: 1.326 on 1 and 10 DF,  p-value: 0.2764
```

```r
# Particle-Associated EXPONENTIAL Shannon vs fractional production per cell 
partprod_percell_MLotu_expshannon <- lm(log10(fracprod_per_cell) ~ exp(mean), 
                                     data = filter(filter(ML_otu_shannon_stats, fraction == "WholePart" & fracprod_per_cell != Inf)))
summary(partprod_percell_MLotu_expshannon)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ exp(mean), data = filter(filter(ML_otu_shannon_stats, 
##     fraction == "WholePart" & fracprod_per_cell != Inf)))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.37450 -0.20313 -0.12798  0.08723  0.71521 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -7.337215   0.196190  -37.40 3.47e-11 ***
## exp(mean)    0.005403   0.001456    3.71  0.00484 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.353 on 9 degrees of freedom
## Multiple R-squared:  0.6047,	Adjusted R-squared:  0.5607 
## F-statistic: 13.77 on 1 and 9 DF,  p-value: 0.004844
```

```r
# Plot 
expshannon_fracprod_percell <- ggplot(filter(ML_otu_shannon_stats, fracprod_per_cell != Inf),
       aes(x= exp(mean), y=log10(fracprod_per_cell), color = fraction)) + 
  geom_point(size = 3.5) + #geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") + 
  xlab("Exponential Shannon") +
  #scale_x_continuous(limits = c(3.4, 5.85), breaks = c(3.5, 4, 4.5, 5, 5.5)) + 
  geom_smooth(data=subset(ML_otu_shannon_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 150, y=-8, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_percell_ML_otu_expshannon)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(freeprod_percell_ML_otu_expshannon)$coefficients[,4][2]), digits = 2))) + 
  annotate("text", x = 200, y=-7, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_percell_MLotu_expshannon)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(partprod_percell_MLotu_expshannon)$coefficients[,4][2]), digits = 3)));


#########################################################  EXPONENTIALSHANNON
combined_expshannon_percell <- lm(log10(fracprod_per_cell) ~ exp(mean), data = filter(filter(ML_otu_shannon_stats, fracprod_per_cell != Inf)))
summary(combined_expshannon_percell)
```

```
## 
## Call:
## lm(formula = log10(fracprod_per_cell) ~ exp(mean), data = filter(filter(ML_otu_shannon_stats, 
##     fracprod_per_cell != Inf)))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -0.92928 -0.21787  0.01261  0.28579  0.84384 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept) -7.813569   0.157920 -49.478  < 2e-16 ***
## exp(mean)    0.007648   0.001533   4.989 6.15e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.4352 on 21 degrees of freedom
## Multiple R-squared:  0.5424,	Adjusted R-squared:  0.5206 
## F-statistic: 24.89 on 1 and 21 DF,  p-value: 6.153e-05
```

```r
# Plot 
expshannon_combined <- ggplot(filter(ML_otu_shannon_stats, fracprod_per_cell != Inf),
       aes(x= exp(mean), y=log10(fracprod_per_cell))) + 
  geom_point(size = 3.5) + 
  ylab("log10(Fraction Production/Cell)\n (μgC/cell/hr)") + 
  xlab("Exponential Shannon") +
  #scale_x_continuous(limits = c(3.4, 5.85), breaks = c(3.5, 4, 4.5, 5, 5.5)) + 
  geom_smooth(data= ML_otu_shannon_stats, method='lm', color = "black") + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 150, y=-8, color = "black", fontface = "bold",
           label = paste("R2 =", round(summary(combined_expshannon_percell)$adj.r.squared, digits = 2), "\n", 
                         "p =", round(unname(summary(combined_expshannon_percell)$coefficients[,4][2]), digits = 5))) 

plot_grid(expshannon_fracprod, expshannon_fracprod_percell, expshannon_combined,
          labels = c("A", "B", "C"), ncol = 3)
```

<img src="Rarefied_Figures/exp-shannon-1.png" style="display: block; margin: auto;" />



## Congruency between fractions 

```r
rich_df <- ML_otu_rich_stats %>%
  dplyr::select(norep_filter_name, mean, sd, measure, fraction) %>%
  mutate(norep_water_name = paste(substr(norep_filter_name, 1,4), substr(norep_filter_name, 6,9), sep = "")) %>%
  spread(fraction, mean) 

rich_part_df <- dplyr::select(rich_df, norep_water_name, WholePart) %>%
  distinct() %>%  filter(!is.na(WholePart))

rich_free_df <- dplyr::select(rich_df, norep_water_name, WholeFree) %>%
  distinct() %>%  filter(!is.na(WholeFree))

combined_rich_df <- left_join(rich_free_df, rich_part_df, by = "norep_water_name") %>%
  rename(names = norep_water_name) %>%
  make_metadata_norep() %>%
  dplyr::select(-c(year, fraction, month, season, nuc_acid_type, project))


summary(lm(WholePart ~ WholeFree, data = combined_rich_df))
```

```
## 
## Call:
## lm(formula = WholePart ~ WholeFree, data = combined_rich_df)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -122.70 -106.94  -29.13   45.48  355.39 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept) 271.8734   157.9006   1.722    0.116
## WholeFree     0.6495     0.5407   1.201    0.257
## 
## Residual standard error: 146 on 10 degrees of freedom
## Multiple R-squared:  0.1261,	Adjusted R-squared:  0.03868 
## F-statistic: 1.443 on 1 and 10 DF,  p-value: 0.2574
```

```r
ggplot(combined_rich_df, aes(x=WholePart, y=WholeFree)) + 
  geom_point(size = 3.5) +
  geom_abline(intercept = 0, slope = 1)
```

<img src="Rarefied_Figures/frac-congruency-1.png" style="display: block; margin: auto;" />


<!-- # Total Production  --> 



<!-- # Prefiltered Fraction Diversity-Production Analysis --> 



<!-- # Does DNA extraction concentration influcence diversity? --> 



