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

poster_rich1 <- ggplot(filter(otu_alphadiv, 
                          fraction %in% c("WholePart", "WholeFree") & measure == "Richness"), 
       aes(y = mean, x = fraction)) +
  scale_color_manual(values = fraction_colors) + 
  scale_fill_manual(values = fraction_colors) +
  geom_jitter(size = 3, aes(color = fraction, fill = fraction)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  xlab("") + ylab("Observed Richness") +  theme(legend.position = "none") +
  geom_path(data = nums1, aes(x = a, y = b), linetype = 1, color = "gray40") 

poster_rich <- poster_rich1 + 
  annotate("text", x=1.5, y=800, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("***\np =", round(rich_wilcox$p.value, digits = 3))) 



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
  scale_color_manual(values = fraction_colors) + 
  scale_fill_manual(values = fraction_colors) +
  ylab("Shannon Entropy") + xlab("") +
  geom_jitter(size = 3, aes(color = fraction, fill = fraction)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction))  +
  theme(legend.position = "none") + 
  geom_path(data = nums2, aes(x = a, y = b), linetype = 1, color = "gray40") 

poster_shannon <- poster_shannon1 + 
  annotate("text", x=1.5, y=5.8, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("***\np =", round(shannon_wilcox$p.value, digits = 3))) 


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
  scale_color_manual(values = fraction_colors) + 
  scale_fill_manual(values = fraction_colors) +
  ylab("Inverse Simpson") + xlab("") +
  geom_jitter(size = 3, aes(color = fraction, fill = fraction)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  geom_path(data = nums3, aes(x = a, y = b), linetype = "dotted", color = "gray40") 


  
poster_invsimps <- poster_invsimps1 +
  annotate("text", x=1.5, y=90, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("NS\np =", round(simpson_wilcox$p.value, digits = 3))) +
  theme(legend.position = "none")

plot_grid(poster_rich, poster_shannon, poster_invsimps,
          labels = c("A", "B", "C"),
          ncol = 3)
```

<img src="Rarefied_Figures/direct-comparison-1.png" style="display: block; margin: auto;" />


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
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  ylab("Production (μgC/L/hr)") + xlab("Observed Richness") +
  geom_smooth(data=subset(ML_otu_rich_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.15,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 500, y=45, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_ML_otu_rich)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(freeprod_ML_otu_rich)$coefficients[,4][2]), digits = 4))) + 
  annotate("text", x = 650, y=5, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_MLotu_rich)$adj.r.squared, digits = 4), "\n", 
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
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  ylab("Production (μgC/L/hr)") + xlab("Shannon Entropy") +
  geom_smooth(data=subset(ML_otu_shannon_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.15,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 4.5, y=45, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_ML_otu_shannon)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(freeprod_ML_otu_shannon)$coefficients[,4][2]), digits = 4))) + 
  annotate("text", x = 5.35, y=5, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_MLotu_shannon)$adj.r.squared, digits = 4), "\n", 
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
  geom_point(size = 3.5) + geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd)) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) + 
  #scale_y_continuous(limits = c(0,70),expand = c(0,0)) + 
  ylab("Production (μgC/L/hr)") + xlab("Inverse Simpson") +
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
  ggtitle("OTU: Vegan") +
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  scale_x_continuous(expand = c(0,0), limits=c(0, 0.15)) + 
  ylab("Production (μgC/L/hr)") + xlab("Simpson's Evenness") +
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
```

<img src="Rarefied_Figures/plot-vegan-fracprod-1.png" style="display: block; margin: auto;" />

```r
#ggsave("../Figures/veganRARE_otu_alpha_vs_prod.png", otu_vegan, dpi = 600, units = "in", width = 10, height = 8)
```




```r
plot_grid(poster_rich1 + xlab("Fraction") +theme(legend.position = c(0.85, 0.85), 
                              axis.text.y = element_blank(), 
                              legend.title = element_blank()) +coord_flip() +
                annotate("text", x=1.5, y=700, fontface = "bold",  size = 3.5, color = "gray40",
                          label= paste("***\np =", round(rich_wilcox$p.value, digits = 3))), 
          poster_shannon1 + xlab("Fraction") + theme(legend.position = c(0.85, 0.85), 
                              axis.text.y = element_blank(), 
                              legend.title = element_blank()) +coord_flip() +
                annotate("text", x=1.5, y=5.5, fontface = "bold",  size = 3.5, color = "gray40",
                          label= paste("***\np =", round(shannon_wilcox$p.value, digits = 3))),  
          poster_invsimps1 + xlab("Fraction") + theme(legend.position = c(0.85, 0.85), 
                              axis.text.y = element_blank(), 
                              legend.title = element_blank()) +coord_flip() +
                annotate("text", x=1.5, y=80, fontface = "bold",  size = 3.5, color = "gray40",
                          label= paste("NS\np =", round(simpson_wilcox$p.value, digits = 3))), 
          otu_rich_vegan, otu_shannon_vegan, otu_invsimps_vegan,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 3, nrow = 2)
```

<img src="Rarefied_Figures/unnamed-chunk-1-1.png" style="display: block; margin: auto;" />




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
  scale_color_manual(values = fraction_colors) + xlab("") +
  scale_fill_manual(values = fraction_colors) +
  geom_jitter(size = 3, aes(color = fraction, fill = fraction)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  ylab("Log10(Bacterial Cells/mL)") +
  #####  WHOLE PARTICLE VS WHOLE FREE CELL ABUNDANCES
  geom_path(data = dat1, aes(x = a, y = b), linetype = 1, color = "gray40") +
  annotate("text", x=1.5, y=6.5, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("***\np =", round(frac_abund_wilcox$p.value, digits = 7))) +
  theme(legend.position = "none")



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
       aes(y = frac_bacprod, x = fraction)) + xlab("") +
  scale_color_manual(values = fraction_colors) + 
  scale_fill_manual(values = fraction_colors) +
  geom_jitter(size = 3, aes(color = fraction, fill = fraction)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  ylab("Total Fraction production\n (ug C/L/hr)") +
  #####  WHOLE PARTICLE VS WHOLE FREE TOTAL PRODUCTION 
  geom_path(data = dat2, aes(x = a, y = b), linetype = 1, color = "gray40") +
  annotate("text", x=1.5, y=68, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("***\np =", round(totprod_wilcox$p.value, digits = 3))) +
  theme(legend.position = "none")



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
  geom_jitter(size = 3, aes(color = fraction, fill = fraction)) + 
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(color = fraction, fill = fraction)) +
  ylim(c(-8.5, -5)) + xlab("") +
  ylab("log10(Fraction Production per Cell) \n(ug C/cell/hr)") +
  #####  WHOLE PARTICLE VS WHOLE FREE PER CELL PRODUCTION 
  geom_path(data = dat3, aes(x = a, y = b), linetype = 1, color = "gray40") +
  annotate("text", x=1.5, y=-5, fontface = "bold",  size = 3.5, color = "gray40",
           label= paste("***\np =", round(percellprod_wilcox$p.value, digits = 5))) +
  theme(legend.position = "none")

plot_grid(poster_a, poster_b, poster_c,
          labels = c("A", "B", "C"),
          ncol = 3)
```

<img src="Rarefied_Figures/per-cell-1.png" style="display: block; margin: auto;" />


# Diversity vs Total Production 

```r
######################################################### RICHNESS
# Total Bacterial Production vs Free-Living Richness 
summary(lm(tot_bacprod ~ mean, data = filter(ML_otu_rich_stats, fraction == "WholeFree")))
```

```
## 
## Call:
## lm(formula = tot_bacprod ~ mean, data = filter(ML_otu_rich_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -28.114 -12.944   0.594   8.159  41.427 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept) -9.95411   23.04539  -0.432    0.675  
## mean         0.15434    0.07892   1.956    0.079 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 21.31 on 10 degrees of freedom
## Multiple R-squared:  0.2766,	Adjusted R-squared:  0.2043 
## F-statistic: 3.824 on 1 and 10 DF,  p-value: 0.07901
```

```r
# Total Bacterial Production vs Particle-Associated Richness 
summary(lm(tot_bacprod ~ mean, data = filter(ML_otu_rich_stats, fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = tot_bacprod ~ mean, data = filter(ML_otu_rich_stats, 
##     fraction == "WholePart"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -31.914 -15.358  -0.829   7.469  41.379 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept) -6.13872   20.29427  -0.302    0.768  
## mean         0.08714    0.04259   2.046    0.068 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 21.03 on 10 degrees of freedom
## Multiple R-squared:  0.2951,	Adjusted R-squared:  0.2246 
## F-statistic: 4.186 on 1 and 10 DF,  p-value: 0.06797
```

```r
# Plot it 
ggplot(ML_otu_rich_stats, aes(x = mean, y = tot_bacprod)) + 
  geom_point(size = 3) + xlab("Richness") +
  facet_grid(.~fraction)
```

<img src="Rarefied_Figures/diversity-totalprod-1.png" style="display: block; margin: auto;" />

```r
######################################################### SHANNON ENTROPY
# Total Bacterial Production vs Free-Living Shannon 
summary(lm(tot_bacprod ~ mean, data = filter(ML_otu_shannon_stats, fraction == "WholeFree")))
```

```
## 
## Call:
## lm(formula = tot_bacprod ~ mean, data = filter(ML_otu_shannon_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -30.573 -11.102  -2.162   7.367  50.363 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)   -85.67      86.81  -0.987    0.347
## mean           29.72      21.59   1.376    0.199
## 
## Residual standard error: 22.97 on 10 degrees of freedom
## Multiple R-squared:  0.1593,	Adjusted R-squared:  0.07521 
## F-statistic: 1.895 on 1 and 10 DF,  p-value: 0.1987
```

```r
# Total Bacterial Production vs Particle-Associated Shannon 
summary(lm(tot_bacprod ~ mean, data = filter(ML_otu_shannon_stats, fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = tot_bacprod ~ mean, data = filter(ML_otu_shannon_stats, 
##     fraction == "WholePart"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -29.412 -15.260   0.078   5.226  42.319 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept)   -71.03      49.21  -1.443   0.1795  
## mean           22.94      10.72   2.140   0.0581 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 20.75 on 10 degrees of freedom
## Multiple R-squared:  0.3141,	Adjusted R-squared:  0.2455 
## F-statistic: 4.578 on 1 and 10 DF,  p-value: 0.05806
```

```r
# Plot it 
ggplot(ML_otu_shannon_stats, aes(x = mean, y = tot_bacprod)) + 
  geom_point(size = 3) + xlab("Shannon Entropy") +
  geom_smooth(method = "lm", se = FALSE, data = filter(ML_otu_shannon_stats, fraction == "WholePart")) +
  facet_grid(.~fraction)
```

<img src="Rarefied_Figures/diversity-totalprod-2.png" style="display: block; margin: auto;" />

```r
######################################################### INVERSE SIMPSON
# Total Bacterial Production vs Free-Living Inverse Simpson  
summary(lm(tot_bacprod ~ mean, data = filter(ML_otu_invsimps_stats, fraction == "WholeFree"))) 
```

```
## 
## Call:
## lm(formula = tot_bacprod ~ mean, data = filter(ML_otu_invsimps_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -26.050 -11.343  -4.267   4.427  51.512 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)   5.4207    20.1418   0.269    0.793
## mean          1.1289     0.7662   1.473    0.171
## 
## Residual standard error: 22.71 on 10 degrees of freedom
## Multiple R-squared:  0.1783,	Adjusted R-squared:  0.09618 
## F-statistic: 2.171 on 1 and 10 DF,  p-value: 0.1714
```

```r
# Total Bacterial Production vs Particle-Associated Inverse Simpson 
summary(lm(tot_bacprod ~ mean, data = filter(ML_otu_invsimps_stats, fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = tot_bacprod ~ mean, data = filter(ML_otu_invsimps_stats, 
##     fraction == "WholePart"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -31.643  -9.674  -2.796   5.421  30.770 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept)   9.8929     9.4884   1.043   0.3217  
## mean          0.6288     0.2105   2.986   0.0137 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 18.21 on 10 degrees of freedom
## Multiple R-squared:  0.4714,	Adjusted R-squared:  0.4185 
## F-statistic: 8.918 on 1 and 10 DF,  p-value: 0.01366
```

```r
ggplot(ML_otu_invsimps_stats, aes(x = mean, y = tot_bacprod)) + 
  geom_point(size = 3) + xlab("Inverse Simpson") +
  geom_smooth(method = "lm", data = filter(ML_otu_invsimps_stats, fraction == "WholePart")) +
  facet_grid(.~fraction)
```

<img src="Rarefied_Figures/diversity-totalprod-3.png" style="display: block; margin: auto;" />

```r
######################################################### SIMPSONS EVENNESS
# Total Bacterial Production vs Free-Living Simpsons Evenness 
summary(lm(tot_bacprod ~ mean, data = filter(ML_otu_simpseven_stats, fraction == "WholeFree")))
```

```
## 
## Call:
## lm(formula = tot_bacprod ~ mean, data = filter(ML_otu_simpseven_stats, 
##     fraction == "WholeFree"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -24.340 -16.571  -5.009   9.425  61.363 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)
## (Intercept)    25.85      32.87   0.786    0.450
## mean           86.07     362.02   0.238    0.817
## 
## Residual standard error: 24.98 on 10 degrees of freedom
## Multiple R-squared:  0.00562,	Adjusted R-squared:  -0.09382 
## F-statistic: 0.05652 on 1 and 10 DF,  p-value: 0.8169
```

```r
# Total Bacterial Production vs Particle-Associated Simpsons Evenness  
summary(lm(tot_bacprod ~ mean, data = filter(ML_otu_simpseven_stats, fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = tot_bacprod ~ mean, data = filter(ML_otu_simpseven_stats, 
##     fraction == "WholePart"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -19.776 -11.062  -4.714   6.358  34.205 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)  
## (Intercept)    -5.91      13.92  -0.425   0.6801  
## mean          523.76     171.64   3.051   0.0122 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 18.03 on 10 degrees of freedom
## Multiple R-squared:  0.4822,	Adjusted R-squared:  0.4304 
## F-statistic: 9.311 on 1 and 10 DF,  p-value: 0.01222
```

```r
# Plot it 
ggplot(ML_otu_simpseven_stats, aes(x = mean, y = tot_bacprod)) + 
  geom_point(size = 3) + xlab("Simpson's Evenness") +
  geom_smooth(method = "lm", data = filter(ML_otu_simpseven_stats, fraction == "WholePart")) +
  facet_grid(.~fraction)
```

<img src="Rarefied_Figures/diversity-totalprod-4.png" style="display: block; margin: auto;" />


# Per cell Fraction Production vs Diversity

```r
# Can fractional production be predicted by richness? 
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
# PARTICLE
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
  ggtitle("OTU: Vegan") +
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  ylab("log10(Fraction Production per Cell) \n(μg C/cell/hr)") +
  xlab("Observed Richness") +
  geom_smooth(data=subset(ML_otu_rich_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.15,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 500, y=-8, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_percell_ML_otu_rich)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(freeprod_percell_ML_otu_rich)$coefficients[,4][2]), digits = 4))) + 
  annotate("text", x = 650, y=-7, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_percell_MLotu_rich)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(partprod_percell_MLotu_rich)$coefficients[,4][2]), digits = 4)));


###### INVERSE SIMPSON
# Can fractional production be predicted by invsimpsness? 
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
## PARTICLE
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
  ggtitle("OTU: Vegan") +
  scale_color_manual(values = c("firebrick3","cornflowerblue"), limits = c("WholePart", "WholeFree")) +
  scale_x_continuous(limits = c(0,100), expand = c(0,0)) + 
  ylab("log10(Fraction Production per Cell) \n(μg C/cell/hr)") + 
  xlab("Inverse Simpson") +
  geom_smooth(data=subset(ML_otu_invsimps_stats, fraction == "WholePart"), method='lm') + 
  theme(legend.position=c(0.2,0.9),        
        legend.title=element_blank()) +
  annotate("text", x = 40, y=-8, color = "cornflowerblue", fontface = "bold",
           label = paste("R2 =", round(summary(freeprod_percell_ML_otu_invsimps)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(freeprod_percell_ML_otu_invsimps)$coefficients[,4][2]), digits = 4))) + 
  annotate("text", x = 75, y=-7, color = "firebrick3", fontface = "bold",
           label = paste("R2 =", round(summary(partprod_percell_MLotu_invsimps)$adj.r.squared, digits = 4), "\n", 
                         "p =", round(unname(summary(partprod_percell_MLotu_invsimps)$coefficients[,4][2]), digits = 4))); 


plot_grid(rich_vs_fracprod_percell, invsimps_vs_fracprod_percell,
          labels = c("A", "B"), 
          ncol = 2)
```

<img src="Rarefied_Figures/fracprod_percell-vs-div-1.png" style="display: block; margin: auto;" />




# Prefiltered Fraction Diversity-Production Analysis 




# Does DNA extraction concentration influcence diversity?



