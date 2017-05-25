# OTU Removal Analysis
Marian L. Schmidt  
May 1st, 2017  
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
library(sandwich) # for vcovHC function in post-hoc test
library(MASS) # For studres in plot_residuals function
library(boot) # For cross validation
source("Muskegon_functions.R")
source("set_colors.R")
```


# Load data

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
# Productivity measurements are reliable only up to 1 decimal
df1 <- sample_data(otu_merged_musk_pruned) %>% 
  dplyr::mutate(tot_bacprod = round(tot_bacprod, digits = 1),
                SD_tot_bacprod = round(SD_tot_bacprod, digits = 1),
                frac_bacprod = round(frac_bacprod, digits = 1),
                SD_frac_bacprod = round(SD_frac_bacprod, digits = 1),
                fraction_bac_abund = as.numeric(fraction_bac_abund),
                fracprod_per_cell = frac_bacprod/(1000*fraction_bac_abund),
                fracprod_per_cell_noinf = ifelse(fracprod_per_cell == Inf, NA, fracprod_per_cell)) %>%
  dplyr::select(norep_filter_name, lakesite, limnion, fraction, year, season, tot_bacprod, SD_tot_bacprod, frac_bacprod, SD_frac_bacprod, fraction_bac_abund, fracprod_per_cell, fracprod_per_cell_noinf)
```

```
## Warning in class(x) <- c("tbl_df", "tbl", "data.frame"): Setting class(x) to multiple strings ("tbl_df", "tbl", ...); result will no longer be an S4 object
```

```r
row.names(df1) = df1$norep_filter_name
# Add new sample data back into phyloseq object 
sample_data(otu_merged_musk_pruned) <- df1

# Remove MOTHJ715 and MBRHP715 because of low sequencing depth 
otu_merged_musk_pruned_noMOTHJ715_MBRHP715 <- subset_samples(otu_merged_musk_pruned, norep_filter_name != "MOTHJ715" & norep_filter_name != "MBRHP715")
otu_merged_musk_pruned_noMOTHJ715_MBRHP715
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 52980 taxa and 161 samples ]
## sample_data() Sample Data:       [ 161 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 52980 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 52980 tips and 52978 internal nodes ]
```

```r
# Subset only the surface samples for the current study!!  
musk_surface <- subset_samples(otu_merged_musk_pruned_noMOTHJ715_MBRHP715, 
                               limnion == "Top" & year == "2015" & 
                                 fraction %in% c("WholePart","WholeFree")) # Surface samples, 2015, and WholePart/WholeFree samples only!
musk_surface_pruned <- prune_taxa(taxa_sums(musk_surface) > 0, musk_surface) 
musk_surface_pruned
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 7806 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 7806 taxa by 8 taxonomic ranks ]
## phy_tree()    Phylogenetic Tree: [ 7806 tips and 7804 internal nodes ]
```

```r
# Remove tree
notree_musk_surface_pruned <- phyloseq(tax_table(musk_surface_pruned), otu_table(musk_surface_pruned), sample_data(musk_surface_pruned))

# Remove singletons!
musk_surface_pruned_rm1 <- prune_taxa(taxa_sums(notree_musk_surface_pruned) > 1, notree_musk_surface_pruned) 
musk_surface_pruned_rm1
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 4522 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 4522 taxa by 8 taxonomic ranks ]
```

```r
# Remove the tree for less computationally intensive steps
notree_musk_surface_pruned_rm1 <- phyloseq(tax_table(musk_surface_pruned_rm1), otu_table(musk_surface_pruned_rm1), sample_data(musk_surface_pruned_rm1))

# If taxa with 2 counts are removed 
prune_taxa(taxa_sums(notree_musk_surface_pruned_rm1) > 2, notree_musk_surface_pruned_rm1) 
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 2979 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 2979 taxa by 8 taxonomic ranks ]
```

```r
# If taxa with 5 counts are removed 
notree_musk_surface_pruned_rm5 <- prune_taxa(taxa_sums(notree_musk_surface_pruned) > 5, notree_musk_surface_pruned) 
notree_musk_surface_pruned_rm5
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1655 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 1655 taxa by 8 taxonomic ranks ]
```

```r
# If taxa with 10 counts are removed 
notree_musk_surface_pruned_rm10 <- prune_taxa(taxa_sums(notree_musk_surface_pruned) > 10, notree_musk_surface_pruned) 
notree_musk_surface_pruned_rm10
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 1059 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 1059 taxa by 8 taxonomic ranks ]
```

```r
# If taxa with 30 counts are removed 
notree_musk_surface_pruned_rm30 <- prune_taxa(taxa_sums(notree_musk_surface_pruned) > 30, notree_musk_surface_pruned) 
notree_musk_surface_pruned_rm30
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 531 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 531 taxa by 8 taxonomic ranks ]
```

```r
# If taxa with 60 counts are removed 
notree_musk_surface_pruned_rm60 <- prune_taxa(taxa_sums(notree_musk_surface_pruned) > 60, notree_musk_surface_pruned) 
notree_musk_surface_pruned_rm60
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 372 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 372 taxa by 8 taxonomic ranks ]
```

```r
# If taxa with 90 counts are removed 
notree_musk_surface_pruned_rm90 <- prune_taxa(taxa_sums(notree_musk_surface_pruned) > 90, notree_musk_surface_pruned) 
notree_musk_surface_pruned_rm90
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 293 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 293 taxa by 8 taxonomic ranks ]
```

```r
# If taxa with 150 counts are removed 
notree_musk_surface_pruned_rm150 <- prune_taxa(taxa_sums(notree_musk_surface_pruned) > 150, notree_musk_surface_pruned) 
notree_musk_surface_pruned_rm150
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 212 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 212 taxa by 8 taxonomic ranks ]
```

```r
# If taxa with 300 counts are removed 
notree_musk_surface_pruned_rm225 <- prune_taxa(taxa_sums(notree_musk_surface_pruned) > 225, notree_musk_surface_pruned) 
notree_musk_surface_pruned_rm225
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 156 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 156 taxa by 8 taxonomic ranks ]
```

```r
# If taxa with 300 counts are removed 
notree_musk_surface_pruned_rm300 <- prune_taxa(taxa_sums(notree_musk_surface_pruned) > 300, notree_musk_surface_pruned) 
notree_musk_surface_pruned_rm300
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 132 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 132 taxa by 8 taxonomic ranks ]
```


# Calculate Diversity 


```r
set.seed(777)

################## Remove singltons 
alpha_rm1 <- calc_alpha_diversity(physeq = notree_musk_surface_pruned_rm1)
min(sample_sums(notree_musk_surface_pruned_rm1)) - 1
```

```
## [1] 6589
```

```r
otu_alphadiv_rm1 <- calc_mean_alphadiv(physeq = notree_musk_surface_pruned_rm1,
                           richness_df = alpha_rm1$Richness, 
                           evenness_df = alpha_rm1$Inverse_Simpson, 
                           shannon_df = alpha_rm1$Shannon) %>%
    mutate(fraction = factor(fraction, levels = c("WholePart", "Particle", "WholeFree", "Free")),
         lakesite = factor(lakesite,  levels = c("MOT", "MDP", "MBR", "MIN")),
         measure = factor(measure, levels = c("Richness", "Simpsons_Evenness", "Shannon_Entropy", "Inverse_Simpson")),
         Removed = "1-tons")  

################## Remove 5-tons 
alpha_rm5 <- calc_alpha_diversity(physeq = notree_musk_surface_pruned_rm5)
min(sample_sums(notree_musk_surface_pruned_rm5)) - 1
```

```
## [1] 6376
```

```r
otu_alphadiv_rm5 <- calc_mean_alphadiv(physeq = notree_musk_surface_pruned_rm5,
                           richness_df = alpha_rm5$Richness, 
                           evenness_df = alpha_rm5$Inverse_Simpson, 
                           shannon_df = alpha_rm5$Shannon) %>%
    mutate(fraction = factor(fraction, levels = c("WholePart", "Particle", "WholeFree", "Free")),
         lakesite = factor(lakesite,  levels = c("MOT", "MDP", "MBR", "MIN")),
         measure = factor(measure, levels = c("Richness", "Simpsons_Evenness", "Shannon_Entropy", "Inverse_Simpson")),
         Removed = "5-tons")  


################## Remove 10-tons 
alpha_rm10 <- calc_alpha_diversity(physeq = notree_musk_surface_pruned_rm10)
min(sample_sums(notree_musk_surface_pruned_rm10)) - 1
```

```
## [1] 6246
```

```r
otu_alphadiv_rm10 <- calc_mean_alphadiv(physeq = notree_musk_surface_pruned_rm10,
                           richness_df = alpha_rm10$Richness, 
                           evenness_df = alpha_rm10$Inverse_Simpson, 
                           shannon_df = alpha_rm10$Shannon) %>%
    mutate(fraction = factor(fraction, levels = c("WholePart", "Particle", "WholeFree", "Free")),
         lakesite = factor(lakesite,  levels = c("MOT", "MDP", "MBR", "MIN")),
         measure = factor(measure, levels = c("Richness", "Simpsons_Evenness", "Shannon_Entropy", "Inverse_Simpson")),
         Removed = "10-tons") 

################## Remove 30-tons 
alpha_rm30 <- calc_alpha_diversity(physeq = notree_musk_surface_pruned_rm30)
min(sample_sums(notree_musk_surface_pruned_rm30)) - 1
```

```
## [1] 5945
```

```r
otu_alphadiv_rm30 <- calc_mean_alphadiv(physeq = notree_musk_surface_pruned_rm30,
                           richness_df = alpha_rm30$Richness, 
                           evenness_df = alpha_rm30$Inverse_Simpson, 
                           shannon_df = alpha_rm30$Shannon) %>%
    mutate(fraction = factor(fraction, levels = c("WholePart", "Particle", "WholeFree", "Free")),
         lakesite = factor(lakesite,  levels = c("MOT", "MDP", "MBR", "MIN")),
         measure = factor(measure, levels = c("Richness", "Simpsons_Evenness", "Shannon_Entropy", "Inverse_Simpson")),
         Removed = "30-tons")  



################## Remove 60-tons 
alpha_rm60 <- calc_alpha_diversity(physeq = notree_musk_surface_pruned_rm60)
min(sample_sums(notree_musk_surface_pruned_rm60)) - 1
```

```
## [1] 5802
```

```r
otu_alphadiv_rm60 <- calc_mean_alphadiv(physeq = notree_musk_surface_pruned_rm60,
                           richness_df = alpha_rm60$Richness, 
                           evenness_df = alpha_rm60$Inverse_Simpson, 
                           shannon_df = alpha_rm60$Shannon) %>%
    mutate(fraction = factor(fraction, levels = c("WholePart", "Particle", "WholeFree", "Free")),
         lakesite = factor(lakesite,  levels = c("MOT", "MDP", "MBR", "MIN")),
         measure = factor(measure, levels = c("Richness", "Simpsons_Evenness", "Shannon_Entropy", "Inverse_Simpson")),
         Removed = "60-tons")  


################## Remove 90-tons 
alpha_rm90 <- calc_alpha_diversity(physeq = notree_musk_surface_pruned_rm90)
min(sample_sums(notree_musk_surface_pruned_rm90)) - 1
```

```
## [1] 5662
```

```r
otu_alphadiv_rm90 <- calc_mean_alphadiv(physeq = notree_musk_surface_pruned_rm90,
                           richness_df = alpha_rm90$Richness, 
                           evenness_df = alpha_rm90$Inverse_Simpson, 
                           shannon_df = alpha_rm90$Shannon) %>%
    mutate(fraction = factor(fraction, levels = c("WholePart", "Particle", "WholeFree", "Free")),
         lakesite = factor(lakesite,  levels = c("MOT", "MDP", "MBR", "MIN")),
         measure = factor(measure, levels = c("Richness", "Simpsons_Evenness", "Shannon_Entropy", "Inverse_Simpson")),
         Removed = "90-tons")  


################## Remove 90-tons 
alpha_rm150 <- calc_alpha_diversity(physeq = notree_musk_surface_pruned_rm150)
min(sample_sums(notree_musk_surface_pruned_rm150)) - 1
```

```
## [1] 5425
```

```r
otu_alphadiv_rm150 <- calc_mean_alphadiv(physeq = notree_musk_surface_pruned_rm150,
                           richness_df = alpha_rm150$Richness, 
                           evenness_df = alpha_rm150$Inverse_Simpson, 
                           shannon_df = alpha_rm150$Shannon) %>%
    mutate(fraction = factor(fraction, levels = c("WholePart", "Particle", "WholeFree", "Free")),
         lakesite = factor(lakesite,  levels = c("MOT", "MDP", "MBR", "MIN")),
         measure = factor(measure, levels = c("Richness", "Simpsons_Evenness", "Shannon_Entropy", "Inverse_Simpson")),
         Removed = "150-tons")  

################## Remove 300-tons 
alpha_rm225 <- calc_alpha_diversity(physeq = notree_musk_surface_pruned_rm225)
min(sample_sums(notree_musk_surface_pruned_rm225)) - 1
```

```
## [1] 5233
```

```r
otu_alphadiv_rm225 <- calc_mean_alphadiv(physeq = notree_musk_surface_pruned_rm225,
                           richness_df = alpha_rm225$Richness, 
                           evenness_df = alpha_rm225$Inverse_Simpson, 
                           shannon_df = alpha_rm225$Shannon) %>%
    mutate(fraction = factor(fraction, levels = c("WholePart", "Particle", "WholeFree", "Free")),
         lakesite = factor(lakesite,  levels = c("MOT", "MDP", "MBR", "MIN")),
         measure = factor(measure, levels = c("Richness", "Simpsons_Evenness", "Shannon_Entropy", "Inverse_Simpson")),
         Removed = "225-tons")  


################## Remove 300-tons 
alpha_rm300 <- calc_alpha_diversity(physeq = notree_musk_surface_pruned_rm300)
min(sample_sums(notree_musk_surface_pruned_rm300)) - 1
```

```
## [1] 5061
```

```r
otu_alphadiv_rm300 <- calc_mean_alphadiv(physeq = notree_musk_surface_pruned_rm300,
                           richness_df = alpha_rm300$Richness, 
                           evenness_df = alpha_rm300$Inverse_Simpson, 
                           shannon_df = alpha_rm300$Shannon) %>%
    mutate(fraction = factor(fraction, levels = c("WholePart", "Particle", "WholeFree", "Free")),
         lakesite = factor(lakesite,  levels = c("MOT", "MDP", "MBR", "MIN")),
         measure = factor(measure, levels = c("Richness", "Simpsons_Evenness", "Shannon_Entropy", "Inverse_Simpson")),
         Removed = "300-tons")  
```


```r
min_seqs <- c(min(sample_sums(notree_musk_surface_pruned_rm1)) - 1, min(sample_sums(notree_musk_surface_pruned_rm5)) - 1, 
  min(sample_sums(notree_musk_surface_pruned_rm10)) - 1, min(sample_sums(notree_musk_surface_pruned_rm30)) - 1,
  min(sample_sums(notree_musk_surface_pruned_rm60)) - 1, min(sample_sums(notree_musk_surface_pruned_rm90)) - 1,
  min(sample_sums(notree_musk_surface_pruned_rm150)) - 1, min(sample_sums(notree_musk_surface_pruned_rm225)) - 1,
  min(sample_sums(notree_musk_surface_pruned_rm300)) - 1)


num_otus <- c(ncol(otu_table(notree_musk_surface_pruned_rm1)), ncol(otu_table(notree_musk_surface_pruned_rm5)), ncol(otu_table(notree_musk_surface_pruned_rm10)), 
              ncol(otu_table(notree_musk_surface_pruned_rm30)), ncol(otu_table(notree_musk_surface_pruned_rm60)), ncol(otu_table(notree_musk_surface_pruned_rm90)),
              ncol(otu_table(notree_musk_surface_pruned_rm150)), ncol(otu_table(notree_musk_surface_pruned_rm225)), ncol(otu_table(notree_musk_surface_pruned_rm300)))

Removed <- c("1-tons","5-tons", "10-tons", "30-tons", "60-tons", "90-tons", "150-tons", "225-tons","300-tons")


statz <- data.frame(cbind(as.numeric(min_seqs), as.numeric(num_otus), Removed)) %>%
         mutate(Removed = factor(Removed, levels = c("1-tons","5-tons", "10-tons", "30-tons", "60-tons", 
                                              "90-tons", "150-tons", "225-tons","300-tons"))) %>%
  mutate(num_otus = as.numeric(num_otus))

p1 <- ggplot(statz, aes(x = Removed, y = min_seqs, fill = Removed)) +
  geom_bar(stat = "identity") + ylab("Minimum Number of Sequences") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 8000)) +
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        legend.position = "none")
  
p2 <-ggplot(statz, aes(x = Removed, y = num_otus, fill = Removed)) +
  geom_bar(stat = "identity") + ylab("Richness") +
  scale_y_continuous(expand = c(0,0)) +
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        legend.position = c(0.9, 0.8), legend.title = element_blank())

plot_grid(p1, p2, align = "v", labels = c("A", "B"), nrow = 2, ncol =1)
```

<img src="OTU_Removal_Analysis_Figs/min-seqs-plots-1.png" style="display: block; margin: auto;" />






# Combine all data 

```r
# Combine all div metrics 
all_divs <- bind_rows(otu_alphadiv_rm1, otu_alphadiv_rm5, otu_alphadiv_rm10, otu_alphadiv_rm30, 
                      otu_alphadiv_rm60, otu_alphadiv_rm90, otu_alphadiv_rm150, 
                      otu_alphadiv_rm225, otu_alphadiv_rm300) %>%
  dplyr::filter(fraction %in% c("WholePart", "WholeFree") & year == "2015") %>%
  mutate(Removed = factor(Removed, levels = c("1-tons","5-tons", "10-tons", "30-tons", "60-tons", 
                                              "90-tons", "150-tons", "225-tons","300-tons")))
```


# Richness 

```r
### PLOT
ggplot(dplyr::filter(all_divs, measure == "Richness"), 
       aes(y = mean, x = Removed, color = Removed, fill = Removed)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) + geom_point(size = 3, position = position_jitter(w = 0.1)) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(.~fraction) +
  ylab("Mean Richness") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
```

<img src="OTU_Removal_Analysis_Figs/rich-plots1-1.png" style="display: block; margin: auto;" />

```r
ggplot(dplyr::filter(all_divs, measure == "Richness"), 
       aes(y = mean, x = fraction, color = Removed, fill = Removed)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) + geom_point(size = 3, position = position_jitter(w = 0.1)) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(.~Removed) +
  ylab("Mean Richness") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
```

<img src="OTU_Removal_Analysis_Figs/rich-plots1-2.png" style="display: block; margin: auto;" />

```r
summary(lm(frac_bacprod ~ mean + fraction, data = dplyr::filter(all_divs, measure == "Richness" & 
                                                       Removed == "1-tons")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean + fraction, data = dplyr::filter(all_divs, 
##     measure == "Richness" & Removed == "1-tons"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -20.797  -6.772  -0.585   6.031  33.411 
## 
## Coefficients:
##                    Estimate Std. Error t value Pr(>|t|)   
## (Intercept)       -13.35169   10.71466  -1.246   0.2264   
## mean                0.03843    0.01663   2.311   0.0311 * 
## fractionWholeFree  23.18701    6.44851   3.596   0.0017 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 12.52 on 21 degrees of freedom
## Multiple R-squared:  0.3815,	Adjusted R-squared:  0.3226 
## F-statistic: 6.476 on 2 and 21 DF,  p-value: 0.006446
```

```r
summary(lm(frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, measure == "Richness" & fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, 
##     measure == "Richness" & fraction == "WholePart"))
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -9.928 -5.493 -1.521  3.624 22.370 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## (Intercept)          -10.286510   5.124149  -2.007 0.047454 *  
## mean                   0.033280   0.008436   3.945 0.000150 ***
## mean:Removed5-tons     0.009687   0.006221   1.557 0.122657    
## mean:Removed10-tons    0.016918   0.007742   2.185 0.031268 *  
## mean:Removed30-tons    0.036689   0.012813   2.863 0.005127 ** 
## mean:Removed60-tons    0.052790   0.017170   3.075 0.002732 ** 
## mean:Removed90-tons    0.067999   0.021474   3.167 0.002057 ** 
## mean:Removed150-tons   0.092711   0.028425   3.262 0.001525 ** 
## mean:Removed225-tons   0.126222   0.037525   3.364 0.001098 ** 
## mean:Removed300-tons   0.150614   0.044357   3.395 0.000991 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 7.652 on 98 degrees of freedom
## Multiple R-squared:  0.1469,	Adjusted R-squared:  0.06859 
## F-statistic: 1.875 on 9 and 98 DF,  p-value: 0.06447
```

```r
# Are the slopes of Removed different from each other?
summary(lm(frac_bacprod ~ mean*Removed, data = dplyr::filter(all_divs, measure == "Richness" & fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean * Removed, data = dplyr::filter(all_divs, 
##     measure == "Richness" & fraction == "WholePart"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -12.7683  -4.2203  -0.9146   3.5894  20.3023 
## 
## Coefficients:
##                       Estimate Std. Error t value Pr(>|t|)   
## (Intercept)           -9.61507    7.03808  -1.366  0.17530   
## mean                   0.03227    0.01107   2.916  0.00448 **
## Removed5-tons         -8.36171   12.16940  -0.687  0.49378   
## Removed10-tons       -16.40576   14.92731  -1.099  0.27468   
## Removed30-tons       -27.68454   27.15250  -1.020  0.31065   
## Removed60-tons        -9.09932   33.71779  -0.270  0.78788   
## Removed90-tons        38.97631   34.72069   1.123  0.26461   
## Removed150-tons       69.89030   33.46261   2.089  0.03957 * 
## Removed225-tons       65.49632   38.04611   1.721  0.08860 . 
## Removed300-tons       86.75844   40.58064   2.138  0.03523 * 
## mean:Removed5-tons     0.02601    0.02307   1.128  0.26247   
## mean:Removed10-tons    0.05519    0.03347   1.649  0.10264   
## mean:Removed30-tons    0.12965    0.09025   1.437  0.15428   
## mean:Removed60-tons    0.08942    0.14010   0.638  0.52494   
## mean:Removed90-tons   -0.13008    0.17141  -0.759  0.44992   
## mean:Removed150-tons  -0.35003    0.20646  -1.695  0.09346 . 
## mean:Removed225-tons  -0.39790    0.29742  -1.338  0.18431   
## mean:Removed300-tons  -0.65000    0.36712  -1.771  0.08002 . 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 7.32 on 90 degrees of freedom
## Multiple R-squared:  0.2831,	Adjusted R-squared:  0.1477 
## F-statistic: 2.091 on 17 and 90 DF,  p-value: 0.01365
```



```r
# Linear Model output
lm_fraction_output(dataframe = dplyr::filter(all_divs,  measure == "Richness"))
```

```
## $free
##    Removed Adj_R2  pval  fraction
## 1   1-tons   0.05 0.241 WholeFree
## 2   5-tons  -0.08 0.666 WholeFree
## 3  10-tons  -0.10 0.952 WholeFree
## 4  30-tons  -0.08 0.702 WholeFree
## 5  60-tons  -0.05 0.504 WholeFree
## 6  90-tons   0.04 0.263 WholeFree
## 7 150-tons   0.08 0.191 WholeFree
## 8 225-tons   0.05 0.234 WholeFree
## 9 300-tons   0.06 0.223 WholeFree
## 
## $part
##    Removed Adj_R2  pval  fraction
## 1   1-tons   0.57 0.003 WholePart
## 2   5-tons   0.55 0.003 WholePart
## 3  10-tons   0.50 0.006 WholePart
## 4  30-tons   0.16 0.111 WholePart
## 5  60-tons  -0.04 0.466 WholePart
## 6  90-tons  -0.07 0.635 WholePart
## 7 150-tons   0.09 0.182 WholePart
## 8 225-tons   0.02 0.296 WholePart
## 9 300-tons   0.12 0.141 WholePart
```

```r
sig_rich_lms <- c("1-tons", "5-tons", "10-tons")

ggplot(dplyr::filter(all_divs, measure == "Richness"), 
       aes(y = frac_bacprod, x = mean, color = Removed, fill = Removed)) +
  geom_point(size = 3) + 
  xlab("Richness") +
  ylab("Bacterial Production by Fraction") +
  geom_smooth(method = "lm", data = filter(all_divs, 
                                           measure == "Richness" & fraction == "WholePart" & Removed %in% sig_rich_lms)) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(fraction~Removed, scales = "free") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

<img src="OTU_Removal_Analysis_Figs/richness-plots-1.png" style="display: block; margin: auto;" />


# Shannon Entropy

```r
### PLOT
ggplot(dplyr::filter(all_divs, measure == "Shannon_Entropy"), 
       aes(y = mean, x = Removed, color = Removed, fill = Removed)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) + geom_point(size = 3, position = position_jitter(w = 0.1)) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(.~fraction) +
  ylab("Mean Shannon_Entropy") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
```

<img src="OTU_Removal_Analysis_Figs/shannon-plots1-1.png" style="display: block; margin: auto;" />

```r
ggplot(dplyr::filter(all_divs, measure == "Shannon_Entropy"), 
       aes(y = mean, x = fraction, color = Removed, fill = Removed)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) + geom_point(size = 3, position = position_jitter(w = 0.1)) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(.~Removed) +
  ylab("Mean Shannon_Entropy") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
```

<img src="OTU_Removal_Analysis_Figs/shannon-plots1-2.png" style="display: block; margin: auto;" />

```r
summary(lm(frac_bacprod ~ mean + fraction, data = dplyr::filter(all_divs, measure == "Shannon_Entropy" & 
                                                       Removed == "1-tons")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean + fraction, data = dplyr::filter(all_divs, 
##     measure == "Shannon_Entropy" & Removed == "1-tons"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -22.152  -5.408  -1.649   4.042  36.152 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)   
## (Intercept)        -40.068     27.334  -1.466  0.15749   
## mean                10.912      5.906   1.848  0.07878 . 
## fractionWholeFree   20.301      6.280   3.232  0.00399 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 13 on 21 degrees of freedom
## Multiple R-squared:  0.3327,	Adjusted R-squared:  0.2691 
## F-statistic: 5.234 on 2 and 21 DF,  p-value: 0.01431
```

```r
summary(lm(frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, measure == "Shannon_Entropy" & fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, 
##     measure == "Shannon_Entropy" & fraction == "WholePart"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -10.869  -2.977  -0.974   1.779  13.872 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)          -52.8965     6.0044  -8.810 4.62e-14 ***
## mean                  13.6646     1.3384  10.210  < 2e-16 ***
## mean:Removed5-tons     0.4437     0.5120   0.867 0.388311    
## mean:Removed10-tons    0.7346     0.5206   1.411 0.161433    
## mean:Removed30-tons    1.4101     0.5462   2.582 0.011318 *  
## mean:Removed60-tons    1.8617     0.5671   3.283 0.001425 ** 
## mean:Removed90-tons    2.2584     0.5878   3.842 0.000217 ***
## mean:Removed150-tons   2.8748     0.6233   4.612 1.21e-05 ***
## mean:Removed225-tons   3.5678     0.6669   5.350 5.77e-07 ***
## mean:Removed300-tons   4.0108     0.6972   5.753 1.00e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.68 on 98 degrees of freedom
## Multiple R-squared:   0.53,	Adjusted R-squared:  0.4869 
## F-statistic: 12.28 on 9 and 98 DF,  p-value: 8.15e-13
```

```r
# Are the slopes of Removed different from each other?
summary(lm(frac_bacprod ~ mean*Removed, data = dplyr::filter(all_divs, measure == "Shannon_Entropy" & fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean * Removed, data = dplyr::filter(all_divs, 
##     measure == "Shannon_Entropy" & fraction == "WholePart"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -11.039  -2.590  -1.495   2.027  13.228 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)           -38.634     13.894  -2.781 0.006606 ** 
## mean                   10.599      3.008   3.523 0.000672 ***
## Removed5-tons          -6.146     20.897  -0.294 0.769335    
## Removed10-tons        -11.074     21.938  -0.505 0.614943    
## Removed30-tons        -21.134     24.202  -0.873 0.384849    
## Removed60-tons        -25.396     25.210  -1.007 0.316463    
## Removed90-tons        -28.758     26.264  -1.095 0.276457    
## Removed150-tons       -28.546     26.759  -1.067 0.288922    
## Removed225-tons       -19.226     24.830  -0.774 0.440766    
## Removed300-tons       -17.370     24.774  -0.701 0.485018    
## mean:Removed5-tons      1.706      4.606   0.370 0.712021    
## mean:Removed10-tons     3.076      4.903   0.627 0.531996    
## mean:Removed30-tons     6.110      5.607   1.090 0.278734    
## mean:Removed60-tons     7.656      5.985   1.279 0.204080    
## mean:Removed90-tons     8.970      6.377   1.407 0.162975    
## mean:Removed150-tons    9.674      6.706   1.442 0.152657    
## mean:Removed225-tons    7.984      6.374   1.253 0.213625    
## mean:Removed300-tons    7.944      6.486   1.225 0.223880    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.841 on 90 degrees of freedom
## Multiple R-squared:  0.5435,	Adjusted R-squared:  0.4573 
## F-statistic: 6.303 on 17 and 90 DF,  p-value: 1.967e-09
```



```r
# Linear Model output
lm_fraction_output(dataframe = dplyr::filter(all_divs,  measure == "Shannon_Entropy"))
```

```
## $free
##    Removed Adj_R2  pval  fraction
## 1   1-tons  -0.05 0.503 WholeFree
## 2   5-tons  -0.09 0.745 WholeFree
## 3  10-tons  -0.10 0.844 WholeFree
## 4  30-tons  -0.10 0.953 WholeFree
## 5  60-tons  -0.10 0.987 WholeFree
## 6  90-tons  -0.10 0.927 WholeFree
## 7 150-tons  -0.10 0.950 WholeFree
## 8 225-tons  -0.09 0.832 WholeFree
## 9 300-tons  -0.09 0.771 WholeFree
## 
## $part
##    Removed Adj_R2  pval  fraction
## 1   1-tons   0.52 0.005 WholePart
## 2   5-tons   0.52 0.005 WholePart
## 3  10-tons   0.53 0.005 WholePart
## 4  30-tons   0.53 0.005 WholePart
## 5  60-tons   0.53 0.005 WholePart
## 6  90-tons   0.51 0.006 WholePart
## 7 150-tons   0.47 0.008 WholePart
## 8 225-tons   0.45 0.010 WholePart
## 9 300-tons   0.42 0.013 WholePart
```

```r
ggplot(dplyr::filter(all_divs, measure == "Shannon_Entropy"), 
       aes(y = frac_bacprod, x = mean, color = Removed, fill = Removed)) +
  geom_point(size = 3) + 
  xlab("Shannon_Entropy") +
  ylab("Bacterial Production by Fraction") +
  geom_smooth(method = "lm", data = filter(all_divs, measure == "Shannon_Entropy" & fraction == "WholePart" )) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(fraction~Removed, scales = "free") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

<img src="OTU_Removal_Analysis_Figs/shannon-plots-1.png" style="display: block; margin: auto;" />


# Inverse Simpson

```r
### PLOT
ggplot(dplyr::filter(all_divs, measure == "Inverse_Simpson"), 
       aes(y = mean, x = Removed, color = Removed, fill = Removed)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) + geom_point(size = 3, position = position_jitter(w = 0.1)) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(.~fraction) +
  ylab("Mean Inverse_Simpson") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
```

<img src="OTU_Removal_Analysis_Figs/invsimps-plots1-1.png" style="display: block; margin: auto;" />

```r
ggplot(dplyr::filter(all_divs, measure == "Inverse_Simpson"), 
       aes(y = mean, x = fraction, color = Removed, fill = Removed)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) + geom_point(size = 3, position = position_jitter(w = 0.1)) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(.~Removed) +
  ylab("Mean Inverse_Simpson") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
```

<img src="OTU_Removal_Analysis_Figs/invsimps-plots1-2.png" style="display: block; margin: auto;" />

```r
summary(lm(frac_bacprod ~ mean + fraction, data = dplyr::filter(all_divs, measure == "Inverse_Simpson" & 
                                                       Removed == "1-tons")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean + fraction, data = dplyr::filter(all_divs, 
##     measure == "Inverse_Simpson" & Removed == "1-tons"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -20.233  -6.159  -0.140   2.733  37.288 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)   
## (Intercept)        -1.2879     6.4374  -0.200  0.84335   
## mean                0.3068     0.1443   2.127  0.04545 * 
## fractionWholeFree  17.8422     5.4821   3.255  0.00379 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 12.72 on 21 degrees of freedom
## Multiple R-squared:  0.3617,	Adjusted R-squared:  0.3009 
## F-statistic: 5.949 on 2 and 21 DF,  p-value: 0.008974
```

```r
summary(lm(frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, measure == "Inverse_Simpson" & fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, 
##     measure == "Inverse_Simpson" & fraction == "WholePart"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -10.3694  -2.0780  -0.7986   1.4728   9.7613 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)          -2.03446    0.96421  -2.110 0.037408 *  
## mean                  0.31228    0.03624   8.618 1.20e-13 ***
## mean:Removed5-tons    0.03872    0.04697   0.824 0.411761    
## mean:Removed10-tons   0.06390    0.04898   1.304 0.195140    
## mean:Removed30-tons   0.11820    0.05378   2.198 0.030323 *  
## mean:Removed60-tons   0.15147    0.05701   2.657 0.009211 ** 
## mean:Removed90-tons   0.18186    0.06018   3.022 0.003206 ** 
## mean:Removed150-tons  0.22972    0.06528   3.519 0.000658 ***
## mean:Removed225-tons  0.27327    0.07001   3.903 0.000174 ***
## mean:Removed300-tons  0.30889    0.07410   4.169 6.62e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 4.748 on 98 degrees of freedom
## Multiple R-squared:  0.6716,	Adjusted R-squared:  0.6415 
## F-statistic: 22.27 on 9 and 98 DF,  p-value: < 2.2e-16
```

```r
# Are the slopes of Removed different from each other?
summary(lm(frac_bacprod ~ mean*Removed, data = dplyr::filter(all_divs, measure == "Inverse_Simpson" & fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean * Removed, data = dplyr::filter(all_divs, 
##     measure == "Inverse_Simpson" & fraction == "WholePart"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -10.3907  -2.1335  -0.5506   1.7252   9.4510 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)          -0.22434    2.58861  -0.087   0.9311    
## mean                  0.27781    0.05901   4.708 9.02e-06 ***
## Removed5-tons        -0.87070    3.76140  -0.231   0.8175    
## Removed10-tons       -1.44156    3.83739  -0.376   0.7081    
## Removed30-tons       -2.35080    3.97971  -0.591   0.5562    
## Removed60-tons       -2.62043    4.04640  -0.648   0.5189    
## Removed90-tons       -2.79312    4.11453  -0.679   0.4990    
## Removed150-tons      -2.83775    4.17493  -0.680   0.4984    
## Removed225-tons      -2.38590    4.15132  -0.575   0.5669    
## Removed300-tons      -2.38496    4.19628  -0.568   0.5712    
## mean:Removed5-tons    0.05273    0.09128   0.578   0.5649    
## mean:Removed10-tons   0.08966    0.09736   0.921   0.3596    
## mean:Removed30-tons   0.16763    0.11166   1.501   0.1368    
## mean:Removed60-tons   0.21038    0.12084   1.741   0.0851 .  
## mean:Removed90-tons   0.24831    0.13027   1.906   0.0598 .  
## mean:Removed150-tons  0.30128    0.14388   2.094   0.0391 *  
## mean:Removed225-tons  0.33024    0.15296   2.159   0.0335 *  
## mean:Removed300-tons  0.36740    0.16401   2.240   0.0275 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 4.928 on 90 degrees of freedom
## Multiple R-squared:  0.6752,	Adjusted R-squared:  0.6138 
## F-statistic:    11 on 17 and 90 DF,  p-value: 2.107e-15
```


```r
# Linear Model output
lm_fraction_output(dataframe = dplyr::filter(all_divs,  measure == "Inverse_Simpson"))
```

```
## $free
##    Removed Adj_R2  pval  fraction
## 1   1-tons  -0.02 0.392 WholeFree
## 2   5-tons  -0.04 0.457 WholeFree
## 3  10-tons  -0.04 0.481 WholeFree
## 4  30-tons  -0.05 0.499 WholeFree
## 5  60-tons  -0.05 0.501 WholeFree
## 6  90-tons  -0.05 0.514 WholeFree
## 7 150-tons  -0.05 0.501 WholeFree
## 8 225-tons  -0.03 0.418 WholeFree
## 9 300-tons  -0.02 0.388 WholeFree
## 
## $part
##    Removed Adj_R2  pval  fraction
## 1   1-tons   0.69 0.000 WholePart
## 2   5-tons   0.70 0.000 WholePart
## 3  10-tons   0.70 0.000 WholePart
## 4  30-tons   0.69 0.001 WholePart
## 5  60-tons   0.67 0.001 WholePart
## 6  90-tons   0.63 0.001 WholePart
## 7 150-tons   0.60 0.002 WholePart
## 8 225-tons   0.56 0.003 WholePart
## 9 300-tons   0.54 0.004 WholePart
```

```r
ggplot(dplyr::filter(all_divs, measure == "Inverse_Simpson"), 
       aes(y = frac_bacprod, x = mean, color = Removed, fill = Removed)) +
  geom_point(size = 3) + 
  xlab("Inverse_Simpson") +
  ylab("Bacterial Production by Fraction") +
  geom_smooth(method = "lm", data = filter(all_divs, measure == "Inverse_Simpson" & fraction == "WholePart" )) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(fraction~Removed, scales = "free") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

<img src="OTU_Removal_Analysis_Figs/inverse-simpson-plots-1.png" style="display: block; margin: auto;" />


# Simpson's Evenness

```r
### PLOT
ggplot(dplyr::filter(all_divs, measure == "Simpsons_Evenness"), 
       aes(y = mean, x = Removed, color = Removed, fill = Removed)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) + geom_point(size = 3, position = position_jitter(w = 0.1)) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(.~fraction) +
  ylab("Mean Simpsons_Evenness") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
```

<img src="OTU_Removal_Analysis_Figs/simpseven-plots1-1.png" style="display: block; margin: auto;" />

```r
ggplot(dplyr::filter(all_divs, measure == "Simpsons_Evenness"), 
       aes(y = mean, x = fraction, color = Removed, fill = Removed)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) + geom_point(size = 3, position = position_jitter(w = 0.1)) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(.~Removed) +
  ylab("Mean Simpsons_Evenness") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
```

<img src="OTU_Removal_Analysis_Figs/simpseven-plots1-2.png" style="display: block; margin: auto;" />

```r
summary(lm(frac_bacprod ~ mean + fraction, data = dplyr::filter(all_divs, measure == "Simpsons_Evenness" & 
                                                       Removed == "1-tons")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean + fraction, data = dplyr::filter(all_divs, 
##     measure == "Simpsons_Evenness" & Removed == "1-tons"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -19.541  -7.275  -1.774   0.993  40.759 
## 
## Coefficients:
##                   Estimate Std. Error t value Pr(>|t|)  
## (Intercept)         -0.263      8.889  -0.030   0.9767  
## mean               184.151    143.925   1.279   0.2147  
## fractionWholeFree   12.045      5.742   2.098   0.0482 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 13.5 on 21 degrees of freedom
## Multiple R-squared:  0.2803,	Adjusted R-squared:  0.2117 
## F-statistic: 4.089 on 2 and 21 DF,  p-value: 0.03164
```

```r
summary(lm(frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, measure == "Simpsons_Evenness" & fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, 
##     measure == "Simpsons_Evenness" & fraction == "WholePart"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -11.0505  -2.5797  -0.7464   1.9396  13.1164 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)            -2.856      1.154  -2.476 0.015020 *  
## mean                  233.840     30.756   7.603 1.76e-11 ***
## mean:Removed5-tons    -35.054     32.836  -1.068 0.288350    
## mean:Removed10-tons   -56.983     31.573  -1.805 0.074180 .  
## mean:Removed30-tons   -98.211     29.871  -3.288 0.001402 ** 
## mean:Removed60-tons  -118.576     29.378  -4.036 0.000108 ***
## mean:Removed90-tons  -132.604     29.170  -4.546 1.56e-05 ***
## mean:Removed150-tons -146.394     29.063  -5.037 2.16e-06 ***
## mean:Removed225-tons -157.873     29.044  -5.436 4.00e-07 ***
## mean:Removed300-tons -164.232     29.062  -5.651 1.57e-07 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.186 on 98 degrees of freedom
## Multiple R-squared:  0.6082,	Adjusted R-squared:  0.5722 
## F-statistic:  16.9 on 9 and 98 DF,  p-value: < 2.2e-16
```

```r
# Are the slopes of Removed different from each other?
summary(lm(frac_bacprod ~ mean*Removed, data = dplyr::filter(all_divs, measure == "Simpsons_Evenness" & fraction == "WholePart")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean * Removed, data = dplyr::filter(all_divs, 
##     measure == "Simpsons_Evenness" & fraction == "WholePart"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -11.0044  -2.1792  -0.6513   1.8458  12.7042 
## 
## Coefficients:
##                       Estimate Std. Error t value Pr(>|t|)    
## (Intercept)            -4.0388     4.1759  -0.967 0.336048    
## mean                  252.1779    69.7898   3.613 0.000497 ***
## Removed5-tons           0.1749     5.7318   0.031 0.975725    
## Removed10-tons          0.4974     5.6325   0.088 0.929827    
## Removed30-tons          0.9554     5.4741   0.175 0.861836    
## Removed60-tons          1.2497     5.4425   0.230 0.818910    
## Removed90-tons          1.5366     5.4093   0.284 0.777010    
## Removed150-tons         1.7644     5.3936   0.327 0.744335    
## Removed225-tons         1.8233     5.4107   0.337 0.736919    
## Removed300-tons         1.9021     5.4256   0.351 0.726727    
## mean:Removed5-tons    -40.3843    88.9879  -0.454 0.651052    
## mean:Removed10-tons   -67.5423    84.1917  -0.802 0.424525    
## mean:Removed30-tons  -114.6162    77.4198  -1.480 0.142245    
## mean:Removed60-tons  -137.3961    75.2410  -1.826 0.071154 .  
## mean:Removed90-tons  -153.1747    73.9150  -2.072 0.041097 *  
## mean:Removed150-tons -167.9025    72.8577  -2.305 0.023494 *  
## mean:Removed225-tons -179.2587    72.1670  -2.484 0.014845 *  
## mean:Removed300-tons -185.7204    71.8293  -2.586 0.011328 *  
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.403 on 90 degrees of freedom
## Multiple R-squared:  0.6095,	Adjusted R-squared:  0.5357 
## F-statistic: 8.263 on 17 and 90 DF,  p-value: 3.986e-12
```


```r
# Linear Model output
lm_fraction_output(dataframe = dplyr::filter(all_divs,  measure == "Simpsons_Evenness"))
```

```
## $free
##    Removed Adj_R2  pval  fraction
## 1   1-tons  -0.10 0.912 WholeFree
## 2   5-tons  -0.05 0.520 WholeFree
## 3  10-tons  -0.02 0.393 WholeFree
## 4  30-tons   0.03 0.278 WholeFree
## 5  60-tons   0.05 0.233 WholeFree
## 6  90-tons   0.08 0.186 WholeFree
## 7 150-tons   0.08 0.192 WholeFree
## 8 225-tons   0.08 0.196 WholeFree
## 9 300-tons   0.07 0.200 WholeFree
## 
## $part
##    Removed Adj_R2  pval  fraction
## 1   1-tons   0.46 0.009 WholePart
## 2   5-tons   0.53 0.004 WholePart
## 3  10-tons   0.56 0.003 WholePart
## 4  30-tons   0.62 0.001 WholePart
## 5  60-tons   0.62 0.002 WholePart
## 6  90-tons   0.61 0.002 WholePart
## 7 150-tons   0.60 0.002 WholePart
## 8 225-tons   0.58 0.003 WholePart
## 9 300-tons   0.56 0.003 WholePart
```

```r
ggplot(dplyr::filter(all_divs, measure == "Simpsons_Evenness"), 
       aes(y = frac_bacprod, x = mean, color = Removed, fill = Removed)) +
  geom_point(size = 3) + 
  xlab("Simpsons_Evenness") +
  ylab("Bacterial Production by Fraction") +
  geom_smooth(method = "lm", data = filter(all_divs, measure == "Simpsons_Evenness" & fraction == "WholePart" )) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(fraction~Removed, scales = "free") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

<img src="OTU_Removal_Analysis_Figs/simps-evenness-plots-1.png" style="display: block; margin: auto;" />











## Subset Diversity Data 

```r
######################################################### RICHNESS
# Subset only richness data 
ML_rm1_rich_stats <- filter(otu_alphadiv, measure == "Richness" & 
                              fraction %in% c("WholePart", "WholeFree") & year == "2015")

######################################################### SHANNON ENTROPY
# Subset only Shannon_Entropy data 
ML_otu_shannon_stats <- filter(otu_alphadiv, measure == "Shannon_Entropy" & 
                                 fraction %in% c("WholePart", "WholeFree") & year == "2015")

######################################################### INVERSE SIMPSON
# Subset only Inverse_Simpson data 
ML_otu_invsimps_stats <- filter(otu_alphadiv, measure == "Inverse_Simpson" & 
                                  fraction %in% c("WholePart", "WholeFree") & year == "2015")

######################################################### SIMPSON'S EVENNESS
# Subset only Simpsons_Evenness data 
ML_otu_simpseven_stats <- filter(otu_alphadiv, measure == "Simpsons_Evenness" & 
                                   fraction %in% c("WholePart", "WholeFree") & year == "2015")
```



# Diversity vs Fraction Production 

```r
######################################################### RICHNESS
# Free-Living Richness vs fractional production 
freeprod_ML_otu_rich <- lm(frac_bacprod ~ mean, data = filter(ML_otu_rich_stats, fraction == "WholeFree"))
summary(freeprod_ML_otu_rich) 

# Particle-Associated Richness vs fractional production 
partprod_MLotu_rich <- lm(frac_bacprod ~ mean, data = filter(ML_otu_rich_stats, fraction == "WholePart"))
summary(partprod_MLotu_rich) 

# Both fractions together
summary(lm(frac_bacprod ~ mean, data = ML_otu_rich_stats))
 
# Plot 
ggplot(ML_otu_rich_stats, aes(x=mean, y=frac_bacprod, color = fraction)) + 
  geom_point(size = 3.5) + 
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), alpha = 0.7) + # X-axis errorbars
  # Y-axis errorbars
  geom_errorbar(aes(ymin = frac_bacprod - SD_frac_bacprod, ymax = frac_bacprod + SD_frac_bacprod),  alpha = 0.5) + 
  scale_color_manual(values = c("firebrick3","cornflowerblue"), 
                     limits = c("WholePart", "WholeFree"),
                     breaks=c("WholePart", "WholeFree"),
                     labels=c("Particle", "Free")) + 
  ylab("Heterotrophic Production (μgC/L/hr)") + xlab("Observed Richness") +
  geom_smooth(data=subset(ML_otu_rich_stats, fraction == "WholePart"), method='lm') + 
  #scale_x_continuous(limits = c(180, 810), breaks = c(200, 400, 600, 800)) + 
  theme(legend.position=c(0.15,0.9),        
        legend.title=element_blank()) 
  #annotate("text", x = 250, y=55, color = "cornflowerblue", fontface = "bold",
  #         label = paste("R2 =", round(summary(freeprod_ML_otu_rich)$adj.r.squared, digits = 2), "\n", 
  #                       "p =", round(unname(summary(freeprod_ML_otu_rich)$coefficients[,4][2]), digits = 2))) + 
  #annotate("text", x = 650, y=3, color = "firebrick3", fontface = "bold",
  #         label = paste("R2 =", round(summary(partprod_MLotu_rich)$adj.r.squared, digits = 2), "\n", 
  #                       "p =", round(unname(summary(partprod_MLotu_rich)$coefficients[,4][2]), digits = 4)));


######################################################### SHANNON ENTROPY
# Free-Living Shannon vs fractional production 
freeprod_ML_otu_shannon <- lm(frac_bacprod ~ mean, data = filter(ML_otu_shannon_stats, fraction == "WholeFree"))
summary(freeprod_ML_otu_shannon)

# Particle-Associated Shannon vs fractional production 
partprod_MLotu_shannon <- lm(frac_bacprod ~ mean, data = filter(ML_otu_shannon_stats, fraction == "WholePart"))
summary(partprod_MLotu_shannon)

# Both fractions together
summary(lm(frac_bacprod ~ mean, data = ML_otu_shannon_stats))

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

# Particle-Associated Inverse Simpson vs fractional production 
partprod_MLotu_invsimps <- lm(frac_bacprod ~ mean, data = filter(ML_otu_invsimps_stats, fraction == "WholePart"))
summary(partprod_MLotu_invsimps)

# Both fractions together
summary(lm(frac_bacprod ~ mean, data = ML_otu_invsimps_stats))

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

# Particle-Associated Simpson's Evenness vs fractional production 
partprod_MLotu_simpseven <- lm(frac_bacprod ~ mean, data = filter(ML_otu_simpseven_stats, fraction == "WholePart"))
summary(partprod_MLotu_simpseven)

# Both fractions together
summary(lm(frac_bacprod ~ mean, data = ML_otu_simpseven_stats))

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

otu_vegan <- plot_grid(otu_rich_vegan, otu_simpseven_vegan, otu_shannon_vegan, otu_invsimps_vegan,
                       labels = c("A", "B", "C", "D"), 
                       align = "h", nrow = 2, ncol = 2)
otu_vegan
```





