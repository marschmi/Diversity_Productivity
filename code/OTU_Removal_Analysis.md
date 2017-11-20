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
library(forcats)
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
# If taxa with 20 counts are removed 
notree_musk_surface_pruned_rm20 <- prune_taxa(taxa_sums(notree_musk_surface_pruned) > 20, notree_musk_surface_pruned) 
notree_musk_surface_pruned_rm20
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 694 taxa and 24 samples ]
## sample_data() Sample Data:       [ 24 samples by 13 sample variables ]
## tax_table()   Taxonomy Table:    [ 694 taxa by 8 taxonomic ranks ]
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


################## Remove 20-tons 
alpha_rm20 <- calc_alpha_diversity(physeq = notree_musk_surface_pruned_rm20)
min(sample_sums(notree_musk_surface_pruned_rm20)) - 1
```

```
## [1] 6067
```

```r
otu_alphadiv_rm20 <- calc_mean_alphadiv(physeq = notree_musk_surface_pruned_rm20,
                           richness_df = alpha_rm20$Richness, 
                           evenness_df = alpha_rm20$Inverse_Simpson, 
                           shannon_df = alpha_rm20$Shannon) %>%
    mutate(fraction = factor(fraction, levels = c("WholePart", "Particle", "WholeFree", "Free")),
         lakesite = factor(lakesite,  levels = c("MOT", "MDP", "MBR", "MIN")),
         measure = factor(measure, levels = c("Richness", "Simpsons_Evenness", "Shannon_Entropy", "Inverse_Simpson")),
         Removed = "20-tons")  




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


# Minimum Sequences Plot

```r
min_seqs <- c(min(sample_sums(notree_musk_surface_pruned_rm1)) - 1, min(sample_sums(notree_musk_surface_pruned_rm5)) - 1, 
  min(sample_sums(notree_musk_surface_pruned_rm10)) - 1, min(sample_sums(notree_musk_surface_pruned_rm20)) - 1,
  min(sample_sums(notree_musk_surface_pruned_rm30)) - 1,
  min(sample_sums(notree_musk_surface_pruned_rm60)) - 1, min(sample_sums(notree_musk_surface_pruned_rm90)) - 1,
  min(sample_sums(notree_musk_surface_pruned_rm150)) - 1, min(sample_sums(notree_musk_surface_pruned_rm225)) - 1,
  min(sample_sums(notree_musk_surface_pruned_rm300)) - 1)


num_otus <- c(ncol(otu_table(notree_musk_surface_pruned_rm1)), ncol(otu_table(notree_musk_surface_pruned_rm5)), ncol(otu_table(notree_musk_surface_pruned_rm10)), 
              ncol(otu_table(notree_musk_surface_pruned_rm20)),
              ncol(otu_table(notree_musk_surface_pruned_rm30)), ncol(otu_table(notree_musk_surface_pruned_rm60)), ncol(otu_table(notree_musk_surface_pruned_rm90)),
              ncol(otu_table(notree_musk_surface_pruned_rm150)), ncol(otu_table(notree_musk_surface_pruned_rm225)), ncol(otu_table(notree_musk_surface_pruned_rm300)))

Removed <- c("1-tons","5-tons", "10-tons", "20-tons", "30-tons", "60-tons", "90-tons", "150-tons", "225-tons","300-tons")


statz <- data.frame(cbind(as.numeric(min_seqs), as.numeric(num_otus), Removed)) %>%
         mutate(Removed = factor(Removed, levels = c("1-tons","5-tons", "10-tons", "20-tons", "30-tons", "60-tons", 
                                              "90-tons", "150-tons", "225-tons","300-tons"))) %>%
  mutate(num_otus = as.numeric(num_otus))

p1 <- ggplot(statz, aes(x = Removed, y = min_seqs, fill = Removed)) +
  geom_bar(stat = "identity") + ylab("Minimum # of Sequences") +
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
        legend.position = c(0.8, 0.65), legend.title = element_blank())

plot_grid(p1, p2, align = "v", labels = c("A", "B"), nrow = 2, ncol =1)
```

<img src="OTU_Removal_Analysis_Figs/min-seqs-plots-1.png" style="display: block; margin: auto;" />






# Combine all data 

```r
# Combine all div metrics 
all_divs <- bind_rows(otu_alphadiv_rm1, otu_alphadiv_rm5, otu_alphadiv_rm10, otu_alphadiv_rm20, otu_alphadiv_rm30, 
                      otu_alphadiv_rm60, otu_alphadiv_rm90, otu_alphadiv_rm150, 
                      otu_alphadiv_rm225, otu_alphadiv_rm300) %>%
  dplyr::filter(fraction %in% c("WholePart", "WholeFree") & year == "2015") %>% 
  mutate(fraction = fct_recode(fraction, "Particle" = "WholePart", "Free" = "WholeFree")) %>%
  mutate(Removed = factor(Removed, levels = c("1-tons","5-tons", "10-tons", "20-tons", "30-tons", "60-tons", 
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
##               Estimate Std. Error t value Pr(>|t|)   
## (Intercept)  -13.35169   10.71466  -1.246   0.2264   
## mean           0.03843    0.01663   2.311   0.0311 * 
## fractionFree  23.18701    6.44851   3.596   0.0017 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 12.52 on 21 degrees of freedom
## Multiple R-squared:  0.3815,	Adjusted R-squared:  0.3226 
## F-statistic: 6.476 on 2 and 21 DF,  p-value: 0.006446
```

```r
# Are the linear models for different OTU removals significant? 
summary(lm(frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, measure == "Richness" & fraction == "Particle")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, 
##     measure == "Richness" & fraction == "Particle"))
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -10.089  -5.440  -1.310   3.722  22.488 
## 
## Coefficients:
##                        Estimate Std. Error t value Pr(>|t|)    
## (Intercept)          -11.820593   4.947751  -2.389 0.018608 *  
## mean                   0.035581   0.008183   4.348  3.1e-05 ***
## mean:Removed5-tons     0.010441   0.006144   1.699 0.092104 .  
## mean:Removed10-tons    0.018249   0.007613   2.397 0.018231 *  
## mean:Removed20-tons    0.029796   0.010190   2.924 0.004202 ** 
## mean:Removed30-tons    0.039610   0.012515   3.165 0.002011 ** 
## mean:Removed60-tons    0.056972   0.016733   3.405 0.000928 ***
## mean:Removed90-tons    0.073401   0.020903   3.511 0.000649 ***
## mean:Removed150-tons   0.100057   0.027640   3.620 0.000448 ***
## mean:Removed225-tons   0.136096   0.036460   3.733 0.000303 ***
## mean:Removed300-tons   0.162378   0.043085   3.769 0.000267 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 7.591 on 109 degrees of freedom
## Multiple R-squared:  0.1597,	Adjusted R-squared:  0.08265 
## F-statistic: 2.072 on 10 and 109 DF,  p-value: 0.03274
```



```r
# Linear Model output
lm_fraction_output(dataframe = dplyr::filter(all_divs,  measure == "Richness"))
```

```
## $free
##     Removed Adj_R2  pval fraction
## 1    1-tons   0.05 0.241     Free
## 2    5-tons  -0.08 0.666     Free
## 3   10-tons  -0.10 0.952     Free
## 4   20-tons  -0.09 0.779     Free
## 5   30-tons  -0.08 0.702     Free
## 6   60-tons  -0.05 0.504     Free
## 7   90-tons   0.04 0.263     Free
## 8  150-tons   0.08 0.191     Free
## 9  225-tons   0.05 0.234     Free
## 10 300-tons   0.06 0.223     Free
## 
## $part
##     Removed Adj_R2  pval fraction
## 1    1-tons   0.57 0.003 Particle
## 2    5-tons   0.55 0.003 Particle
## 3   10-tons   0.50 0.006 Particle
## 4   20-tons   0.35 0.025 Particle
## 5   30-tons   0.16 0.111 Particle
## 6   60-tons  -0.04 0.466 Particle
## 7   90-tons  -0.07 0.635 Particle
## 8  150-tons   0.09 0.182 Particle
## 9  225-tons   0.02 0.296 Particle
## 10 300-tons   0.12 0.141 Particle
```

```r
sig_rich_lms <- c("1-tons", "5-tons", "10-tons")

ggplot(dplyr::filter(all_divs, measure == "Richness"), 
       aes(y = frac_bacprod, x = mean, color = Removed, fill = Removed)) +
  geom_point(size = 3) + 
  xlab("Richness") +
  ylab("Bacterial Production by Fraction") +
  geom_smooth(method = "lm", data = filter(all_divs, 
                                           measure == "Richness" & fraction == "Particle" & Removed %in% sig_rich_lms)) + 
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
##              Estimate Std. Error t value Pr(>|t|)   
## (Intercept)   -40.068     27.334  -1.466  0.15749   
## mean           10.912      5.906   1.848  0.07878 . 
## fractionFree   20.301      6.280   3.232  0.00399 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 13 on 21 degrees of freedom
## Multiple R-squared:  0.3327,	Adjusted R-squared:  0.2691 
## F-statistic: 5.234 on 2 and 21 DF,  p-value: 0.01431
```

```r
# Are the linear models for different OTU removals significant? 
summary(lm(frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, measure == "Shannon_Entropy" & fraction == "Particle")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, 
##     measure == "Shannon_Entropy" & fraction == "Particle"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -10.8903  -2.9868  -0.9647   1.8390  13.8331 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)          -53.2811     5.6778  -9.384 1.09e-15 ***
## mean                  13.7473     1.2703  10.822  < 2e-16 ***
## mean:Removed5-tons     0.4465     0.5092   0.877 0.382500    
## mean:Removed10-tons    0.7392     0.5175   1.429 0.155997    
## mean:Removed20-tons    1.1273     0.5305   2.125 0.035839 *  
## mean:Removed30-tons    1.4189     0.5417   2.619 0.010065 *  
## mean:Removed60-tons    1.8733     0.5614   3.337 0.001159 ** 
## mean:Removed90-tons    2.2725     0.5808   3.913 0.000159 ***
## mean:Removed150-tons   2.8927     0.6141   4.711 7.32e-06 ***
## mean:Removed225-tons   3.5898     0.6549   5.482 2.75e-07 ***
## mean:Removed300-tons   4.0355     0.6833   5.906 4.04e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.65 on 109 degrees of freedom
## Multiple R-squared:  0.5344,	Adjusted R-squared:  0.4917 
## F-statistic: 12.51 on 10 and 109 DF,  p-value: 3.071e-14
```



```r
# Linear Model output
lm_fraction_output(dataframe = dplyr::filter(all_divs,  measure == "Shannon_Entropy"))
```

```
## $free
##     Removed Adj_R2  pval fraction
## 1    1-tons  -0.05 0.503     Free
## 2    5-tons  -0.09 0.745     Free
## 3   10-tons  -0.10 0.844     Free
## 4   20-tons  -0.10 0.900     Free
## 5   30-tons  -0.10 0.953     Free
## 6   60-tons  -0.10 0.987     Free
## 7   90-tons  -0.10 0.927     Free
## 8  150-tons  -0.10 0.950     Free
## 9  225-tons  -0.09 0.832     Free
## 10 300-tons  -0.09 0.771     Free
## 
## $part
##     Removed Adj_R2  pval fraction
## 1    1-tons   0.52 0.005 Particle
## 2    5-tons   0.52 0.005 Particle
## 3   10-tons   0.53 0.005 Particle
## 4   20-tons   0.53 0.004 Particle
## 5   30-tons   0.53 0.005 Particle
## 6   60-tons   0.53 0.005 Particle
## 7   90-tons   0.51 0.006 Particle
## 8  150-tons   0.47 0.008 Particle
## 9  225-tons   0.45 0.010 Particle
## 10 300-tons   0.42 0.013 Particle
```

```r
ggplot(dplyr::filter(all_divs, measure == "Shannon_Entropy"), 
       aes(y = frac_bacprod, x = mean, color = Removed, fill = Removed)) +
  geom_point(size = 3) + 
  xlab("Shannon_Entropy") +
  ylab("Bacterial Production by Fraction") +
  geom_smooth(method = "lm", data = filter(all_divs, measure == "Shannon_Entropy" & fraction == "Particle" )) + 
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
##              Estimate Std. Error t value Pr(>|t|)   
## (Intercept)   -1.2879     6.4374  -0.200  0.84335   
## mean           0.3068     0.1443   2.127  0.04545 * 
## fractionFree  17.8422     5.4821   3.255  0.00379 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 12.72 on 21 degrees of freedom
## Multiple R-squared:  0.3617,	Adjusted R-squared:  0.3009 
## F-statistic: 5.949 on 2 and 21 DF,  p-value: 0.008974
```

```r
# Are the linear models for different OTU removals significant? 
summary(lm(frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, measure == "Inverse_Simpson" & fraction == "Particle")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, 
##     measure == "Inverse_Simpson" & fraction == "Particle"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -10.3703  -2.0544  -0.7842   1.4876   9.7477 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)          -2.05965    0.90498  -2.276 0.024806 *  
## mean                  0.31276    0.03543   8.826 2.01e-14 ***
## mean:Removed5-tons    0.03879    0.04654   0.833 0.406450    
## mean:Removed10-tons   0.06401    0.04853   1.319 0.189876    
## mean:Removed20-tons   0.09628    0.05124   1.879 0.062921 .  
## mean:Removed30-tons   0.11841    0.05323   2.224 0.028180 *  
## mean:Removed60-tons   0.15175    0.05640   2.691 0.008254 ** 
## mean:Removed90-tons   0.18220    0.05950   3.062 0.002768 ** 
## mean:Removed150-tons  0.23015    0.06448   3.569 0.000534 ***
## mean:Removed225-tons  0.27377    0.06911   3.961 0.000133 ***
## mean:Removed300-tons  0.30947    0.07310   4.233 4.82e-05 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 4.705 on 109 degrees of freedom
## Multiple R-squared:  0.6772,	Adjusted R-squared:  0.6475 
## F-statistic: 22.86 on 10 and 109 DF,  p-value: < 2.2e-16
```


```r
# Linear Model output
lm_fraction_output(dataframe = dplyr::filter(all_divs,  measure == "Inverse_Simpson"))
```

```
## $free
##     Removed Adj_R2  pval fraction
## 1    1-tons  -0.02 0.392     Free
## 2    5-tons  -0.04 0.457     Free
## 3   10-tons  -0.04 0.481     Free
## 4   20-tons  -0.04 0.484     Free
## 5   30-tons  -0.05 0.499     Free
## 6   60-tons  -0.05 0.501     Free
## 7   90-tons  -0.05 0.514     Free
## 8  150-tons  -0.05 0.501     Free
## 9  225-tons  -0.03 0.418     Free
## 10 300-tons  -0.02 0.388     Free
## 
## $part
##     Removed Adj_R2  pval fraction
## 1    1-tons   0.69 0.000 Particle
## 2    5-tons   0.70 0.000 Particle
## 3   10-tons   0.70 0.000 Particle
## 4   20-tons   0.70 0.000 Particle
## 5   30-tons   0.69 0.001 Particle
## 6   60-tons   0.67 0.001 Particle
## 7   90-tons   0.63 0.001 Particle
## 8  150-tons   0.60 0.002 Particle
## 9  225-tons   0.56 0.003 Particle
## 10 300-tons   0.54 0.004 Particle
```

```r
ggplot(dplyr::filter(all_divs, measure == "Inverse_Simpson"), 
       aes(y = frac_bacprod, x = mean, color = Removed, fill = Removed)) +
  geom_point(size = 3) + 
  xlab("Inverse_Simpson") +
  ylab("Bacterial Production by Fraction") +
  geom_smooth(method = "lm", data = filter(all_divs, measure == "Inverse_Simpson" & fraction == "Particle" )) + 
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
##              Estimate Std. Error t value Pr(>|t|)  
## (Intercept)    -0.263      8.889  -0.030   0.9767  
## mean          184.151    143.925   1.279   0.2147  
## fractionFree   12.045      5.742   2.098   0.0482 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 13.5 on 21 degrees of freedom
## Multiple R-squared:  0.2803,	Adjusted R-squared:  0.2117 
## F-statistic: 4.089 on 2 and 21 DF,  p-value: 0.03164
```

```r
# Are the linear models for different OTU removals significant? 
summary(lm(frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, measure == "Simpsons_Evenness" & fraction == "Particle")))
```

```
## 
## Call:
## lm(formula = frac_bacprod ~ mean/Removed, data = dplyr::filter(all_divs, 
##     measure == "Simpsons_Evenness" & fraction == "Particle"))
## 
## Residuals:
##      Min       1Q   Median       3Q      Max 
## -11.0520  -2.5509  -0.7527   1.9354  13.1082 
## 
## Coefficients:
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)            -2.879      1.090  -2.642  0.00947 ** 
## mean                  234.204     30.116   7.777 4.49e-12 ***
## mean:Removed5-tons    -35.115     32.700  -1.074  0.28526    
## mean:Removed10-tons   -57.081     31.419  -1.817  0.07200 .  
## mean:Removed20-tons   -82.839     30.217  -2.741  0.00715 ** 
## mean:Removed30-tons   -98.376     29.649  -3.318  0.00123 ** 
## mean:Removed60-tons  -118.771     29.112  -4.080 8.61e-05 ***
## mean:Removed90-tons  -132.820     28.867  -4.601 1.14e-05 ***
## mean:Removed150-tons -146.630     28.722  -5.105 1.41e-06 ***
## mean:Removed225-tons -158.125     28.669  -5.515 2.36e-07 ***
## mean:Removed300-tons -164.493     28.669  -5.738 8.72e-08 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 5.167 on 109 degrees of freedom
## Multiple R-squared:  0.6107,	Adjusted R-squared:  0.575 
## F-statistic:  17.1 on 10 and 109 DF,  p-value: < 2.2e-16
```


```r
# Linear Model output
lm_fraction_output(dataframe = dplyr::filter(all_divs,  measure == "Simpsons_Evenness"))
```

```
## $free
##     Removed Adj_R2  pval fraction
## 1    1-tons  -0.10 0.912     Free
## 2    5-tons  -0.05 0.520     Free
## 3   10-tons  -0.02 0.393     Free
## 4   20-tons   0.02 0.298     Free
## 5   30-tons   0.03 0.278     Free
## 6   60-tons   0.05 0.233     Free
## 7   90-tons   0.08 0.186     Free
## 8  150-tons   0.08 0.192     Free
## 9  225-tons   0.08 0.196     Free
## 10 300-tons   0.07 0.200     Free
## 
## $part
##     Removed Adj_R2  pval fraction
## 1    1-tons   0.46 0.009 Particle
## 2    5-tons   0.53 0.004 Particle
## 3   10-tons   0.56 0.003 Particle
## 4   20-tons   0.60 0.002 Particle
## 5   30-tons   0.62 0.001 Particle
## 6   60-tons   0.62 0.002 Particle
## 7   90-tons   0.61 0.002 Particle
## 8  150-tons   0.60 0.002 Particle
## 9  225-tons   0.58 0.003 Particle
## 10 300-tons   0.56 0.003 Particle
```

```r
ggplot(dplyr::filter(all_divs, measure == "Simpsons_Evenness"), 
       aes(y = frac_bacprod, x = mean, color = Removed, fill = Removed)) +
  geom_point(size = 3) + 
  xlab("Simpsons_Evenness") +
  ylab("Bacterial Production by Fraction") +
  geom_smooth(method = "lm", data = filter(all_divs, measure == "Simpsons_Evenness" & fraction == "Particle" )) + 
  scale_color_manual(values = tons_colors) +
  scale_fill_manual(values = tons_colors) +  
  facet_grid(fraction~Removed, scales = "free") +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
```

<img src="OTU_Removal_Analysis_Figs/simps-evenness-plots-1.png" style="display: block; margin: auto;" />
