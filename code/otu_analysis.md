# Diversity-Productivity Relationships in Muskegon Lake
Marian L. Schmidt  
January 2017  



# Load Libraries


## Load Mothur OTU Data 




# Sample Sequencing Read Counts
<img src="Figures/cached/seq-read-count-1.png" style="display: block; margin: auto;" /><img src="Figures/cached/seq-read-count-2.png" style="display: block; margin: auto;" /><img src="Figures/cached/seq-read-count-3.png" style="display: block; margin: auto;" /><img src="Figures/cached/seq-read-count-4.png" style="display: block; margin: auto;" />


# Scale Reads for Between Sample Diversity Analysis 

```
## [1] 1562
```

```
## [1] 1562
```

```
## [1] 1562
```

```
## [1] 6665
```

```
## [1] 7887
```

```
## [1] 12162
```

<img src="Figures/cached/scaled_reads_seq_depth-1.png" style="display: block; margin: auto;" /><img src="Figures/cached/scaled_reads_seq_depth-2.png" style="display: block; margin: auto;" />





##  How linked are the diversities of the fractions?
<img src="Figures/cached/free_vs_wholefree_diversity-1.png" style="display: block; margin: auto;" />

##  How linked are the diversities of the fractions?
<img src="Figures/cached/part_vs_wholepart_diversity-1.png" style="display: block; margin: auto;" />




# Is there a relationship between diversity and productivity?
## For total production?

### D1 Diversity with Total Production Analysis 
<img src="Figures/cached/D1_totalproduction_vs_diversity-1.png" style="display: block; margin: auto;" />

### D2 Diversity with Total Production Analysis 
<img src="Figures/cached/D2_totalproduction_vs_diversity-1.png" style="display: block; margin: auto;" /><img src="Figures/cached/D2_totalproduction_vs_diversity-2.png" style="display: block; margin: auto;" />

# Is there a relationship between HNA cells per uL and Total Production?
<img src="Figures/cached/totalproduction_vs_HNA-LNA-1.png" style="display: block; margin: auto;" /><img src="Figures/cached/totalproduction_vs_HNA-LNA-2.png" style="display: block; margin: auto;" />

```
## [1] 0.4911129
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  free_meta_data$D2 and free_meta_data$HNA_percent
## t = 4.6491, df = 68, p-value = 1.584e-05
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.2895530 0.6509662
## sample estimates:
##       cor 
## 0.4911129
```

```
## [1] 0.4780852
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  free_meta_data$D1 and free_meta_data$HNA_percent
## t = 4.4886, df = 68, p-value = 2.843e-05
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.2738779 0.6410457
## sample estimates:
##       cor 
## 0.4780852
```

```
## [1] 0.7018938
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  free_meta_data$D0 and free_meta_data$HNA_percent
## t = 8.126, df = 68, p-value = 1.294e-11
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.5591364 0.8042289
## sample estimates:
##       cor 
## 0.7018938
```

```
## [1] 0.7084204
```

```
## 
## 	Pearson's product-moment correlation
## 
## data:  free_meta_data$D0_chao and free_meta_data$HNA_percent
## t = 8.2769, df = 68, p-value = 6.883e-12
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.5679950 0.8087668
## sample estimates:
##       cor 
## 0.7084204
```



# Is there a relationship between diversity and *fractionated* bacterial production?

```
## 
## Call:
## lm(formula = frac_bacprod ~ D2, data = wholeparticle_2015_df)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -7.2971 -2.1906 -0.2354  1.3639  7.6078 
## 
## Coefficients:
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  0.02044    2.33018   0.009 0.993173    
## D2           0.26298    0.05084   5.173 0.000417 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 4.517 on 10 degrees of freedom
## Multiple R-squared:  0.7279,	Adjusted R-squared:  0.7007 
## F-statistic: 26.76 on 1 and 10 DF,  p-value: 0.0004174
```

```
## Error in eval(expr, envir, enclos): could not find function "rlm"
```

```
## Error in summary(rlm_part_otu_D2_stats): object 'rlm_part_otu_D2_stats' not found
```

```
## Error in car::Anova(rlm_part_otu_D2_stats): object 'rlm_part_otu_D2_stats' not found
```

<img src="Figures/cached/fractional_production_vs_diversity-1.png" style="display: block; margin: auto;" /><img src="Figures/cached/fractional_production_vs_diversity-2.png" style="display: block; margin: auto;" />



# Phylum production analysis




```
## Error in ggplot(alpha_comb_final, aes(x = Estimate, y = frac_bacprod)): object 'alpha_comb_final' not found
```





```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'alpha_comb_final' not found
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'actinos' not found
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'actinos' not found
```

```
## Error in ggplot(actinos, aes(x = Estimate, y = frac_bacprod)): object 'actinos' not found
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'alpha_comb_final' not found
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'bacteroidetes' not found
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'bacteroidetes' not found
```

```
## Error in ggplot(bacteroidetes, aes(x = Estimate, y = frac_bacprod)): object 'bacteroidetes' not found
```




# How do the samples relate to each other? 
## PCOA
<img src="Figures/cached/PCOA-1.png" style="display: block; margin: auto;" />



#  Stacked Bar plots for taxonomic structure
## OTU  Phylum Stacked Bar 





# Environmental Conditions 







# Calculate the Absolute abundances!



### HNA VS LNA 









