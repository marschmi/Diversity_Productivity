# Community wide
percent_dormant_taxa <- richness %>%
  dplyr::filter(fraction == "Particle") %>%
  mutate(percent_dormant_taxa = 295/mean)

# linear model 
summary(lm(frac_bacprod ~ percent_dormant_taxa, data = percent_dormant_taxa))

# Plot it
ggplot(percent_dormant_taxa, aes(x = percent_dormant_taxa, y = frac_bacprod, fill = fraction)) + 
  geom_smooth(data=filter(percent_dormant_taxa, fraction == "Particle"), method='lm', color = "#FF6600", fill = "#FF6600") + 
  geom_point(aes(shape = season), size = 3.5) + xlab("295/Richness") + ylab("Community-Wide Production") + 
  scale_fill_manual(values = fraction_colors) +
  scale_shape_manual(values = season_shapes) 



# linear model 
summary(lm(fracprod_per_cell_noinf ~ percent_dormant_taxa, data = percent_dormant_taxa))

# Per capita ggplot(percent_dormant_taxa, aes(x = percent_dormant_taxa, y = frac_bacprod, fill = fraction)) + 
ggplot(percent_dormant_taxa, aes(x = percent_dormant_taxa, y = fracprod_per_cell_noinf, fill = fraction)) + 
  geom_smooth(data=filter(percent_dormant_taxa, fraction == "Particle"), method='lm', color = "#FF6600", fill = "#FF6600") + 
  geom_point(aes(shape = season), size = 3.5) + xlab("295/Richness") + ylab("Per-Capita Production") + 
  scale_fill_manual(values = fraction_colors) +
  scale_shape_manual(values = season_shapes) 

