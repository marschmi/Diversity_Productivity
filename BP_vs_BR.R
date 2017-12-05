# Reinthaler et al., 2005
# Relationship between Bacterioplankton Richness, Respiration, and Production in the Southern North Sea

bac_prod <- c(1.75, 0.53, 0.44, 0.22, 0.13, 0.03, 0.40)
bac_resp <- c(2.4, 1.6, 1.7, 1.6, 0.7, 0.8, 1.4)

df <- data.frame(bac_prod, bac_resp)

## 1. Run the linear model
lm_BR_BP <- lm(bac_resp ~ bac_prod, data = df)

## 2. Extract the R2 and p-value from the linear model: 
lm_BR_BP_lab <- paste("atop(R^2 ==", round(summary(lm_BR_BP)$adj.r.squared, digits = 2), ",",
                                   "p ==", round(unname(summary(lm_BR_BP)$coefficients[,4][2]), digits = 3), ")")

ggplot(df, aes(x = bac_prod, bac_resp)) + 
  geom_point() + theme_bw() +
  ggtitle("Reinthaler et al., 2005, AEM") + 
  labs(x = "Bacterial Production \n (mmol C/m3/day)", y = "Bacterial Respiration \n(mmol C/m3/day)") + 
  geom_smooth(method = "lm") +
  annotate("text", x=1.25, y=1, label=lm_BR_BP_lab, parse = TRUE, color = "#424645", size = 4) +
  theme(legend.title = element_blank(), legend.position ="bottom", 
        legend.text = element_text(size = 14))

ggsave("BP_vs_BR.jpeg", width = 3.5, height =3)