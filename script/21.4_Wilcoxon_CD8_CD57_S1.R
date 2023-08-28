Test stat frequence cellules CD8 CD57 + et - à 6 mois









# import datat------------------------------------------------------------------
rm(list = ls())
data <- read.xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/cyto_V1_V4_S1_TEMRA_Agspe.xlsx", sheet = 7)
# TEMRA_pos <- data[,c("Groupe", "CD8_TEMRA_CD57+", "Time")]
# TEMRA_neg <- data[,c("Groupe", "CD8_TEMRA_CD57-", "Time")]
TEMRA <- data[,c("Groupe", "CD8_CD57+", "CD8_CD57-", "Time")]

# TEMRA_pos <- na.omit(TEMRA_pos)
# TEMRA_neg <- na.omit(TEMRA_neg)
TEMRA <- na.omit(TEMRA)

dplyr::count(TEMRA, Time, Groupe)

# CD8 CD57+ --------------------------------------------------------------------

with(TEMRA , wilcox.test(`CD8_CD57+` ~ Groupe))

ggplot(TEMRA , aes(x = Time, 
                   y = `CD8_CD57+`, 
                   fill =  Groupe)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  # geom_jitter() +
  theme_classic() +
  labs(title = "Quantity of CD8 CD57+ all patient at 6 months",
       subtitle = "R 11 / RP 5",
       y =  "Frequence of CD8 CD57+",
       x = "")

# CD8 CD57- --------------------------------------------------------------------

with(TEMRA , wilcox.test(`CD8_CD57-` ~ Groupe))

ggplot(TEMRA , aes(x = Time, 
                   y = `CD8_CD57-`, 
                   fill =  Groupe)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  # geom_jitter() +
  theme_classic() +
  labs(title = "Quantity of CD8 CD57- all patient at 6 months",
       subtitle = "R 11 / RP 5",
       y =  "Frequence of CD8 CD57+",
       x = "")
