Test stat frequence cellules antigène spé chez tout les patients à tout les times points (V4 et 6 mois)







# import datat------------------------------------------------------------------
rm(list = ls())
data <- read.xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/cyto_V1_V4_S1_TEMRA_Agspe.xlsx", sheet = 1)
Ag_spe <- data[,c("Groupe", "Ag_spe", "Time")]

Ag_spe <- na.omit(Ag_spe)

# Kinetics Ag spe----------------------------------------------------------------
dplyr::count(Ag_spe, Time, Groupe)

ggplot(Ag_spe , aes(x = factor(Time, levels = c("VT4", "S1")), 
                   y = Ag_spe, 
                   fill =  Groupe)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  # geom_jitter() +
  theme_classic() +
  labs(title = "Quantity of Ag_spe all patient and time points",
       subtitle = "VT4 (13 R / 5 RP) S1 (11 R / 5 RP)",
       y =  "Frequence of Ag_spe",
       x = "")

#-- Tout les times confondues, tester si un groupe a été significativement plus prélevé que le 2eme groupe
#à certains instants.
chisq.test(table(Ag_spe$Groupe , Ag_spe$Time))
##### Résultat du test négatif. 

#-- Test statistique, non apparié et non paramétrique (transfrmation log non utile dans ce cas).
# Est-ce que la frequence de Ag_spe est différente entre les 2 groupes, R et RP ?
with(Ag_spe , wilcox.test(Ag_spe ~ Groupe))
#####  Résultat du test positif.

dplyr::count(Ag_spe, Groupe)

ggplot(Ag_spe , aes(x =  Groupe , y = Ag_spe))+
  geom_boxplot( outlier.shape = NA) + #  outlier.shape = NA
  # geom_point() +
  geom_jitter() + # aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for Ag_spe frequency all time points and patients", 
       subtitle = "R 24 / RP 10",
       y = "Frequence of Ag_spe") +
  theme_classic() 

# Wilcoxon test of eatch time points--------------------------------------------
T_F <- Ag_spe$Time %in% "VT4"
VT4 <- Ag_spe[T_F == T,]

T_F <- Ag_spe$Time %in% "S1"
S1 <- Ag_spe[T_F == T,]

wilcox_test(VT4,formula(paste0( "Ag_spe", "~ Groupe"))) %>% 
  add_significance()
wilcox_test(S1,formula(paste0( "Ag_spe", "~ Groupe"))) %>% 
  add_significance()
