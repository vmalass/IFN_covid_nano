Test stat frequence TEMRA chez tout les patients à tout les times points







# import datat------------------------------------------------------------------
rm(list = ls())
data <- read.xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/cyto_V1_V4_S1_TEMRA_Agspe.xlsx", sheet = 1)
TEMRA <- data[,c("Groupe", "CD8_TEMRA_totaux", "Time")]

# Kinetics TEMRA----------------------------------------------------------------
TEMRA <- na.omit(TEMRA)
dplyr::count(TEMRA, Time, Groupe)

ggplot(TEMRA , aes(x = factor(Time, levels = c("VT1", "VT2", "VT4", "S1")), 
                          y = CD8_TEMRA, 
                          fill =  Groupe)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  # geom_jitter() +
  theme_classic() +
  labs(title = "Quantity of TEMRA all patients and times points",
       subtitle = "V1 (13 R / 5 RP) V2 (10 R / 3 RP) V4 (10 R / 5 RP) S1 (11 R / 5 RP)",
       y =  "Frequence of TEMRA",
       x = "")

#-- Tout les times confondues, tester si un groupe a été significativement plus prélevé que le 2eme groupe
#à certains instants.
chisq.test(table(TEMRA$Groupe , TEMRA$Time))
##### Résultat du test négatif. 

#-- Test statistique, non apparié et non paramétrique (transfrmation log non utile dans ce cas).
# Est-ce que la frequence de TEMRA est différente entre les 2 groupes, R et RP ?
with(TEMRA , wilcox.test(CD8_TEMRA ~ Groupe))
#####  Résultat du test positif.

dplyr::count(TEMRA, Groupe)

ggplot(TEMRA , aes(x =  Groupe , y = CD8_TEMRA))+
  geom_boxplot( outlier.shape = NA) + #  outlier.shape = NA
  # geom_point() +
  geom_jitter() + # aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for TEMRA frequency all time points and patients", 
       subtitle = "R 44 / RP 18",
       y = "Frequence of TEMRA") +
  theme_classic() 

# Wilcoxon test of eatch time points--------------------------------------------
T_F <- TEMRA$Time %in% "VT1"
VT1 <- TEMRA[T_F == T,]

T_F <- TEMRA$Time %in% "VT2"
VT2 <- TEMRA[T_F == T,]

T_F <- TEMRA$Time %in% "VT4"
VT4 <- TEMRA[T_F == T,]

T_F <- TEMRA$Time %in% "S1"
S1 <- TEMRA[T_F == T,]

wilcox_test(VT1,formula(paste0( "CD8_TEMRA", "~ Groupe"))) %>% 
  add_significance()
wilcox_test(VT2,formula(paste0( "CD8_TEMRA", "~ Groupe"))) %>% 
  add_significance()
wilcox_test(VT4,formula(paste0( "CD8_TEMRA", "~ Groupe"))) %>% 
  add_significance()
wilcox_test(S1,formula(paste0( "CD8_TEMRA", "~ Groupe"))) %>% 
  add_significance()

