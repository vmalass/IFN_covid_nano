Test stat frequence TEMRA chez  les patients avec un suivi longitudinal à VT1-2-4 et S1






# import datat------------------------------------------------------------------
rm(list = ls())
data <- read.xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/cyto_V1_V4_S1_TEMRA_Agspe.xlsx", sheet = 2)
TEMRA <- data[,c("Groupe", "CD8_TEMRA", "Time")]

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
  labs(title = "Quantity of TEMRA patient with full longitudinal follow-up (V1-2-4 S1)",
       subtitle = "R = 7 / RP = 3",
       y =  "Frequence of TEMRA",
       x = "")

#-- Tout les times confondues, tester si un groupe a été significativement plus prélevé que le 2eme groupe
#à certains instants.
chisq.test(table(TEMRA$Groupe , TEMRA$Time))
##### Résultat du test négatif. 

#-- Test statistique, non apparié et non paramétrique (transfrmation log non utile dans ce cas).
# Est-ce que la frequence de TEMRA est différente entre les 2 groupes, R et RP ?
with(TEMRA , t.test(CD8_TEMRA ~ Groupe))
#####  Résultat du test positif.

dplyr::count(TEMRA, Groupe)

ggplot(TEMRA , aes(x =  Groupe , y = CD8_TEMRA))+
  geom_boxplot( outlier.shape = NA) + #  outlier.shape = NA
  # geom_point() +
  geom_jitter() + # aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for TEMRA frequency all time points and patients", 
       subtitle = "R 28 / RP 12",
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

t_test(VT1,formula(paste0( "CD8_TEMRA", "~ Groupe"))) %>% 
  add_significance()
t_test(VT2,formula(paste0( "CD8_TEMRA", "~ Groupe"))) %>% 
  add_significance()
t_test(VT4,formula(paste0( "CD8_TEMRA", "~ Groupe"))) %>% 
  add_significance()
t_test(S1,formula(paste0( "CD8_TEMRA", "~ Groupe"))) %>% 
  add_significance()

