Test stat frequence cellules TEMRA CD57 + et - à V1 et V4








# import datat------------------------------------------------------------------
rm(list = ls())
data <- read.xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/cyto_V1_V4_S1_TEMRA_Agspe.xlsx", sheet = 1)
# TEMRA_pos <- data[,c("Groupe", "CD8_TEMRA_CD57+", "Time")]
# TEMRA_neg <- data[,c("Groupe", "CD8_TEMRA_CD57-", "Time")]
TEMRA <- data[,c("Groupe", "CD8_TEMRA_CD57+", "CD8_TEMRA_CD57-", "Time")]

# TEMRA_pos <- na.omit(TEMRA_pos)
# TEMRA_neg <- na.omit(TEMRA_neg)
TEMRA <- na.omit(TEMRA)

dplyr::count(TEMRA, Time, Groupe)

# CD8 TEMRA CD57+ --------------------------------------------------------------

ggplot(TEMRA , aes(x = Time, 
                   y = `CD8_TEMRA_CD57+`, 
                   color =  Groupe)) +
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
  theme_classic() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(y_position = c(20, 10), 
              xmin = c(0.8, 1.8), 
              xmax = c(1.2, 2.2),
              annotation = c("0.11", "0.07")) + 
  labs(title = "Quantity of CD8 TEMRA CD57+ all patient and time points",
       subtitle = "VT1 (13 R / 5 RP) VT4 (10 R / 5 RP)",
       y =  "Frequence of CD8 TEMRA CD57+",
       x = "")

#-- Tout les times confondues, tester si un groupe a été significativement plus prélevé que le 2eme groupe
#à certains instants.
chisq.test(table(TEMRA$Groupe , TEMRA$Time))
##### Résultat du test négatif. 

#-- Test statistique, non apparié et non paramétrique (transfrmation log non utile dans ce cas).
# Est-ce que la frequence de TEMRA est différente entre les 2 groupes, R et RP ?
with(TEMRA , wilcox.test(`CD8_TEMRA_CD57+` ~ Groupe))
#####  Résultat du test positif.

dplyr::count(TEMRA, Groupe)

ggplot(TEMRA , aes(x =  Groupe , y = `CD8_TEMRA_CD57+`, color = Groupe))+
  geom_boxplot( outlier.shape = NA) + #  outlier.shape = NA
  # geom_point() +
  geom_jitter() + # aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for CD8_TEMRA_CD57+ frequency all time points and patients", 
       subtitle = "R 23 / RP 10",
       y = "Frequence of CD8_TEMRA_CD57+",
       x = "") +
  theme_classic() 

T_F <- TEMRA$Time %in% "VT4"
VT1 <- TEMRA[T_F == F,]
VT4 <- TEMRA[T_F == T,]

with(VT1 , wilcox.test(`CD8_TEMRA_CD57+` ~ Groupe))
with(VT4 , wilcox.test(`CD8_TEMRA_CD57+` ~ Groupe))


# CD8 TEMRA CD57- --------------------------------------------------------------
ggplot(TEMRA , aes(x = Time, 
                   y = `CD8_TEMRA_CD57-`, 
                   fill =  Groupe)) +
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  # geom_jitter() +
  theme_classic() +
  labs(title = "Quantity of CD8 TEMRA CD57- all patient and time points",
       subtitle = "VT4 (13 R / 5 RP) S1 (10 R / 5 RP)",
       y =  "Frequence of CD8 TEMRA CD57-",
       x = "")

#-- Test statistique, non apparié et non paramétrique (transfrmation log non utile dans ce cas).
# Est-ce que la frequence de TEMRA est différente entre les 2 groupes, R et RP ?
with(TEMRA , wilcox.test(`CD8_TEMRA_CD57-` ~ Groupe))
#####  Résultat du test positif.

dplyr::count(TEMRA, Groupe)

ggplot(TEMRA , aes(x =  Groupe , y = `CD8_TEMRA_CD57-`))+
  geom_boxplot( outlier.shape = NA) + #  outlier.shape = NA
  # geom_point() +
  geom_jitter() + # aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for CD8_TEMRA_CD57- frequency all time points and patients", 
       subtitle = "R 23 / RP 10",
       y = "Frequence of CD8_TEMRA_CD57-") +
  theme_classic() 

with(VT1 , wilcox.test(`CD8_TEMRA_CD57-` ~ Groupe, paired = F))
with(VT4 , wilcox.test(`CD8_TEMRA_CD57-` ~ Groupe, paired = F))

