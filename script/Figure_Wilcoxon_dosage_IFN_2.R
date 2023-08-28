Test stat dosage IFN unpaired et unparametric


library(tidyverse)
library(ggsignif)
library(openxlsx)
library(reshape2)

#-- Importer les données
rm(list = ls())
ifn <- read.xlsx( "/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/Dosage_IFN.xlsx",3)

#-- Sélectionner les colonnes pertientes dans la table.
ifn <- ifn[ , c('IFN_l1_IL_29' , 'IFN_gamma' , 'IFN_beta' ,'Reponse' , 'temps')]
ifn <- dplyr::rename(ifn, "IFN_lambda" = "IFN_l1_IL_29")

#-- Empiler la table et supprimer
melted <- reshape2::melt(ifn, id.vars = c('Reponse','temps' ))
colnames(melted) <- c("groupe" ,"temps", "cytokine","concentration")

#-- Toutes cytokines confondues, tester si un groupe a été significativement plus prélevé que le 2eme groupe
#à certains instants.
chisq.test(table(melted$groupe , melted$temps))
##### Résultat du test négatif. 

#-- Tous groupes confondus, tester si une cytokine a été significativement plus mesurée que les autres
#à certains instants.
chisq.test(table(melted$cytokine , melted$temps))
#####  Résultat du test négatif. 

#-- Test statistique, non apparié et non paramétrique (transfrmation log non utile dans ce cas).
# Est-ce que la concentration cytokinique (IFN_l1_IL_29, IFN_gamma & IFN_beta confondues) est 
# différente entre les 2 groupes, R et RP ?
with(melted , wilcox.test(concentration ~ groupe))
#####  Résultat du test positif.

dplyr::count(melted, cytokine, groupe)

ggplot(melted , aes(x = factor(cytokine, levels = c("IFN_beta", "IFN_gamma", "IFN_lambda")) , y = log(concentration),  color = groupe))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
  theme_classic() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(y_position = c(5.5, 8.2, 5.5), 
              xmin = c(0.8, 1.8, 2.8), 
              xmax = c(1.2, 2.2, 3.2),
              annotation = c("0.21", "**", "0.16")) + 
  labs(title = "Wilcoxon test for IFN quantity",
       subtitle = "23 R / 12 RP",
       y = "IFN in log (pg/mL)",
       x = "") + 
  theme(legend.position = c(0.9,0.9))



dplyr::count(melted, groupe)

ggplot(melted , aes(x =  groupe , y = log(concentration),  color = groupe))+
  geom_boxplot( outlier.shape = NA) + #  outlier.shape = NA
  # geom_point() +
  geom_jitter() + # aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for IFN quantity", 
       subtitle = "69 R / 36 RP",
       y = "IFN (pg/mL)",
       x = "") +
  theme_classic() + 
  theme(legend.position = c(0.9,0.99))




# Table par cytokine------------------------------------------------------------
T_F <- melted$cytokine %in% "IFN_lambda"
IFN_lambda <- melted[T_F == T,]

T_F <- melted$cytokine %in% "IFN_gamma"
IFN_gamma <- melted[T_F == T,]

T_F <- melted$cytokine %in% "IFN_beta"
IFN_beta <- melted[T_F == T,]


with(IFN_lambda , wilcox.test(concentration ~ groupe))
with(IFN_gamma , wilcox.test(concentration ~ groupe))
with(IFN_beta , wilcox.test(concentration ~ groupe))


ggplot(IFN_lambda , aes(x = groupe , y = log(concentration)))+
  geom_boxplot()+
  geom_jitter() +
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  theme_classic() +
  labs(title = "Wilcoxon test for IFN l1 IL29 quantity", 
       y = "IFN l1 IL29 (pg/mL)")


ggplot(IFN_gamma , aes(x = groupe , y = log(concentration)))+
  geom_boxplot()+
  geom_jitter() +
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  theme_classic() +
  labs(title = "Wilcoxon test for IFN gamma quantity", 
       y = "IFN gamma (pg/mL)",
       x = "")


ggplot(IFN_beta , aes(x = groupe , y = log(concentration)))+
  geom_boxplot()+
  geom_jitter() +
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  theme_classic() +
  labs(title = "Wilcoxon test for IFN beta quantity", 
       y = "IFN beta (pg/mL)")
