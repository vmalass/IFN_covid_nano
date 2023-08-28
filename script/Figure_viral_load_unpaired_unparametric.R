Test stat non parametrique et non apparié charge virale

library(tidyverse)
library(ggsignif)
library(openxlsx)
library(reshape2)

#-- Importer les données
rm(list = ls())
ifn <- read.xlsx( "/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/Dosage_IFN.xlsx",4)

viral <- ifn[,1:5]
viral <- na.omit(viral)
count(viral, time, Reponse)

#-- Toutes cytokines confondues, tester si un groupe a été significativement plus prélevé que le 2eme groupe
#à certains instants.
chisq.test(table(viral$Reponse , viral$time))
##### Résultat du test négatif. 

#-- Tous groupes confondus, tester si une cytokine a été significativement plus mesurée que les autres
#à certains instants.
chisq.test(table(viral$charge_virale , viral$time))
#####  Résultat du test négatif. 


with(viral , wilcox.test(charge_virale ~ Reponse))

ggplot(viral , aes(x =  Reponse , y = charge_virale))+
  geom_boxplot( outlier.shape = NA)+ #  outlier.shape = NA
  geom_jitter()  +# aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for viral load", 
       y =  "viral load (log10copies /1 000 000 cellules)") +
  theme_classic() 

  