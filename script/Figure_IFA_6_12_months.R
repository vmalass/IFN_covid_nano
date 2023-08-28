Etude du titre en IFA 6 et 12 mois

library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(rstatix)
library(cowplot)

# import datat------------------------------------------------------------------
rm(list = ls())
data <- read.xlsx("data/titre_virale_IFI27_Ac.xlsx", sheet = 8)
data$numero_patient<- as.factor(data$numero_patient)

T_F <- data$Groupe %in% c("R", "RP")
data_2 <- data[T_F == T,]

## IFA.RBD----------------------------------------------------------------------

ggplot(data = data_2, 
       aes(x = time_point, 
           y = IFA.RBD, 
           color = Groupe))+ 
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  labs(title = "Wilcoxon test for title IFA RBD (spike)", 
       x = "Time point", 
       y = "quantity of IFA RBD")+
  theme_classic() + 
  theme(legend.position = c(0.95,0.8))

T_F <- data_2$time_point %in% "S1"
data_VT1 <- data_2[T_F==T,]

T_F <- data_2$time_point %in% "S2"
data_VT2 <- data_2[T_F==T,]


data_VT1<- data_VT1%>% 
  reorder_levels(Groupe, order = c( "R", "RP"))

# Test stat S1
res.wilcox <- data_VT1 %>%
  wilcox_test(IFA.RBD ~ Groupe) %>%
  add_significance()
res.wilcox
res.wilcox <- data_VT1 %>%
  wilcox_test(IFA.MBGP ~ Groupe) %>%
  add_significance()
res.wilcox
res.wilcox <- data_VT1 %>%
  wilcox_test(IFA.NC ~ Groupe) %>%
  add_significance()
res.wilcox

data_VT2<- data_VT2%>% 
  reorder_levels(Groupe, order = c( "R", "RP"))

# Test stat S2
res.wilcox <- data_VT2 %>%
  wilcox_test(IFA.RBD ~ Groupe) %>%
  add_significance()
res.wilcox
res.wilcox <- data_VT2 %>%
  wilcox_test(IFA.MBGP ~ Groupe) %>%
  add_significance()
res.wilcox
res.wilcox <- data_VT2 %>%
  wilcox_test(IFA.NC ~ Groupe) %>%
  add_significance()
res.wilcox

# Test stat S1 with S2
res.wilcox <- data_2 %>%
  wilcox_test(IFA.RBD ~ Groupe) %>%
  add_significance()
res.wilcox
res.wilcox <- data_2 %>%
  wilcox_test(IFA.MBGP ~ Groupe) %>%
  add_significance()
res.wilcox
res.wilcox <- data_2 %>%
  wilcox_test(IFA.NC ~ Groupe) %>%
  add_significance()
res.wilcox

## IFA.MBGP---------------------------------------------------------------------

ggplot(data = data_2, 
       aes(x = time_point, 
           y = IFA.MBGP, 
           color = Groupe))+ 
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  labs(title = "Wilcoxon test for title IFA MBGP (membrane)", 
       x = "Time point", 
       y = "quantity of IFA MBGP")+
  theme_classic() + 
  theme(legend.position = c(0.95,0.8))

# Test stat
res.wilcox <- data_VT1 %>%
  wilcox_test(IFA.MBGP ~ Groupe) %>%
  add_significance()
res.wilcox

# Test stat
res.wilcox <- data_VT2 %>%
  wilcox_test(IFA.MBGP ~ Groupe) %>%
  add_significance()
res.wilcox

# Test stat
res.wilcox <- data_2 %>%
  wilcox_test(IFA.MBGP ~ Groupe) %>%
  add_significance()
res.wilcox
## IFA.NC-----------------------------------------------------------------------

ggplot(data = data_2, 
       aes(x = time_point, 
           y = IFA.NC, 
           color = Groupe))+ 
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  labs(title = "Wilcoxon test for title IFA NC (nucleocapsid)", 
       x = "Time point", 
       y = "quantity of IFA NC")+
  theme_classic() + 
  theme(legend.position = c(0.95,0.8))

# Test stat
res.wilcox <- data_VT1 %>%
  wilcox_test(IFA.NC ~ Groupe) %>%
  add_significance()
res.wilcox

# Test stat
res.wilcox <- data_VT2 %>%
  wilcox_test(IFA.NC ~ Groupe) %>%
  add_significance()
res.wilcox

# Test stat
res.wilcox <- data_2 %>%
  wilcox_test(IFA.NC ~ Groupe) %>%
  add_significance()
res.wilcox

# IFA 6 et 12 mois ensemble-----------------------------------------------------
# test unpaired and unparametric------------------------------------------------

IFA <- data_2[,c("Groupe", "time_point", "IFA.RBD", "IFA.MBGP", "IFA.NC")]

melted <- reshape2::melt(IFA, id.vars = c('Groupe','time_point' ))
colnames(melted) <- c("Groupe" ,"time_point", "IFA","concentration")

dplyr::count(melted, IFA, Groupe)

ggplot(melted , aes(x = factor(IFA, levels = c("IFA.RBD", "IFA.MBGP", "IFA.NC")), 
                    y = log(concentration), 
                    color =  Groupe)) +
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
  theme_classic() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(y_position = c(2.5, 2.5, 2.5), 
              xmin = c(0.8, 1.8, 2.8), 
              xmax = c(1.2, 2.2, 3.2),
              annotation = c("NS", "NS", "NS")) + 
  theme_classic() +
  labs(title = "All time points (6, 12 months) and IFA (RBD, MBGP, NC)", 
       subtitle = "33 R / 14 RP",
       y =  "quantity of IFA",
       x = "")

## Tous les IFA 6 and 12 months-------------------------------------------------

#-- Toutes cytokines confondues, tester si un groupe a été significativement plus prélevé que le 2eme groupe
#à certains instants.
chisq.test(table(melted$Groupe , melted$time_point))
#Résultat du test négatif. 

#-- Tous groupes confondus, tester si une IFA a été significativement plus mesurée que les autres
#à certains instants.
chisq.test(table(melted$IFA , melted$time_point))
#Résultat du test négatif. 

with(melted , wilcox.test(concentration ~ Groupe))

dplyr::count(melted, Groupe)

ggplot(melted , aes(x =  Groupe , y = concentration, color = Groupe))+
  geom_boxplot( outlier.shape = NA)+ #  outlier.shape = NA
  geom_jitter()  +# aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for IFA (RBD, MBGP, NC)", 
       subtitle = "99 R / 42 RP",
       y =  "quantity of IFA",
       x = "") +
  theme_classic() 

## IFA MGBP & NC 6 and 12 months------------------------------------------------

T_F <- melted$IFA %in% "IFA.RBD"
melted_2 <- melted[T_F == F,]

#-- Toutes cytokines confondues, tester si un groupe a été significativement plus prélevé que le 2eme groupe
#à certains instants.
chisq.test(table(melted_2$Groupe , melted_2$time_point))
#Résultat du test négatif. 

#-- Tous groupes confondus, tester si une IFA a été significativement plus mesurée que les autres
#à certains instants.
chisq.test(table(melted_2$IFA , melted_2$time_point))
#Résultat du test négatif. 

with(melted_2 , wilcox.test(concentration ~ Groupe))

dplyr::count(melted_2, Groupe)

ggplot(melted_2 , aes(x =  Groupe , y = concentration, color = Groupe))+
  geom_boxplot( outlier.shape = NA)+ #  outlier.shape = NA
  geom_jitter()  +# aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for IFA (MBGP, NC)", 
       subtitle = "66 R / 28 RP",
       y =  "quantity of IFA",
       x = "") +
  theme_classic() 


# IFA 6 mois--------------------------------------------------------------------
T_F <- melted$IFA %in% "IFA.RBD"
melted_RBD <- melted[T_F==F,]
T_F <- melted_RBD$time_point %in% "S1"
melted_S1 <- melted_RBD[T_F==T,]
melted_S2 <- melted_RBD[T_F==F,]

dplyr::count(melted_S1, IFA, time_point, Groupe)

ggplot(melted_S1 , aes(x = factor(IFA, levels = c("IFA.MBGP", "IFA.NC")), 
                    y = log(concentration), 
                    color =  Groupe)) +
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
  theme_classic() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(y_position = c(2.5, 2.5), 
              xmin = c(0.8, 1.8), 
              xmax = c(1.2, 2.2),
              annotation = c("0.08", "0.16")) + 
  theme_classic() +
  labs(title = "IFA (MBGP, NC) at 6 months", 
       subtitle = "19 R and 7 RP",
       y =  "quantity of IFA",
       x = "")


wilcox_test(melted_S1, IFA.RBD ~ Groupe)
#-- Tous groupes confondus, tester si une IFA a été significativement plus mesurée que les autres
#à certains instants.
chisq.test(table(melted_S1$IFA , melted_S1$time_point))
#Résultat du test négatif. 

with(melted_S1 , wilcox.test(concentration ~ Groupe))

ggplot(melted_S1 , aes(x =  Groupe , y = concentration))+
  geom_boxplot( outlier.shape = NA)+ #  outlier.shape = NA
  geom_jitter()  +# aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for IFA (RBD, MBGP, NC) at 6 months", 
       y =  "quantity of IFA",
       x = "") +
  theme_classic() 

## Avec IFA MBGP & NC-----------------------------------------------------------
T_F <- melted_S1$IFA %in% "IFA.RBD"
melted_S1 <- melted_S1[T_F == F,]

#-- Tous groupes confondus, tester si une IFA a été significativement plus mesurée que les autres
#à certains instants.
chisq.test(table(melted_S1$IFA , melted_S1$time_point))
#Résultat du test positif??????? 
dplyr::count(melted_S1, IFA, time_point, Groupe)

with(melted_S1 , wilcox.test(concentration ~ Groupe))

ggplot(melted_S1 , aes(x =  Groupe , y = concentration, color = Groupe))+
  geom_boxplot( outlier.shape = NA)+ #  outlier.shape = NA
  geom_jitter()  +# aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  # geom_boxplot() +
  # geom_point(position = position_dodge(width=0.75)) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for IFA (MBGP, NC) at 6 months", 
       subtitle = "19 R and 7 RP",
       y =  "quantity of IFA",
       x = "") +
  theme_classic() + 
  theme(legend.position = c(0.95,0.8))



# IFA 12 mois-------------------------------------------------------------------

ggplot(melted_S2 , aes(x = factor(IFA, levels = c("IFA.RBD", "IFA.MBGP", "IFA.NC")), 
                       y = log(concentration), 
                       fill =  Groupe)) +
  geom_boxplot()+
  geom_jitter() +
  theme_classic() +
  labs(title = "IFA (RBD, MBGP, NC) at 12 months", 
       y =  "quantity of IFA",
       x = "")

#-- Tous groupes confondus, tester si une IFA a été significativement plus mesurée que les autres
#à certains instants.
chisq.test(table(melted_S2$IFA , melted_S2$time_point))
#Résultat du test négatif. 

with(melted_S2 , wilcox.test(concentration ~ Groupe))

ggplot(melted_S2 , aes(x =  Groupe , y = concentration))+
  geom_boxplot( outlier.shape = NA)+ #  outlier.shape = NA
  geom_jitter()  +# aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for IFA (RBD, MBGP, NC) at 12 months", 
       y =  "quantity of IFA",
       x = "") +
  theme_classic() 

## Avec IFA MBGP & NC-----------------------------------------------------------
T_F <- melted_S2$IFA %in% "IFA.RBD"
melted_S2 <- melted_S2[T_F == F,]

#-- Tous groupes confondus, tester si une IFA a été significativement plus mesurée que les autres
#à certains instants.
chisq.test(table(melted_S2$IFA , melted_S2$time_point))
#Résultat du test positif??????? 
dplyr::count(melted_S2, IFA, time_point)

with(melted_S2 , wilcox.test(concentration ~ Groupe))

ggplot(melted_S2 , aes(x =  Groupe , y = concentration))+
  geom_boxplot( outlier.shape = NA)+ #  outlier.shape = NA
  geom_jitter()  +# aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for IFA (MBGP, NC) at 6 months", 
       y =  "quantity of IFA",
       x = "") +
  theme_classic() 

