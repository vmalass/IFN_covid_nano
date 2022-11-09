library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(rstatix)

# import datat----
rm(list = ls())
data <- read.xlsx("data/titre_virale_IFI27_Ac.xlsx", sheet = 4)
data_months <- read.xlsx("data/titre_virale_IFI27_Ac.xlsx", sheet = 5)
data$numero_patient<- as.factor(data$numero_patient)

# creation des sous matrice pour les groupes-----
new <- data$Groupe %in% "NR"
groupe_NR <- data[new==T,]
new <- data$Groupe %in% "R"
groupe_R <- data[new==T,]
new <- data$Groupe %in% "RP"
groupe_RP <- data[new==T,]

new <- data_months$months %in% "S1"
data_S1 <- data_months[new==T,]
new <- data_months$months %in% "S2"
data_S2 <- data_months[new==T,]

numgp_NR <- groupe_NR$numero_patient
numgp_R <- groupe_R$numero_patient
numgp_RP <- groupe_RP$numero_patient

# visualisation----
ggplot(data = groupe_NR, 
       aes_string(x="day", 
                  y= "Vidas", 
                  color = "numero_patient ", 
                  group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("Titre en anti-corps du groupe NR")+
  scale_x_continuous(name = "jours post infection",
                     limits = c(0, 75))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15))   

ggplot(data = groupe_R, 
       aes_string(x="day", 
                  y= "Vidas", 
                  color = "numero_patient ", 
                  group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("Titre en anti-corps du groupe R")+
  scale_x_continuous(name = "jours post infection",
                     limits = c(0, 75))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15))   

ggplot(data = groupe_RP, 
       aes_string(x="day", 
                  y= "Vidas", 
                  color = "numero_patient ", 
                  group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("Titre en anti-corps du groupe RP")+
  scale_x_continuous(name = "jours post infection",
                     limits = c(0, 75))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15))   


# Bornage 0 à 20 jours-----
ggplot(data = groupe_NR, 
       aes_string(x="day", 
                  y= "Vidas", 
                  color = "numero_patient ", 
                  group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("Titre en anti-corps du groupe NR")+
  scale_x_continuous(name = "jorus post infection",
                     limits = c(0, 20))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15))   

ggplot(data = groupe_R, 
       aes_string(x="day", 
                  y= "Vidas", 
                  color = "numero_patient ", 
                  group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("Titre en anti-corps du groupe R")+
  scale_x_continuous(name = "jorus post infection",
                     limits = c(0, 20))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15))   

ggplot(data = groupe_RP, 
       aes_string(x="day", 
                  y= "Vidas", 
                  color = "numero_patient ", 
                  group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("Titre en anti-corps du groupe RP")+
  scale_x_continuous(name = "jorus post infection",
                     limits = c(0, 20))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15))   

# Bornage pour nouveau time point ----
data <- arrange(data, day)
data_VT1 <- data[1:34,]
data_VT2 <- data[35:76,]
data_VT3 <- data[77:119,]
data_VT4 <- data[120:157,]
data_VT5 <- data[158:195,]
data_VT6 <- data[196:232,]
data_VT7 <- data[233:266,]
data_VT8 <- data[267:279,]

## VT1----
datat_VT1_NR<- data_VT1$Groupe %in% "NR"
datat_VT1_NR<- data_VT1[datat_VT1_NR==T,]
datat_VT1_R<- data_VT1$Groupe %in% "R"
datat_VT1_R<- data_VT1[datat_VT1_R==T,]
datat_VT1_RP<- data_VT1$Groupe %in% "RP"
datat_VT1_RP<- data_VT1[datat_VT1_RP==T,]

mean_VT1_NR<-mean(datat_VT1_NR$Vidas)
mean_VT1_R<-mean(datat_VT1_R$Vidas)
mean_VT1_RP<-mean(datat_VT1_RP$Vidas)

## VT2----
datat_VT2_NR<- data_VT2$Groupe %in% "NR"
datat_VT2_NR<- data_VT2[datat_VT2_NR==T,]
datat_VT2_R<- data_VT2$Groupe %in% "R"
datat_VT2_R<- data_VT2[datat_VT2_R==T,]
datat_VT2_RP<- data_VT2$Groupe %in% "RP"
datat_VT2_RP<- data_VT2[datat_VT2_RP==T,]

mean_VT2_NR<-mean(datat_VT2_NR$Vidas)
mean_VT2_R<-mean(datat_VT2_R$Vidas)
mean_VT2_RP<-mean(datat_VT2_RP$Vidas)

## VT3----
datat_VT3_NR<- data_VT3$Groupe %in% "NR"
datat_VT3_NR<- data_VT3[datat_VT3_NR==T,]
datat_VT3_R<- data_VT3$Groupe %in% "R"
datat_VT3_R<- data_VT3[datat_VT3_R==T,]
datat_VT3_RP<- data_VT3$Groupe %in% "RP"
datat_VT3_RP<- data_VT3[datat_VT3_RP==T,]

mean_VT3_NR<-mean(datat_VT3_NR$Vidas)
mean_VT3_R<-mean(datat_VT3_R$Vidas)
mean_VT3_RP<-mean(datat_VT3_RP$Vidas)


resume<-t(data.frame(mean_VT1_NR,mean_VT1_R,mean_VT1_RP,mean_VT2_NR,mean_VT2_R,mean_VT2_RP,mean_VT3_NR,mean_VT3_R,mean_VT3_RP))

# write.csv(resume, "resume_mean_ab_covid.csv")

# Stat----
## VT1----
data_VT1 %>% sample_n_by(Groupe, size = 1)
data_VT1 <- data_VT1 %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_VT1)

print(data_VT1 %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_VT1 %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_VT1 %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_VT1 %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_VT1, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("VT1")) +
  ylab("Titre en anti-corps")+
        theme(text = element_text(size = 15))   

## VT2----
data_VT2 %>% sample_n_by(Groupe, size = 1)
data_VT2 <- data_VT2 %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_VT2)

print(data_VT2 %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_VT2 %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_VT2 %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_VT2 %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_VT2, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("VT2")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 

## VT3----
data_VT3 %>% sample_n_by(Groupe, size = 1)
data_VT3 <- data_VT3 %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_VT3)

print(data_VT3 %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_VT3 %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_VT3 %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_VT3 %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_VT3, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("VT3")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 

## VT4----
data_VT4 %>% sample_n_by(Groupe, size = 1)
data_VT4 <- data_VT4 %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_VT4)

print(data_VT4 %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_VT4 %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_VT4 %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_VT4 %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_VT4, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("VT4")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 

## VT5----
data_VT5 %>% sample_n_by(Groupe, size = 1)
data_VT5 <- data_VT5 %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_VT5)

print(data_VT5 %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_VT5 %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_VT5 %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_VT5 %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_VT5, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("VT5")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 

## VT6----
data_VT6 %>% sample_n_by(Groupe, size = 1)
data_VT6 <- data_VT6 %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_VT6)

print(data_VT6 %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_VT6 %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_VT6 %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_VT6 %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_VT6, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("VT6")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 

## VT7----
data_VT7 %>% sample_n_by(Groupe, size = 1)
data_VT7 <- data_VT7 %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_VT7)

print(data_VT7 %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_VT7 %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_VT7 %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_VT7 %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_VT7, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("VT7")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 

## VT8----
data_VT8 %>% sample_n_by(Groupe, size = 1)
data_VT8 <- data_VT8 %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_VT8)

print(data_VT8 %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_VT8 %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_VT8 %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_VT8 %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_VT8, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("VT8")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 


## S1----
data_S1 %>% sample_n_by(Groupe, size = 1)
data_S1 <- data_S1 %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_S1)

print(data_S1 %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_S1 %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_S1 %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_S1 %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_S1, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("S1")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 

## S2----
data_S2 %>% sample_n_by(Groupe, size = 1)
data_S2 <- data_S2 %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_S2)

print(data_S2 %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_S2 %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_S2 %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_S2 %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_S2, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("S2")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 


# Titre en Ac day à 50% du titre

data_final<-
data %>%
  arrange(numero_patient)%>%
  # head(14) %>%
  group_by(numero_patient)%>%
  mutate(max_v = max(Vidas), 
         max_v2=max_v/2,
         test= Vidas>= max_v2)%>%
  dplyr::filter(test == T)%>%
  top_n(1,desc(day))%>%
  select(-max_v,-max_v2, -test)

data_final<-arrange(data_final, Groupe)
data_final<- as.data.frame(data_final)

# data_final ----
data_final %>% sample_n_by(Groupe, size = 1)
data_final <- data_final %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_final)

print(data_final %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_final %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_final %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_final %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values

print(ggboxplot(data_final, x = "Groupe", y = "Vidas") +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("Titre en anticorps à 50% du titre max")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 




print(ggboxplot(data_final, x = "Groupe", y = 'Vidas') +
        geom_point(colour = "gray70")+
        ggtitle("Titre en anticorps à 50% du titre max")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 

# Test de Friedman-----

data_1<- rbind(data_VT1, data_VT2)
data_2 <- rbind(data_VT3, data_VT4,data_VT5,data_VT6, data_VT7,data_VT8)

## TEST data_1 Kruskal + Dunn----
data_1 %>% sample_n_by(Groupe, size = 1)
data_1 <- data_1 %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_1)

print(data_1 %>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_1 %>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_1 %>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_1 %>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_1, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("VT1+VT2")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 


## TEST data_2Kruskal + Dunn----
data_2%>% sample_n_by(Groupe, size = 1)
data_2<- data_2%>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data_1)

print(data_2%>% 
        group_by(Groupe) %>%
        get_summary_stats(`Vidas`, type = "common"))
# Test kruskal
res.kruskal <- data_2%>% kruskal_test(formula(paste0( 'Vidas'," ~ Groupe")))
print(res.kruskal)
# Taille de l’effet
print(data_2%>% kruskal_effsize(formula(paste0('Vidas'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data_2%>% 
  dunn_test(formula(paste0('Vidas'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_2, x = "Groupe", y = 'Vidas') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)
        )+
        geom_point(colour = "gray70")+
        ggtitle("VT3+VT4...")) +
  ylab("Titre en anti-corps")+
  theme(text = element_text(size = 15)) 

## data_1 Friedman----
fried<- data_1[,c(4,6,10,11)]
fried %>%
  group_by(Groupe) %>%
  get_summary_stats(Vidas_N, type = "common")

fried %>%
  group_by(real_time_point) %>%
  get_summary_stats(Vidas_N, type = "common")
# fried<-as.data.frame(fried)
fried$real_time_point<-as.factor(fried$real_time_point)
res.fried <- fried %>% friedman_test(Vidas_N ~ real_time_point |Groupe  )
res.fried

write.xlsx(fried, "result/data_victor_friedman.xlsx")







## data_2 Friedman----


data_2 %>%
  group_by(Groupe) %>%
  get_summary_stats(Vidas_N, type = "common") 
data_2 %>%
  group_by(real_time_point) %>%
  get_summary_stats(Vidas_N, type = "common")

data_2$real_time_point<-as.factor(data_2$real_time_point)
res.fried <- data_2 %>% friedman_test(Vidas_N ~ real_time_point |Groupe  )
res.fried







