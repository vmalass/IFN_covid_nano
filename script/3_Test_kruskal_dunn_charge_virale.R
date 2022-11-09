# Test stat kruskal Willis puis DUNN sur les groupes en fct de la charge virale

library(openxlsx)
library(dplyr)
library(ggpubr)
library(rstatix)

rm(list = ls())

data<- read.xlsx("data/charge_virale_groupe.xlsx")

data %>% sample_n_by(Groupe, size = 1)
data <- data %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data)

print(data %>% 
        group_by(Groupe) %>%
        get_summary_stats(`charge_virale`, type = "common"))

# Test kruskal
res.kruskal <- data %>% kruskal_test(formula(paste0( 'charge_virale'," ~ Groupe")))
print(res.kruskal)

# Taille de l’effet
print(data %>% kruskal_effsize(formula(paste0('charge_virale'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data %>% 
  dunn_test(formula(paste0('charge_virale'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data, x = "Groupe", y = 'charge_virale') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc))) 

### Avec data log 10 / 1 000 000 copie
rm(list = ls())

data<- read.xlsx("data/charge_virale_groupe.xlsx")

data %>% sample_n_by(Groupe, size = 1)
data <- data %>%
  reorder_levels(Groupe, order = c( "NR","R", "RP"))
summary(data)

print(data %>% 
        group_by(Groupe) %>%
        get_summary_stats(`charge_virale_log10`, type = "common"))

# Test kruskal
res.kruskal <- data %>% kruskal_test(formula(paste0( 'charge_virale_log10'," ~ Groupe")))
print(res.kruskal)

# Taille de l’effet
print(data %>% kruskal_effsize(formula(paste0('charge_virale_log10'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data %>% 
  dunn_test(formula(paste0('charge_virale_log10'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data, x = "Groupe", y = 'charge_virale_log10') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)))



### Avec data R RP +3J after symptome
rm(list = ls())

data<- read.xlsx("data/charge_virale_groupe.xlsx")
data<- data[data$jours_prelevement != "1",] #supprime tout les jours 1
data<- data[data$jours_prelevement != "2",] #supprime tout les jours 1
data<- data[data$Groupe != "NR",] #supprime tout les jours 1


data %>% sample_n_by(Groupe, size = 1)
data <- data %>%
  reorder_levels(Groupe, order = c("R", "RP"))
summary(data)

print(data %>% 
        group_by(Groupe) %>%
        get_summary_stats(`charge_virale_log10`, type = "common"))

# Test kruskal
res.kruskal <- data %>% kruskal_test(formula(paste0( 'charge_virale_log10'," ~ Groupe")))
print(res.kruskal)

# Taille de l’effet
print(data %>% kruskal_effsize(formula(paste0('charge_virale_log10'," ~ Groupe"))))

# Comparaisons par paires test de Dunn
pwc <- data %>% 
  dunn_test(formula(paste0('charge_virale_log10'," ~ Groupe")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Groupe")
print(ggboxplot(data, x = "Groupe", y = 'charge_virale_log10') +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)))
