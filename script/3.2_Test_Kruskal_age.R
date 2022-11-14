# Test stat kruskal Willis puis DUNN sur les Reponses en fct de l'age

library(openxlsx)
library(dplyr)
library(ggpubr)
library(rstatix)
library(ggplot2)

rm(list = ls())

data<- read.xlsx("data/unique_patient.xlsx")
data<- data[data$Reponse != "NC",] #supprime tout les NC

data %>% sample_n_by(Reponse, size = 1)
data <- data %>%
  reorder_levels(Reponse, order = c( "NR","R", "RP"))
summary(data)

print(data %>% 
        group_by(Reponse) %>%
        get_summary_stats(`Age`, type = "common"))

# Test kruskal
res.kruskal <- data %>% kruskal_test(formula(paste0( 'Age'," ~ Reponse")))
print(res.kruskal)

# Taille de l’effet
print(data %>% kruskal_effsize(formula(paste0('Age'," ~ Reponse"))))

# Comparaisons par paires test de Dunn
pwc <- data %>% 
  dunn_test(formula(paste0('Age'," ~ Reponse")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Reponse")
print(ggboxplot(data, x = "Reponse", y = 'Age') +
        geom_point()+
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(
          subtitle = get_test_label(res.kruskal, detailed = TRUE),
          caption = get_pwc_label(pwc)))

