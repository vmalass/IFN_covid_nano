# 1- Library----
library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)
library(ggplot2)
library(dplyr)

# 2- import data----
setwd("~/Documents/JM/NanoString/Data_brut")
rm(list = ls())
data_V4 <- read_xlsx("Covid_IFN-avec_data_V5.xlsx",10)

# 3- Nettoyage----
data_V4<- data_V4[data_V4$Stimulation!="no stim",] #supprime tout les no stim
data_V4<- data_V4[data_V4$response != "NC",] #supprime tout les no stim
data_V4<- data_V4[data_V4$response != "NR",] #supprime tout les no stim
max(colSums(is.na(data_V4))) # verif valleur manquante

nom_col <- names(data_V4[5:38])

data_V4 %>% sample_n_by(response, size = 1)
data_V4 <- data_V4 %>%
  reorder_levels(response, order = c("NC", "NR", "R", "RP"))
summary(data_V4)

# 4-Test Kruskal Wallis + Dunn ----
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4[i])
  variable$response<-data_V4$response
  variable<-na.omit(variable)
  
  # Repartition dans la variable
  print(variable %>% 
          group_by(response) %>%
          get_summary_stats(`i`, type = "common"))
  
  # Calculs
  res.kruskal <- variable %>% kruskal_test(formula(paste0("`" , i , "`"," ~ response")))
  print(res.kruskal)
  
  # Taille de l’effet
  print(variable %>% kruskal_effsize(formula(paste0("`" , i , "`"," ~ response"))))
  
  # Comparaisons par paires test de Dunn
  pwc <- variable %>% 
    dunn_test(formula(paste0("`" , i , "`"," ~ response")), p.adjust.method = "bonferroni") 
  print(pwc)
  
  # Visualisation : Boxplots avec p-values
  pwc <- pwc %>% add_xy_position(x = "response")
  print(ggboxplot(variable, x = "response", y = i) +geom_point(colour = "gray70") +
          stat_pvalue_manual(pwc, hide.ns = TRUE) +
          labs(title = formula(paste0("`" , i , "`"," ~ response")),
            subtitle = get_test_label(res.kruskal, detailed = TRUE),
            caption = get_pwc_label(pwc)
          ) ) 
}
