# 1 Library---------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)
library(ggplot2)
library(dplyr)
library(cowplot)
library(stringr)

# 2 Import data-----------------------------------------------------------------
rm(list = ls())
data_V4_CD4_8_Agspe <- read_xlsx("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/nano_freq_V4_T_cell_omiq.xlsx",2)
data_V4_T_CD4_8<- read_xlsx("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/nano_freq_V4_T_cell_omiq.xlsx",1)
data_6m <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",8)

data_V4_CD4_8_Agspe_T_CD4_8[,9:62]<- lapply(data_V4_CD4_8_Agspe_T_CD4_8[,9:62], as.numeric) 

# data_V4_CD4_8_Agspe_T_CD4_8 <- data_V4_CD4_8_Agspe_T_CD4_8[data_V4_CD4_8_Agspe_T_CD4_8$numero_patient != 50,] ## supprime le patient 50

data_True_false <- data_V4_CD4_8_Agspe_T_CD4_8$Stimulation %in% "no stim"
data_V4_CD4_8_Agspe_T_CD4_8_no_stim <- data_V4_CD4_8_Agspe_T_CD4_8[data_True_false==T,]
data_V4_CD4_8_Agspe_T_CD4_8_stim <- data_V4_CD4_8_Agspe_T_CD4_8[data_True_false==F,]

colnames(data_6m) <- str_replace_all(colnames(data_6m), pattern = " ", replacement = "_")

# 3 Test stat Kruskal Wallis----------------------------------------------------
## 3.1 data_V4_CD4_8_Agspe------------------------------------------------------

max(colSums(is.na(data_V4_CD4_8_Agspe))) # verif valleur manquante

nom_col <- names(data_V4_CD4_8_Agspe[9:40])

data_V4_CD4_8_Agspe %>% sample_n_by(Responder, size = 1)
data_V4_CD4_8_Agspe <- data_V4_CD4_8_Agspe %>%
  reorder_levels(Responder, order = c("NA", "NR", "R", "RP-", "RP"))

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_CD4_8_Agspe[i])
  variable$Responder<-data_V4_CD4_8_Agspe$Responder
  variable<-na.omit(variable)
  
  # Repartition dans la variable
  # print(variable %>% 
  #         group_by(Responder) %>%
  #         get_summary_stats(`i`, type = "common"))
  
  # Calculs
  res.kruskal <- variable %>% kruskal_test(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.kruskal)
  
  # Taille de l’effet
  res.dunn <- variable %>% kruskal_effsize(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.dunn)
  # Comparaisons par paires test de Dunn
  pwc <- variable %>% 
    dunn_test(formula(paste0("`" , i , "`"," ~ Responder")), p.adjust.method = "bonferroni") 
  # print(pwc)
  
  # Visualisation : Boxplots avec p-values
  pwc <- pwc %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, x = "Responder", y = i) +geom_point(colour = "gray70") +
          stat_pvalue_manual(pwc, hide.ns = TRUE) +
          labs(title = formula(paste0("`" , i , "`"," ~ Responder")),
               subtitle = get_test_label(res.kruskal, detailed = TRUE),
               caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_V4_CD4_8_Agspe.pdf", plot = p_all, width = 20, height = 50, limitsize = F )


## 3.2 data_V4_CD4_8_Agspe sans les NC------------------------------------------
T_F <- data_V4_CD4_8_Agspe$Responder %in% "NA"
data_V4_CD4_8_Agspe <- data_V4_CD4_8_Agspe[T_F == F,]
max(colSums(is.na(data_V4_CD4_8_Agspe))) # verif valleur manquante

nom_col <- names(data_V4_CD4_8_Agspe[9:40])

data_V4_CD4_8_Agspe %>% sample_n_by(Responder, size = 1)
data_V4_CD4_8_Agspe <- data_V4_CD4_8_Agspe %>%
  reorder_levels(Responder, order = c("NA", "NR", "R", "RP-", "RP"))

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_CD4_8_Agspe[i])
  variable$Responder<-data_V4_CD4_8_Agspe$Responder
  variable<-na.omit(variable)
  
  # Repartition dans la variable
  # print(variable %>% 
  #         group_by(Responder) %>%
  #         get_summary_stats(`i`, type = "common"))
  
  # Calculs
  res.kruskal <- variable %>% kruskal_test(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.kruskal)
  
  # Taille de l’effet
  res.dunn <- variable %>% kruskal_effsize(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.dunn)
  # Comparaisons par paires test de Dunn
  pwc <- variable %>% 
    dunn_test(formula(paste0("`" , i , "`"," ~ Responder")), p.adjust.method = "bonferroni") 
  # print(pwc)
  
  # Visualisation : Boxplots avec p-values
  pwc <- pwc %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, x = "Responder", y = i) +geom_point(colour = "gray70") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(title = formula(paste0("`" , i , "`"," ~ Responder")),
         subtitle = get_test_label(res.kruskal, detailed = TRUE),
         caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_V4_CD4_8_Agspe_sans_NA.pdf", plot = p_all, width = 20, height = 50, limitsize = F )


## 3.3 data_V4_CD4_8_Agspe_T_CD4_8_no_stim--------------------------------------
max(colSums(is.na(data_V4_CD4_8_Agspe_T_CD4_8_no_stim))) # verif valleur manquante

nom_col <- names(data_V4_CD4_8_Agspe_T_CD4_8_no_stim[9:62])

data_V4_CD4_8_Agspe_T_CD4_8_no_stim %>% sample_n_by(Responder, size = 1)
data_V4_CD4_8_Agspe_T_CD4_8_no_stim <- data_V4_CD4_8_Agspe_T_CD4_8_no_stim %>%
  reorder_levels(Responder, order = c("NA", "NR", "R", "RP-", "RP"))

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_CD4_8_Agspe_T_CD4_8_no_stim[i])
  variable$Responder<-data_V4_CD4_8_Agspe_T_CD4_8_no_stim$Responder
  variable<-na.omit(variable)
  
  # Repartition dans la variable
  # print(variable %>% 
  #         group_by(Responder) %>%
  #         get_summary_stats(`i`, type = "common"))
  
  # Calculs
  res.kruskal <- variable %>% kruskal_test(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.kruskal)
  
  # Taille de l’effet
  res.dunn <- variable %>% kruskal_effsize(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.dunn)
  # Comparaisons par paires test de Dunn
  pwc <- variable %>% 
    dunn_test(formula(paste0("`" , i , "`"," ~ Responder")), p.adjust.method = "bonferroni") 
  # print(pwc)
  
  # Visualisation : Boxplots avec p-values
  pwc <- pwc %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, x = "Responder", y = i) +geom_point(colour = "gray70") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(title = formula(paste0("`" , i , "`"," ~ Responder")),
         subtitle = get_test_label(res.kruskal, detailed = TRUE),
         caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_data_V4_CD4_8_Agspe_T_CD4_8_no_stim.pdf", plot = p_all, width = 20, height = 75, limitsize = F )


## 3.4 data_V4_CD4_8_Agspe_T_CD4_8_no_stim sans NA------------------------------
T_F <- data_V4_CD4_8_Agspe_T_CD4_8_no_stim$Responder %in% "NA"
data_V4_CD4_8_Agspe_T_CD4_8_no_stim <- data_V4_CD4_8_Agspe_T_CD4_8_no_stim[T_F == F,]
max(colSums(is.na(data_V4_CD4_8_Agspe_T_CD4_8_no_stim))) # verif valleur manquante

nom_col <- names(data_V4_CD4_8_Agspe_T_CD4_8_no_stim[9:62])

data_V4_CD4_8_Agspe_T_CD4_8_no_stim %>% sample_n_by(Responder, size = 1)
data_V4_CD4_8_Agspe_T_CD4_8_no_stim <- data_V4_CD4_8_Agspe_T_CD4_8_no_stim %>%
  reorder_levels(Responder, order = c("NA", "NR", "R", "RP-", "RP"))

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_CD4_8_Agspe_T_CD4_8_no_stim[i])
  variable$Responder<-data_V4_CD4_8_Agspe_T_CD4_8_no_stim$Responder
  variable<-na.omit(variable)
  
  # Repartition dans la variable
  # print(variable %>% 
  #         group_by(Responder) %>%
  #         get_summary_stats(`i`, type = "common"))
  
  # Calculs
  res.kruskal <- variable %>% kruskal_test(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.kruskal)
  
  # Taille de l’effet
  res.dunn <- variable %>% kruskal_effsize(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.dunn)
  # Comparaisons par paires test de Dunn
  pwc <- variable %>% 
    dunn_test(formula(paste0("`" , i , "`"," ~ Responder")), p.adjust.method = "bonferroni") 
  # print(pwc)
  
  # Visualisation : Boxplots avec p-values
  pwc <- pwc %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, x = "Responder", y = i) +geom_point(colour = "gray70") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(title = formula(paste0("`" , i , "`"," ~ Responder")),
         subtitle = get_test_label(res.kruskal, detailed = TRUE),
         caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_data_V4_CD4_8_Agspe_T_CD4_8_no_stim_sans_NA.pdf", plot = p_all, width = 20, height = 75, limitsize = F )


## 3.5 data_V4_CD4_8_Agspe_T_CD4_8_stim-----------------------------------------
max(colSums(is.na(data_V4_CD4_8_Agspe_T_CD4_8_stim))) # verif valleur manquante

nom_col <- names(data_V4_CD4_8_Agspe_T_CD4_8_stim[9:62])

data_V4_CD4_8_Agspe_T_CD4_8_stim %>% sample_n_by(Responder, size = 1)
data_V4_CD4_8_Agspe_T_CD4_8_stim <- data_V4_CD4_8_Agspe_T_CD4_8_stim %>%
  reorder_levels(Responder, order = c("NA", "NR", "R", "RP-", "RP"))

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_CD4_8_Agspe_T_CD4_8_stim[i])
  variable$Responder<-data_V4_CD4_8_Agspe_T_CD4_8_stim$Responder
  variable<-na.omit(variable)
  
  # Repartition dans la variable
  # print(variable %>% 
  #         group_by(Responder) %>%
  #         get_summary_stats(`i`, type = "common"))
  
  # Calculs
  res.kruskal <- variable %>% kruskal_test(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.kruskal)
  
  # Taille de l’effet
  res.dunn <- variable %>% kruskal_effsize(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.dunn)
  # Comparaisons par paires test de Dunn
  pwc <- variable %>% 
    dunn_test(formula(paste0("`" , i , "`"," ~ Responder")), p.adjust.method = "bonferroni") 
  # print(pwc)
  
  # Visualisation : Boxplots avec p-values
  pwc <- pwc %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, x = "Responder", y = i) +geom_point(colour = "gray70") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(title = formula(paste0("`" , i , "`"," ~ Responder")),
         subtitle = get_test_label(res.kruskal, detailed = TRUE),
         caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_V4_CD4_8_Agspe_T_CD4_8_stim.pdf", plot = p_all, width = 20, height = 75, limitsize = F )


## 3.6 data_V4_CD4_8_Agspe_T_CD4_8_stim sans NA------------------------------
T_F <- data_V4_CD4_8_Agspe_T_CD4_8_stim$Responder %in% "NA"
data_V4_CD4_8_Agspe_T_CD4_8_stim <- data_V4_CD4_8_Agspe_T_CD4_8_stim[T_F == F,]
max(colSums(is.na(data_V4_CD4_8_Agspe_T_CD4_8_stim))) # verif valleur manquante

nom_col <- names(data_V4_CD4_8_Agspe_T_CD4_8_stim[9:62])

data_V4_CD4_8_Agspe_T_CD4_8_stim %>% sample_n_by(Responder, size = 1)
data_V4_CD4_8_Agspe_T_CD4_8_stim <- data_V4_CD4_8_Agspe_T_CD4_8_stim %>%
  reorder_levels(Responder, order = c("NA", "NR", "R", "RP-", "RP"))

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_CD4_8_Agspe_T_CD4_8_stim[i])
  variable$Responder<-data_V4_CD4_8_Agspe_T_CD4_8_stim$Responder
  variable<-na.omit(variable)
  
  # Repartition dans la variable
  # print(variable %>% 
  #         group_by(Responder) %>%
  #         get_summary_stats(`i`, type = "common"))
  
  # Calculs
  res.kruskal <- variable %>% kruskal_test(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.kruskal)
  
  # Taille de l’effet
  res.dunn <- variable %>% kruskal_effsize(formula(paste0("`" , i , "`"," ~ Responder")))
  # print(res.dunn)
  # Comparaisons par paires test de Dunn
  pwc <- variable %>% 
    dunn_test(formula(paste0("`" , i , "`"," ~ Responder")), p.adjust.method = "bonferroni") 
  print(pwc)
  
  # Visualisation : Boxplots avec p-values
  pwc <- pwc %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, x = "Responder", y = i) +geom_point(colour = "gray70") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(title = formula(paste0("`" , i , "`"," ~ Responder")),
         subtitle = get_test_label(res.kruskal, detailed = TRUE),
         caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_V4_CD4_8_Agspe_T_CD4_8_stim_sans_NA.pdf", plot = p_all, width = 20, height = 75, limitsize = F )


## 3.7 data_6m-----------------------------------------
max(colSums(is.na(data_6m))) # verif valleur manquante

nom_col <- names(data_6m[14:99])

data_6m %>% sample_n_by(REPONSE, size = 1)
data_6m <- data_6m %>%
  reorder_levels(REPONSE, order = c("A", "NR", "NR-", "R", "RP-", "RP",  "T"))

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_6m[i])
  variable$REPONSE<-data_6m$REPONSE
  variable<-na.omit(variable)
  
  # Repartition dans la variable
  # print(variable %>% 
  #         group_by(REPONSE) %>%
  #         get_summary_stats(`i`, type = "common"))
  
  # Calculs
  res.kruskal <- variable %>% kruskal_test(formula(paste0("`" , i , "`"," ~ REPONSE")))
  print(res.kruskal)
  
  # Taille de l’effet
  res.dunn <- variable %>% kruskal_effsize(formula(paste0("`" , i , "`"," ~ REPONSE")))
  print(res.dunn)
  
  # Comparaisons par paires test de Dunn
  pwc <- variable %>% 
    dunn_test(formula(paste0("`" , i , "`"," ~ REPONSE")), p.adjust.method = "bonferroni") 
  print(pwc)
  
  # Visualisation : Boxplots avec p-values
  pwc <- pwc %>% add_xy_position(x = "REPONSE")
  p_i <- ggboxplot(variable, x = "REPONSE", y = paste0("`" , i , "`")) +
    geom_point(colour = "gray70") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(title = formula(paste0("`" , i , "`"," ~ REPONSE")),
         subtitle = get_test_label(res.kruskal, detailed = TRUE),
         caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_6_months.pdf", plot = p_all, width = 20, height = 100, limitsize = F )










