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

data_V1 <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",7)
data_V4 <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",10)
data_V4_T_CD4_8<- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",11)
data_V4_CD4_8_Agspe <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",12)
data_6m <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",13)

T_F <- data_V1$REPONSE %in% "R"
data <- data_V1[T_F == T,]
T_F <- data_V1$REPONSE %in% "RP"
data_V1 <- rbind(data, data_V1[T_F == T,])

T_F <- data_V4$Responder %in% "R"
data <- data_V4[T_F == T,]
T_F <- data_V4$Responder %in% "RP"
data_V4 <- rbind(data, data_V4[T_F == T,])

data_True_false <- data_V4$Stimulation %in% "no stim"
data_V4_no_stim <- data_V4[data_True_false==T,]
data_V4_stim <- data_V4[data_True_false==F,]

T_F <- data_V4_CD4_8_Agspe$Responder %in% "R"
data <- data_V4_CD4_8_Agspe[T_F == T,]
T_F <- data_V4_CD4_8_Agspe$Responder %in% "RP"
data_V4_CD4_8_Agspe <- rbind(data, data_V4_CD4_8_Agspe[T_F == T,])

T_F <- data_V4_T_CD4_8$Responder %in% "R"
data <- data_V4_T_CD4_8[T_F == T,]
T_F <- data_V4_T_CD4_8$Responder %in% "RP"
data_V4_T_CD4_8 <- rbind(data, data_V4_T_CD4_8[T_F == T,])

T_F <- data_6m$REPONSE %in% "R"
data <- data_6m[T_F == T,]
T_F <- data_6m$REPONSE %in% "RP"
data_6m <- rbind(data, data_6m[T_F == T,])
# T_F <- data_6m$REPONSE %in% "T"
# data_6m <- rbind(data, data_6m[T_F == T,])


data_V4_T_CD4_8[,9:62]<- lapply(data_V4_T_CD4_8[,9:62], as.numeric) 

# data_V4_T_CD4_8 <- data_V4_T_CD4_8[data_V4_T_CD4_8$numero_patient != 50,] ## supprime le patient 50

data_True_false <- data_V4_T_CD4_8$Stimulation %in% "no stim"
data_V4_T_CD4_8_no_stim <- data_V4_T_CD4_8[data_True_false==T,]
data_V4_T_CD4_8_stim <- data_V4_T_CD4_8[data_True_false==F,]

# colnames(data_6m) <- str_replace_all(colnames(data_6m), pattern = " ", replacement = "_")


# 3 Test stat Kruskal Wallis----------------------------------------------------
## 3.0 data_V4_no_stim supervise---------------------------------------------------------- 
max(colSums(is.na(data_V4_no_stim))) # verif valleur manquante

nom_col <- names(data_V4_no_stim[8:62])

data_V4_no_stim %>% sample_n_by(Responder, size = 1)
data_V4_no_stim <- data_V4_no_stim %>%
  reorder_levels(Responder, order = c("R", "RP"))

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_no_stim[i])
  variable$Responder<-data_V4_no_stim$Responder
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
  p_i <- ggboxplot(variable, x = "Responder", y = paste0("`" , i , "`")) +geom_point(colour = "gray70") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(title = formula(paste0("`" , i , "`"," ~ Responder")),
         subtitle = get_test_label(res.kruskal, detailed = TRUE),
         caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_V4_no_stim_R_RP_sup.pdf", plot = p_all, width = 20, height = 75, limitsize = F )


## 3.0 data_V4_stim supervise-------------------------------------------------------------
max(colSums(is.na(data_V4_stim))) # verif valleur manquante

nom_col <- names(data_V4_stim[8:62])

data_V4_stim %>% sample_n_by(Responder, size = 1)
data_V4_stim <- data_V4_stim %>%
  reorder_levels(Responder, order = c("R", "RP"))

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_stim[i])
  variable$Responder<-data_V4_stim$Responder
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
  p_i <- ggboxplot(variable, x = "Responder", y = paste0("`" , i , "`")) +geom_point(colour = "gray70") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(title = formula(paste0("`" , i , "`"," ~ Responder")),
         subtitle = get_test_label(res.kruskal, detailed = TRUE),
         caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_V4_stim_R_RP_sup.pdf", plot = p_all, width = 20, height = 75, limitsize = F )



## 3.1 data_V4_CD4_8_Agspe non supervise------------------------------------------------------
max(colSums(is.na(data_V4_CD4_8_Agspe))) # verif valleur manquante

nom_col <- names(data_V4_CD4_8_Agspe[9:40])

data_V4_CD4_8_Agspe %>% sample_n_by(Responder, size = 1)
data_V4_CD4_8_Agspe <- data_V4_CD4_8_Agspe %>%
  reorder_levels(Responder, order = c("R", "RP"))

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
ggsave(filename = "result/kruskal_wallis_V4_CD4_8_Agspe_R_RP_non_sup.pdf", plot = p_all, width = 20, height = 40, limitsize = F )

## 3.2 data_V4_T_CD4_8_no_stim non supervise------------------------------------------------------
max(colSums(is.na(data_V4_T_CD4_8_no_stim))) # verif valleur manquante

nom_col <- names(data_V4_T_CD4_8_no_stim[9:62])

data_V4_T_CD4_8_no_stim %>% sample_n_by(Responder, size = 1)
data_V4_T_CD4_8_no_stim <- data_V4_T_CD4_8_no_stim %>%
  reorder_levels(Responder, order = c("R", "RP"))

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_T_CD4_8_no_stim[i])
  variable$Responder<-data_V4_T_CD4_8_no_stim$Responder
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
ggsave(filename = "result/kruskal_wallis_V4_T_CD4_8_no_stim_R_RP_non_sup.pdf", plot = p_all, width = 20, height = 70, limitsize = F )


## 3.3 data_V4_T_CD4_8_stim non supervise------------------------------------------------------
max(colSums(is.na(data_V4_T_CD4_8_stim))) # verif valleur manquante

nom_col <- names(data_V4_T_CD4_8_stim[9:62])

data_V4_T_CD4_8_stim %>% sample_n_by(Responder, size = 1)
data_V4_T_CD4_8_stim <- data_V4_T_CD4_8_stim %>%
  reorder_levels(Responder, order = c("R", "RP"))

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_T_CD4_8_stim[i])
  variable$Responder<-data_V4_T_CD4_8_stim$Responder
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
ggsave(filename = "result/kruskal_wallis_V4_T_CD4_8_stim_R_RP_non_sup.pdf", plot = p_all, width = 20, height = 70, limitsize = F )


## 3.4 data_6m supervise--------------------------------------------------------
max(colSums(is.na(data_6m))) # verif valleur manquante

nom_col <- names(data_6m[14:99])

data_6m %>% sample_n_by(REPONSE, size = 1)
data_6m <- data_6m %>%
  reorder_levels(REPONSE, order = c( "R", "RP"))  ## "T",

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
  # print(res.kruskal)
  
  # Taille de l’effet
  res.dunn <- variable %>% kruskal_effsize(formula(paste0("`" , i , "`"," ~ REPONSE")))
  # print(res.dunn)
  # Comparaisons par paires test de Dunn
  pwc <- variable %>% 
    dunn_test(formula(paste0("`" , i , "`"," ~ REPONSE")), p.adjust.method = "bonferroni") 
  # print(pwc)
  
  # Visualisation : Boxplots avec p-values
  pwc <- pwc %>% add_xy_position(x = "REPONSE")
  p_i <- ggboxplot(variable, x = "REPONSE", y = paste0("`" , i , "`")) +geom_point(colour = "gray70") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(title = formula(paste0("`" , i , "`"," ~ REPONSE")),
         subtitle = get_test_label(res.kruskal, detailed = TRUE),
         caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_6_months_R_RP_sup.pdf", plot = p_all, width = 25, height = 140, limitsize = F )



## 3.5 data_V1 non supervisé----------------------------------------------------
max(colSums(is.na(data_V1))) # verif valleur manquante

nom_col <- names(data_V1[9:41])

data_V1 %>% sample_n_by(REPONSE, size = 1)
data_V1 <- data_V1 %>% 
  reorder_levels(REPONSE, order = c( "R", "RP"))  ## "T",

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V1[i])
  variable$REPONSE<-data_V1$REPONSE
  variable<-na.omit(variable)
  
  # Repartition dans la variable
  # print(variable %>% 
  #         group_by(REPONSE) %>%
  #         get_summary_stats(`i`, type = "common"))
  
  # Calculs
  res.kruskal <- variable %>% kruskal_test(formula(paste0("`" , i , "`"," ~ REPONSE")))
  # print(res.kruskal)
  
  # Taille de l’effet
  res.dunn <- variable %>% kruskal_effsize(formula(paste0("`" , i , "`"," ~ REPONSE")))
  # print(res.dunn)
  # Comparaisons par paires test de Dunn
  pwc <- variable %>% 
    dunn_test(formula(paste0("`" , i , "`"," ~ REPONSE")), p.adjust.method = "bonferroni") 
  # print(pwc)
  
  # Visualisation : Boxplots avec p-values
  pwc <- pwc %>% add_xy_position(x = "REPONSE")
  p_i <- ggboxplot(variable, x = "REPONSE", y = paste0("`" , i , "`")) +geom_point(colour = "gray70") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(title = paste0(i),
         subtitle = get_test_label(res.kruskal, detailed = TRUE),
         caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_V1_non_supervise_R_RP_non_sup.pdf", plot = p_all, width = 20, height = 50, limitsize = F )



## 3.6 data_V1 supervisé--------------------------------------------------------
data_V1 <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",6)
max(colSums(is.na(data_V1))) # verif valleur manquante

nom_col <- names(data_V1[11:50])

data_V1 %>% sample_n_by(REPONSE, size = 1)
data_V1 <- data_V1 %>%
  reorder_levels(REPONSE, order = c( "R", "RP"))  ## "T",

### Test Kruskal Wallis + Dunn ----
list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V1[i])
  variable$REPONSE<-data_V1$REPONSE
  variable<-na.omit(variable)
  
  # Repartition dans la variable
  # print(variable %>% 
  #         group_by(REPONSE) %>%
  #         get_summary_stats(`i`, type = "common"))
  
  # Calculs
  res.kruskal <- variable %>% kruskal_test(formula(paste0("`" , i , "`"," ~ REPONSE")))
  # print(res.kruskal)
  
  # Taille de l’effet
  res.dunn <- variable %>% kruskal_effsize(formula(paste0("`" , i , "`"," ~ REPONSE")))
  # print(res.dunn)
  # Comparaisons par paires test de Dunn
  pwc <- variable %>% 
    dunn_test(formula(paste0("`" , i , "`"," ~ REPONSE")), p.adjust.method = "bonferroni") 
  # print(pwc)
  
  # Visualisation : Boxplots avec p-values
  pwc <- pwc %>% add_xy_position(x = "REPONSE")
  p_i <- ggboxplot(variable, x = "REPONSE", y = paste0("`" , i , "`")) +geom_point(colour = "gray70") +
    stat_pvalue_manual(pwc, hide.ns = TRUE) +
    labs(title = paste0(i),
         subtitle = get_test_label(res.kruskal, detailed = TRUE),
         caption = get_pwc_label(pwc))
  
  list_plot[[i]] <- p_i
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/kruskal_wallis_V1_supervise_R_RP_sup.pdf", plot = p_all, width = 20, height = 50, limitsize = F )
