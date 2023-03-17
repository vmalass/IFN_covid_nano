# 1 Library---------------------------------------------------------------------
library(readxl)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(cowplot)


# library(ggplot2)
# library(dplyr)
# library(stringr)




# 2 Import data-----------------------------------------------------------------
rm(list = ls())

set.seed(123)

data_V1_sup <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",6)
data_V1_non_sup <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",7)
data_V4 <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",10)
data_V4_T_CD4_8<- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",11)
data_V4_CD4_8_Agspe <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",12)
data_V4_T <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",13)
data_V4_T_sup <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",14)
data_6m <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",14)

## 2.1 Formatage data R RP------------------------------------------------------
### data_V1_sup ###
T_F <- data_V1_sup$REPONSE %in% "R"
data <- data_V1_sup[T_F == T,]
T_F <- data_V1_sup$REPONSE %in% "RP"
data_V1_sup<- rbind(data, data_V1_sup[T_F == T,])

### data_V1_non_sup ###
T_F <- data_V1_non_sup$REPONSE %in% "R"
data <- data_V1_non_sup[T_F == T,]
T_F <- data_V1_non_sup$REPONSE %in% "RP"
data_V1_non_sup<- rbind(data, data_V1_non_sup[T_F == T,])

### data_V4 ###
T_F <- data_V4$Responder %in% "R"
data <- data_V4[T_F == T,]
T_F <- data_V4$Responder %in% "RP"
data_V4 <- rbind(data, data_V4[T_F == T,])

data_True_false <- data_V4$Stimulation %in% "no stim"
data_V4_no_stim <- data_V4[data_True_false==T,]
data_V4_stim <- data_V4[data_True_false==F,]

### data_V4_CD4_8_Agspe ###
T_F <- data_V4_CD4_8_Agspe$Responder %in% "R"
data <- data_V4_CD4_8_Agspe[T_F == T,]
T_F <- data_V4_CD4_8_Agspe$Responder %in% "RP"
data_V4_CD4_8_Agspe <- rbind(data, data_V4_CD4_8_Agspe[T_F == T,])

### data_V4_T_CD4_8 ###
T_F <- data_V4_T_CD4_8$Responder %in% "R"
data <- data_V4_T_CD4_8[T_F == T,]
T_F <- data_V4_T_CD4_8$Responder %in% "RP"
data_V4_T_CD4_8 <- rbind(data, data_V4_T_CD4_8[T_F == T,])

data_V4_T_CD4_8[,9:62]<- lapply(data_V4_T_CD4_8[,9:62], as.numeric) 

data_True_false <- data_V4_T_CD4_8$Stimulation %in% "no stim"
data_V4_T_CD4_8_no_stim <- data_V4_T_CD4_8[data_True_false==T,]
data_V4_T_CD4_8_stim <- data_V4_T_CD4_8[data_True_false==F,]

### data_V4_T ###
T_F <- data_V4_T$Responder %in% "R"
data <- data_V4_T[T_F == T,]
T_F <- data_V4_T$Responder %in% "RP"
data_V4_T <- rbind(data, data_V4_T[T_F == T,])

### data_V4_T_sup ###
T_F <- data_V4_T_sup$Responder %in% "R"
data <- data_V4_T_sup[T_F == T,]
T_F <- data_V4_T_sup$Responder %in% "RP"
data_V4_T_sup <- rbind(data, data_V4_T_sup[T_F == T,])

### data_6m ###
T_F <- data_6m$REPONSE %in% "R"
data <- data_6m[T_F == T,]
T_F <- data_6m$REPONSE %in% "RP"
data_6m <- rbind(data, data_6m[T_F == T,])



# 3 Test Wilcoxon---------------------------------------------------------------
## 3.1 data_V1_sup--------------------------------------------------------------
max(colSums(is.na(data_V1_sup))) # verif valleur manquante

nom_col <- names(data_V1_sup[8:56])

# data_V1_sup%>% sample_n_by(REPONSE, size = 5)
data_V1_sup<- data_V1_sup%>% 
  reorder_levels(REPONSE, order = c( "R", "RP"))  ## "T",

list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V1_sup[i])
  variable$REPONSE<-data_V1_sup$REPONSE
  variable<-na.omit(variable)
  
  # Test stat
  res.wilcox <- variable %>%
    wilcox_test(formula(paste0("`" , i , "`"," ~ REPONSE"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "REPONSE")
  p_i <- ggboxplot(variable, 
                   x = "REPONSE", 
                   y = paste0("`" , i , "`"), 
                   ylab = i, 
                   xlab = "Groups", 
                   add = "point") + 
    stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
    labs(subtitle = get_test_label(res.wilcox, detailed = T))
  
  # Liste with plot
  list_plot[[i]] <- p_i
}

# plot_grid(plotlist = list_plot[1:6], ncol = 2) 

p_all = plot_grid(plotlist = list_plot, ncol = 3)  # pour orga ne nombre de plot sur la page
ggsave(filename = "result/wilcox_test_V1_sup_R_RP.pdf", plot = p_all, width = 20, height = 75, limitsize = F )




## 3.2 data_V1_non_sup----------------------------------------------------------
max(colSums(is.na(data_V1_non_sup))) # verif valleur manquante

nom_col <- names(data_V1_non_sup[9:41])

# data_V1_non_sup%>% sample_n_by(REPONSE, size = 5)
data_V1_non_sup<- data_V1_non_sup%>% 
  reorder_levels(REPONSE, order = c( "R", "RP"))  ## "T",

list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V1_non_sup[i])
  variable$REPONSE<-data_V1_non_sup$REPONSE
  variable<-na.omit(variable)
  
  # Test stat
  res.wilcox <- variable %>%
    wilcox_test(formula(paste0("`" , i , "`"," ~ REPONSE"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "REPONSE")
  p_i <- ggboxplot(variable, 
                   x = "REPONSE", 
                   y = paste0("`" , i , "`"), 
                   ylab = i, 
                   xlab = "Groups", 
                   add = "point") + 
    stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
    labs(subtitle = get_test_label(res.wilcox, detailed = T))

  # Liste with plot
  list_plot[[i]] <- p_i
}

# plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/wilcox_test_V1_non_sup_R_RP.pdf", plot = p_all, width = 20, height = 50, limitsize = F )



## 3.3 data_V4_no_stim sup------------------------------------------------------
max(colSums(is.na(data_V4_no_stim))) # verif valleur manquante

nom_col <- names(data_V4_no_stim[8:62])

# data_V4_no_stim %>% sample_n_by(Responder, size = 5)
data_V4_no_stim <- data_V4_no_stim %>% 
  reorder_levels(Responder, order = c( "R", "RP"))  ## "T",

list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_no_stim[i])
  variable$Responder<-data_V4_no_stim$Responder
  variable<-na.omit(variable)
  
  # Test stat
  res.wilcox <- variable %>%
    wilcox_test(formula(paste0("`" , i , "`"," ~ Responder"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, 
                   x = "Responder", 
                   y = paste0("`" , i , "`"), 
                   ylab = i, 
                   xlab = "Groups", 
                   add = "point") + 
    stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
    labs(subtitle = get_test_label(res.wilcox, detailed = T))
  
  # Liste with plot
  list_plot[[i]] <- p_i
}

# plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/wilcox_test_V4_non_stim_sup_R_RP.pdf", plot = p_all, width = 20, height = 70, limitsize = F )



## 3.4 data_V4_stim sup---------------------------------------------------------
max(colSums(is.na(data_V4_stim))) # verif valleur manquante

nom_col <- names(data_V4_stim[8:62])

# data_V4_stim %>% sample_n_by(Responder, size = 5)
data_V4_stim <- data_V4_stim %>% 
  reorder_levels(Responder, order = c( "R", "RP"))  ## "T",

list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_stim[i])
  variable$Responder<-data_V4_stim$Responder
  variable<-na.omit(variable)
  
  # Test stat
  res.wilcox <- variable %>%
    wilcox_test(formula(paste0("`" , i , "`"," ~ Responder"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, 
                   x = "Responder", 
                   y = paste0("`" , i , "`"), 
                   ylab = i, 
                   xlab = "Groups", 
                   add = "point") + 
    stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
    labs(subtitle = get_test_label(res.wilcox, detailed = T))
  
  # Liste with plot
  list_plot[[i]] <- p_i
}

# plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/wilcox_test_V4_stim_sup_R_RP.pdf", plot = p_all, width = 20, height = 70, limitsize = F )



## 3.5 data_V4_CD4_8_Agspe stim non sup---------------------------------------------
max(colSums(is.na(data_V4_CD4_8_Agspe))) # verif valleur manquante

nom_col <- names(data_V4_CD4_8_Agspe[9:40])

# data_V4_CD4_8_Agspe %>% sample_n_by(Responder, size = 5)
data_V4_CD4_8_Agspe <- data_V4_CD4_8_Agspe %>% 
  reorder_levels(Responder, order = c( "R", "RP"))  ## "T",

list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_CD4_8_Agspe[i])
  variable$Responder<-data_V4_CD4_8_Agspe$Responder
  variable<-na.omit(variable)
  
  # Test stat
  res.wilcox <- variable %>%
    wilcox_test(formula(paste0("`" , i , "`"," ~ Responder"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, 
                   x = "Responder", 
                   y = paste0("`" , i , "`"), 
                   ylab = i, 
                   xlab = "Groups", 
                   add = "point") + 
    stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
    labs(subtitle = get_test_label(res.wilcox, detailed = T))
  
  # Liste with plot
  list_plot[[i]] <- p_i
}

# plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/wilcox_test_V4_CD4_8_Ag_spe_stim_non_sup_R_RP.pdf", plot = p_all, width = 20, height = 50, limitsize = F )



## 3.6 data_V4_T_CD4_8_no_stim non sup------------------------------------------
max(colSums(is.na(data_V4_T_CD4_8_no_stim))) # verif valleur manquante

nom_col <- names(data_V4_T_CD4_8_no_stim[9:62])

# data_V4_T_CD4_8_no_stim %>% sample_n_by(Responder, size = 5)
data_V4_T_CD4_8_no_stim <- data_V4_T_CD4_8_no_stim %>% 
  reorder_levels(Responder, order = c( "R", "RP"))  ## "T",

list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_T_CD4_8_no_stim[i])
  variable$Responder<-data_V4_T_CD4_8_no_stim$Responder
  variable<-na.omit(variable)
  
  # Test stat
  res.wilcox <- variable %>%
    wilcox_test(formula(paste0("`" , i , "`"," ~ Responder"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, 
                   x = "Responder", 
                   y = paste0("`" , i , "`"), 
                   ylab = i, 
                   xlab = "Groups", 
                   add = "point") + 
    stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
    labs(subtitle = get_test_label(res.wilcox, detailed = T))
  
  # Liste with plot
  list_plot[[i]] <- p_i
}

# plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/wilcox_test_V4_T_CD4_8_no_stim_non_sup_R_RP.pdf", plot = p_all, width = 20, height = 70, limitsize = F )



## 3.7 data_V4_T_CD4_8_stim non sup---------------------------------------------
max(colSums(is.na(data_V4_T_CD4_8_stim))) # verif valleur manquante

nom_col <- names(data_V4_T_CD4_8_stim[9:62])

# data_V4_T_CD4_8_stim %>% sample_n_by(Responder, size = 5)
data_V4_T_CD4_8_stim <- data_V4_T_CD4_8_stim %>% 
  reorder_levels(Responder, order = c( "R", "RP"))  ## "T",

list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_T_CD4_8_stim[i])
  variable$Responder<-data_V4_T_CD4_8_stim$Responder
  variable<-na.omit(variable)
  
  # Test stat
  res.wilcox <- variable %>%
    wilcox_test(formula(paste0("`" , i , "`"," ~ Responder"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, 
                   x = "Responder", 
                   y = paste0("`" , i , "`"), 
                   ylab = i, 
                   xlab = "Groups", 
                   add = "point") + 
    stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
    labs(subtitle = get_test_label(res.wilcox, detailed = T))
  
  # Liste with plot
  list_plot[[i]] <- p_i
}

# plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/wilcox_test_V4_T_CD4_8_stim_non_sup_R_RP.pdf", plot = p_all, width = 20, height = 70, limitsize = F )



## 3.8 data_V4_T non sup------------------------------------------------------------
max(colSums(is.na(data_V4_T))) # verif valleur manquante

nom_col <- names(data_V4_T[7:43])

# data_V4_T %>% sample_n_by(Responder, size = 5)
data_V4_T <- data_V4_T %>% 
  reorder_levels(Responder, order = c( "R", "RP"))  ## "T",

list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_T[i])
  variable$Responder<-data_V4_T$Responder
  variable<-na.omit(variable)
  
  # Test stat
  res.wilcox <- variable %>%
    wilcox_test(formula(paste0("`" , i , "`"," ~ Responder"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, 
                   x = "Responder", 
                   y = paste0("`" , i , "`"), 
                   ylab = i, 
                   xlab = "Groups", 
                   add = "point") + 
    stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
    labs(subtitle = get_test_label(res.wilcox, detailed = T))
  
  # Liste with plot
  list_plot[[i]] <- p_i
}

# plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/wilcox_test_V4_T_stim_non_sup_R_RP.pdf", plot = p_all, width = 20, height = 50, limitsize = F )




## 3.9 data_V4_T_sup------------------------------------------------------------
max(colSums(is.na(data_V4_T_sup))) # verif valleur manquante

nom_col <- names(data_V4_T_sup[7:68])

# data_V4_T_sup %>% sample_n_by(Responder, size = 5)
data_V4_T_sup <- data_V4_T_sup %>% 
  reorder_levels(Responder, order = c( "R", "RP"))  ## "T",

list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V4_T_sup[i])
  variable$Responder<-data_V4_T_sup$Responder
  variable<-na.omit(variable)
  
  # Test stat
  res.wilcox <- variable %>%
    wilcox_test(formula(paste0("`" , i , "`"," ~ Responder"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "Responder")
  p_i <- ggboxplot(variable, 
                   x = "Responder", 
                   y = paste0("`" , i , "`"), 
                   ylab = i, 
                   xlab = "Groups", 
                   add = "point") + 
    stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
    labs(subtitle = get_test_label(res.wilcox, detailed = T))
  
  # Liste with plot
  list_plot[[i]] <- p_i
}

# plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/wilcox_test_V4_T_stim_sup_R_RP.pdf", plot = p_all, width = 20, height = 100, limitsize = F )




## 3.10 data_6m sup--------------------------------------------------------------
max(colSums(is.na(data_6m))) # verif valleur manquante

nom_col <- names(data_6m[14:99])

# data_6m %>% sample_n_by(REPONSE, size = 5)
data_6m <- data_6m %>% 
  reorder_levels(REPONSE, order = c( "R", "RP"))  ## "T",

list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_6m[i])
  variable$REPONSE<-data_6m$REPONSE
  variable<-na.omit(variable)
  
  # Test stat
  res.wilcox <- variable %>%
    wilcox_test(formula(paste0("`" , i , "`"," ~ REPONSE"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "REPONSE")
  p_i <- ggboxplot(variable, 
                   x = "REPONSE", 
                   y = paste0("`" , i , "`"), 
                   ylab = i, 
                   xlab = "Groups", 
                   add = "point") + 
    stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
    labs(subtitle = get_test_label(res.wilcox, detailed = T))
  
  # Liste with plot
  list_plot[[i]] <- p_i
}

# plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/wilcox_test_6months_sup_R_RP.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

















