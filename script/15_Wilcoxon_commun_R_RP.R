Stat Wilcoxon sur R et RP commun soit 12R et 5 RP

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

data_V1_non_sup <- read_xlsx("data/Freq_V1_V4_R_RP.xlsx", 1)
data_V1_sup <- read_xlsx("data/Freq_V1_V4_R_RP.xlsx", 3)

data_V1_non_sup[,6:42] <- NULL
data_V1_sup[,57:111] <- NULL
data_V1_sup[,8:9] <- NULL



# 3 Test Wilcoxon---------------------------------------------------------------
## 3.1 data_V1_sup--------------------------------------------------------------
max(colSums(is.na(data_V1_sup))) # verif valleur manquante

nom_col <- names(data_V1_sup[8:54])

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
ggsave(filename = "result/wilcox_test_V1_sup_R_RP_commun.pdf", plot = p_all, width = 20, height = 75, limitsize = F )




## 3.2 data_V1_non_sup----------------------------------------------------------
max(colSums(is.na(data_V1_non_sup))) # verif valleur manquante

nom_col <- names(data_V1_non_sup[6:38])

# data_V1_non_sup%>% sample_n_by(Reponse, size = 5)
data_V1_non_sup<- data_V1_non_sup%>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

list_plot <- list()
for (i in nom_col) {
  # suppression des NA dans la variable
  variable<-as.data.frame(data_V1_non_sup[i])
  variable$Reponse<-data_V1_non_sup$Reponse
  variable<-na.omit(variable)
  
  # Test stat
  res.wilcox <- variable %>%
    wilcox_test(formula(paste0("`" , i , "`"," ~ Reponse"))) %>%
    add_significance()
  
  # Visualization
  res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
  p_i <- ggboxplot(variable, 
                   x = "Reponse", 
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
ggsave(filename = "result/wilcox_test_V1_non_sup_R_RP_commun.pdf", plot = p_all, width = 20, height = 50, limitsize = F )






