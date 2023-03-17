library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(rstatix)
library(cowplot)

# import datat----
rm(list = ls())
data <- read.xlsx("data/titre_virale_IFI27_Ac.xlsx", sheet = 4)
data_months <- read.xlsx("data/titre_virale_IFI27_Ac.xlsx", sheet = 5)
data$numero_patient<- as.factor(data$numero_patient)

# creation des sous matrice pour les groupes-----
new <- data$Groupe %in% "NR"
groupe_NR <- data[new==T,]
new <- data$Groupe %in% "NR-"
groupe_NRm <- data[new==T,]
new <- data$Groupe %in% "R"
groupe_R <- data[new==T,]
new <- data$Groupe %in% "RP-"
groupe_RPm <- data[new==T,]
new <- data$Groupe %in% "RP"
groupe_RP <- data[new==T,]

new <- data_months$months %in% "S1"
data_S1 <- data_months[new==T,]
new <- data_months$months %in% "S2"
data_S2 <- data_months[new==T,]

numgp_NR <- groupe_NR$numero_patient
numgp_R <- groupe_R$numero_patient
numgp_RP <- groupe_RP$numero_patient

data_NR_R_RP <- rbind(groupe_R, groupe_RP)
data_NR_R_RP <- rbind(data_NR_R_RP, groupe_NR)

a <- ggplot(data = data_NR_R_RP, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "Groupe ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() +
  scale_color_manual(breaks = c("NR", "R", "RP"),
                     values = c("#7570BE","cornflowerblue","brown3")) +
  ggtitle("Titre en anti-corps")+
  scale_x_continuous(name = "jours post infection",
                     limits = c(0, 60))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15),
        legend.position = "none") 

b <- ggplot(data = data_NR_R_RP, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "Groupe ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() +
  scale_color_manual(breaks = c("NR", "R", "RP"),
                     values = c("#7570BE","cornflowerblue","brown3")) +
  ggtitle("Titre en anti-corps")+
  scale_x_continuous(name = "jours post infection",
                     limits = c(0, 20))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15)) 

plot_grid(a, b, labels=c("A", "B"), ncol = 2, nrow = 1)  # Plot avec 3 figures

## visualisation NR/R/RP--------------------------------------------------------
a <- ggplot(data = groupe_NR, 
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

b <- ggplot(data = groupe_R, 
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

c <- ggplot(data = groupe_RP, 
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


plot_grid(a, b, c, labels=c("A", "B", "C"), ncol = 3, nrow = 1)  # Plot avec 3 figures


# Bornage 0 à 20 jours NR/R/RP--------------------------------------------------
a<- ggplot(data = groupe_NR, 
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

b <- ggplot(data = groupe_R, 
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

c <- ggplot(data = groupe_RP, 
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


plot_grid(a, b, c, labels=c("A", "B", "C"), ncol = 3, nrow = 1)  # Plot avec 3 figures









## visualisation NR/NR-/R/RP-/RP------------------------------------------------
a <- ggplot(data = groupe_NR, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "numero_patient ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("NR")+
  scale_x_continuous(name = "jours post infection",
                     limits = c(0, 75))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15),
        legend.position = "none")   

b <- ggplot(data = groupe_NRm, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "numero_patient ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("NR-")+
  scale_x_continuous(name = "jours post infection",
                     limits = c(0, 75))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15),
        legend.position = "none")   

c <- ggplot(data = groupe_R, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "numero_patient ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("R")+
  scale_x_continuous(name = "jours post infection",
                     limits = c(0, 75))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15),
        legend.position = "none") 

d <- ggplot(data = groupe_RPm, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "numero_patient ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("RP-")+
  scale_x_continuous(name = "jours post infection",
                     limits = c(0, 75))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15),
        legend.position = "none")   

e <- ggplot(data = groupe_RP, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "numero_patient ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("RP")+
  scale_x_continuous(name = "jours post infection",
                     limits = c(0, 75))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15),
        legend.position = "none")  


plot_grid(a, b, c, d, e, labels=c("A", "B", "C", "D", "E"), ncol = 5, nrow = 1)  # Plot avec 3 figures


# Bornage 0 à 20 jours-----
f <- ggplot(data = groupe_NR, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "numero_patient ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("NR")+
  scale_x_continuous(name = "jorus post infection",
                     limits = c(0, 20))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15),
        legend.position = "none")   

g <- ggplot(data = groupe_NRm, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "numero_patient ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("NR-")+
  scale_x_continuous(name = "jorus post infection",
                     limits = c(0, 20))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15),
        legend.position = "none") 

h <- ggplot(data = groupe_R, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "numero_patient ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("R")+
  scale_x_continuous(name = "jorus post infection",
                     limits = c(0, 20))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15),
        legend.position = "none")   

i <- ggplot(data = groupe_RPm, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "numero_patient ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("RP-")+
  scale_x_continuous(name = "jorus post infection",
                     limits = c(0, 20))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15),
        legend.position = "none") 

j <- ggplot(data = groupe_RP, 
            aes_string(x="day", 
                       y= "Vidas", 
                       color = "numero_patient ", 
                       group = "numero_patient"))+ 
  geom_line() + 
  geom_point() + 
  ggtitle("RP")+
  scale_x_continuous(name = "jorus post infection",
                     limits = c(0, 20))+
  scale_y_continuous(name = "Titre en Ac",
                     limits = c(-0.01,43))+
  theme(text = element_text(size = 15),
        legend.position = "none") 


plot_grid(a, b, c, d, e, f, g, h, i, j, 
          labels=c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"), 
          ncol = 5, 
          nrow = 2)  # Plot avec 3 figures









































# Bornage pour nouveau time point ----
data2 <- rbind(groupe_R, groupe_RP)

T_F <- data2$real_time_point %in% "VT1"
data_VT1 <- data2[T_F==T,]

T_F <- data2$real_time_point %in% "VT2"
data_VT2 <- data2[T_F==T,]

T_F <- data2$real_time_point %in% "VT3"
data_VT3 <- data2[T_F==T,]

T_F <- data2$real_time_point %in% "VT4"
data_VT4 <- data2[T_F==T,]

T_F <- data2$real_time_point %in% "VT5"
data_VT5 <- data2[T_F==T,]

T_F <- data2$real_time_point %in% "VT6"
data_VT6 <- data2[T_F==T,]

T_F <- data2$real_time_point %in% "VT7"
data_VT7 <- data2[T_F==T,]




data_VT1<- data_VT1%>% 
  reorder_levels(Groupe, order = c( "R", "RP"))

# Test stat
res.wilcox <- data_VT1 %>%
  wilcox_test(Vidas ~ Groupe) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Groupe")
a <- ggboxplot(data_VT1, 
                 x = "Groupe", 
                 y = "Vidas", 
                 ylab = "Vidas", 
                 xlab = "Groups", 
                 add = "point") + 
  stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
  labs(title = "VT1",
       subtitle = get_test_label(res.wilcox, detailed = T))

data_VT2<- data_VT2%>% 
  reorder_levels(Groupe, order = c( "R", "RP"))

# Test stat
res.wilcox <- data_VT2 %>%
  wilcox_test(Vidas ~ Groupe) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Groupe")
b <- ggboxplot(data_VT2, 
               x = "Groupe", 
               y = "Vidas", 
               ylab = "Vidas", 
               xlab = "Groups", 
               add = "point") + 
  stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
  labs(title = "VT2",
       subtitle = get_test_label(res.wilcox, detailed = T))

data_VT3<- data_VT3 %>% 
  reorder_levels(Groupe, order = c( "R", "RP"))

# Test stat
res.wilcox <- data_VT3 %>%
  wilcox_test(Vidas ~ Groupe) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Groupe")
c <- ggboxplot(data_VT3, 
               x = "Groupe", 
               y = "Vidas", 
               ylab = "Vidas", 
               xlab = "Groups", 
               add = "point") + 
  stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
  labs(title = "VT3"
    ,subtitle = get_test_label(res.wilcox, detailed = T))

data_VT4 <- data_VT4 %>% 
  reorder_levels(Groupe, order = c( "R", "RP"))

# Test stat
res.wilcox <- data_VT4 %>%
  wilcox_test(Vidas ~ Groupe) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Groupe")
d <- ggboxplot(data_VT4, 
               x = "Groupe", 
               y = "Vidas", 
               ylab = "Vidas", 
               xlab = "Groups", 
               add = "point") + 
  stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
  labs(title = "VT4",
       subtitle = get_test_label(res.wilcox, detailed = T))

data_VT5 <- data_VT5 %>% 
  reorder_levels(Groupe, order = c( "R", "RP"))

# Test stat
res.wilcox <- data_VT5 %>%
  wilcox_test(Vidas ~ Groupe) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Groupe")
e <- ggboxplot(data_VT5, 
               x = "Groupe", 
               y = "Vidas", 
               ylab = "Vidas", 
               xlab = "Groups", 
               add = "point") + 
  stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
  labs(title = "VT5",
       subtitle = get_test_label(res.wilcox, detailed = T))


data_VT6 <- data_VT6 %>% 
  reorder_levels(Groupe, order = c( "R", "RP"))

# Test stat
res.wilcox <- data_VT6 %>%
  wilcox_test(Vidas ~ Groupe) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Groupe")
f <- ggboxplot(data_VT6, 
               x = "Groupe", 
               y = "Vidas", 
               ylab = "Vidas", 
               xlab = "Groups", 
               add = "point") + 
  stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
  labs(title = "VT6",
       subtitle = get_test_label(res.wilcox, detailed = T))

data_VT7<- data_VT7%>% 
  reorder_levels(Groupe, order = c( "R", "RP"))

# Test stat
res.wilcox <- data_VT7 %>%
  wilcox_test(Vidas ~ Groupe) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Groupe")
g <- ggboxplot(data_VT7, 
               x = "Groupe", 
               y = "Vidas", 
               ylab = "Vidas", 
               xlab = "Groups", 
               add = "point") + 
  stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
  labs(title = "VT7",
       subtitle = get_test_label(res.wilcox, detailed = T))


plot_grid(a, b, c, d, e, f, g, labels=c("A", "B", "C", "D", "E", "F", "G"), ncol = 4, nrow = 2)  # Plot avec 3 figures






