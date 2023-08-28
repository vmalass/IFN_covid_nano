Etude du titre en anti corps 

library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(rstatix)
library(cowplot)

# import datat----
rm(list = ls())
data <- read.xlsx("data/titre_virale_IFI27_Ac.xlsx", sheet = 7)
data$numero_patient<- as.factor(data$numero_patient)

# Visu all time point ----
ggplot(data = data, 
       aes(x = factor(real_time_point, levels = c("VT1", "VT2", "VT3", "VT4", "VT5", "VT6", "VT7", "VT8", "S1", "S2")), 
           y = Vidas, 
           color = Groupe))+ 
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  labs(title = "Wilcoxon test for title antibodies", 
       x = "Time point", 
       y = "quantity of antibodies (XXXX)")+
  theme_classic() + 
  theme(legend.position = c(0.95,0.8))

dplyr::count(data, real_time_point, Groupe)

# Visu VT1 at VT4----
T_F <- data$real_time_point %in% c("VT5", "VT6", "VT7", "VT8", "S1", "S2")
data2 <- data[T_F==F,]
dplyr::count(data2, real_time_point, Groupe)

ggplot(data = data2, 
       aes(x = factor(real_time_point, levels = c("VT1", "VT2", "VT3", "VT4")), 
           y = Vidas, 
           color = Groupe))+ 
  geom_boxplot()+
  geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
  theme_classic() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(y_position = c(2, 15, 32, 32), 
              xmin = c(0.8, 1.8, 2.8, 3.8), 
              xmax = c(1.2, 2.2, 3.2, 4.2),
              annotation = c("0.74", "**", "0.23", "0.3")) + 
  labs(title = "Wilcoxon test for title antibodies", 
       subtitle = "VT1 11R/3RP, VT2 12R/5RP, VT3 12R/5RP, VT4 15R/6RP",
       x = "", 
       y = "RBD-specific IgG (vidas BAU/mL)")+
  theme_classic()  


# Data par time point----
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



# Test stat par time point----

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
               add = "point",
               color = "Groupe",
               palette = c("#F15854", "#882255"),
               legend = c(0.9,0.95)) +
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
               add = "point",
               color = "Groupe",
               palette = c("#F15854", "#882255"),
               legend = c(0.9,0.95)) + 
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
               add = "point",
               color = "Groupe",
               palette = c("#F15854", "#882255"),
               legend = c(0.9,0.95)) + 
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
               add = "point",
               color = "Groupe",
               palette = c("#F15854", "#882255"),
               legend = c(0.9,0.95)) + 
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
               add = "point",
               color = "Groupe",
               palette = c("#F15854", "#882255"),
               legend = c(0.9,0.95)) + 
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
               add = "point",
               color = "Groupe",
               palette = c("#F15854", "#882255"),
               legend = c(0.9,0.95)) + 
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
               add = "point",
               color = "Groupe",
               palette = c("#F15854", "#882255"),
               legend = c(0.9,0.95)) + 
  stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
  labs(title = "VT7",
       subtitle = get_test_label(res.wilcox, detailed = T))


plot_grid(a, b, c, d, e, f, g, labels=c("A", "B", "C", "D", "E", "F", "G"), ncol = 4, nrow = 2)  # Plot avec 3 figures





T_F <- data$real_time_point %in% "S1"
data_S1 <- data[T_F==T,]

T_F <- data$real_time_point %in% "S2"
data_S2 <- data[T_F==T,]


# Test stat
res.wilcox <- data_S1 %>%
  wilcox_test(Vidas ~ Groupe) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Groupe")
a <- ggboxplot(data_S1, 
               x = "Groupe", 
               y = "Vidas", 
               ylab = "Vidas", 
               xlab = "Groups", 
               add = "point",
               color = "Groupe",
               palette = c("#F15854", "#882255"),
               legend = c(0.9,0.95)) + 
  stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
  labs(title = "S1",
       subtitle = get_test_label(res.wilcox, detailed = T))


# Test stat
res.wilcox <- data_S2 %>%
  wilcox_test(Vidas ~ Groupe) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Groupe")
b <- ggboxplot(data_S2, 
               x = "Groupe", 
               y = "Vidas", 
               ylab = "Vidas", 
               xlab = "Groups", 
               add = "point",
               color = "Groupe",
               palette = c("#F15854", "#882255"),
               legend = c(0.9,0.95)) + 
  stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
  labs(title = "S2",
       subtitle = get_test_label(res.wilcox, detailed = T))

plot_grid(a, b, labels=c("A", "B"), ncol = 2, nrow = 1)  # Plot avec 3 figures


#-- Test statistique, non apparié et non paramétrique (transfrmation log non utile dans ce cas).
# Est-ce que la concentration cytokinique (IFN_l1_IL_29, IFN_gamma & IFN_beta confondues) est 
# différente entre les 2 groupes, R et RP ?
with(data2 , wilcox.test(Vidas ~ Groupe))

ggplot(data2 , aes(x =  Groupe , y = Vidas))+
  geom_boxplot( outlier.shape = NA) + #  outlier.shape = NA
  # geom_point() +
  geom_jitter() + # aes(col = cytokine), width = 0.2, size = 7, alpha = 0.5
  geom_signif(comparisons = list(c("R", "RP")),
              map_signif_level = TRUE) +
  labs(title = "Wilcoxon test for IFN quantity", 
       y = "IFN (pg/mL)") +
  theme_classic() 
