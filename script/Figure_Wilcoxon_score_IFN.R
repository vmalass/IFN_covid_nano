Score IFN type I/III issu du papier de Sophie 

library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(rstatix)
library(cowplot)

# import datat------------------------------------------------------------------
rm(list = ls())

data <- read.xlsx("data/IFN_titre.xlsx", sheet = 1)
data_all <- read.xlsx("data/IFN_titre.xlsx", sheet = 2)
V1 <- read.xlsx("data/IFN_titre.xlsx", sheet = 3)
V4 <- read.xlsx("data/IFN_titre.xlsx", sheet = 4)

# data$new_time_point <- as.factor(data$new_time_point)

T_F <- data_all$delai.depuis.1ers.symptomes %in% "X"
data_all <- data_all[T_F == F,]
data_all$delai.depuis.1ers.symptomes <- as.numeric(data_all$delai.depuis.1ers.symptomes)
data_all$IFN <- as.numeric(data_all$IFN)

# Visualisation cinétique IFN---------------------------------------------------
ggplot(data = data, 
       aes(x = new_time_point, 
           y = IFN, 
           color = Reponse))+ 
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  labs(title = "Wilcoxon test for title IFN", 
       x = "Time point", 
       y = "Score type I/III IFN Nanostring")+
  theme_classic() + 
  theme(legend.position = c(0.9,0.8))

### IFN V1 ###
V1 <- data[data$new_time_point == "VT1",]

V1 <- V1%>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V1 %>%
  wilcox_test(formula(paste0("IFN"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V1, 
                title = "Wilcoxon test for title IFN in V1",
                x = "Reponse", 
                y = paste0("IFN"), 
                ylab = "quantity of IFN (fg/mL)",
                xlab = "Groups", 
                add = "point",
                color = "Reponse",
                palette = c("gray50", "#CB2027"),
                legend = c(0.9,0.95)) + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))

### IFN V2 ###
V2 <- data[data$new_time_point == "VT2",]
V2 <- V2 %>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V2 %>%
  wilcox_test(formula(paste0("IFN"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V2, 
                title = "Wilcoxon test for title IFN in V2",
                x = "Reponse", 
                y = paste0("IFN"), 
                ylab = "quantity of IFN (fg/mL)",
                xlab = "Groups", 
                add = "point",
                color = "Reponse",
                palette = c("gray50", "#CB2027"),
                legend = c(0.9,0.95)) + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))


### IFN V3 ###
V3 <- data[data$new_time_point == "VT3",]
V3 <- V3 %>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V3 %>%
  wilcox_test(formula(paste0("IFN"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V3, 
                title = "Wilcoxon test for title IFN in V3",
                x = "Reponse", 
                y = paste0("IFN"), 
                ylab = "quantity of IFN (fg/mL)",
                xlab = "Groups", 
                add = "point",
                color = "Reponse",
                palette = c("gray50", "#CB2027"),
                legend = c(0.9,0.95)) + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))

### IFN V4 ###
V4 <- data[data$new_time_point == "VT4",]
V4 <- V4 %>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V4 %>%
  wilcox_test(formula(paste0("IFN"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V4, 
                title = "Wilcoxon test for title IFN in V4",
                x = "Reponse", 
                y = paste0("IFN"), 
                ylab = "quantity of IFN (fg/mL)",
                xlab = "Groups", 
                add = "point",
                color = "Reponse",
                palette = c("gray50", "#CB2027"),
                legend = c(0.9,0.95)) + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))

