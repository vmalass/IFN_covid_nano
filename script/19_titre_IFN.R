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

# IFN vs temps------------------------------------------------------------------

ggplot(data = V1, 
       aes(x = delai.depuis.1ers.symptomes, 
           y = IFN, 
           color = Reponse, 
           group = numero_patient_victor))+ 
  geom_point(size = 4) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  ggtitle("IFN à V1 vs tps")+
  theme_bw() 

ggplot(data = V4, 
       aes(x = delai.depuis.1ers.symptomes, 
           y = IFN, 
           color = Reponse, 
           group = numero_patient_victor))+ 
  geom_point(size = 4) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  ggtitle("IFN à V4 vs tps")+
  theme_bw() 

ggplot(data = data, 
       aes(x = delai.depuis.1ers.symptomes, 
           y = IFN, 
           color = Reponse, 
           group = numero_patient_victor))+ 
  geom_point(size = 4) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  ggtitle("Cinétique IFN vs tps")+
  theme_bw() 

ggplot(data = data_all, 
       aes(x = delai.depuis.1ers.symptomes, 
           y = IFN, 
           color = Reponse, 
           label = numero_patient_victor))+ 
  geom_point(size = 3) +
  geom_text_repel(show.legend = F,
                  max.overlaps  = Inf) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A", "NC"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB", "black")) +
  ggtitle("Cinétique IFN vs tps")+
  theme_bw()

# Visualisation cinétique IFN---------------------------------------------------
data <- data %>% 
  reorder_levels(new_time_point, order = c( "VT1", "VT2", "VT3", "VT4"))

ggplot(data = data, 
       aes(x = delai.depuis.1ers.symptomes, 
           y = IFN, 
           color = Reponse, 
           group = numero_patient_victor))+ 
  geom_line() +
  geom_point() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  ggtitle("Cinétique IFN ")+
  theme_bw()

ggplot(data = data, 
       aes(x = new_time_point, 
           y = IFN, 
           color = Reponse))+ 
  geom_boxplot() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  labs(title = "IFN", 
       subtitle = "Wilcoxon test",
       x = "Time point", 
       y = "quantity of IFN")+
  theme_bw() +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 14))
# Test stat---------------------------------------------------------------------
summary(data)
print(data %>% 
        group_by(new_time_point) %>%
        get_summary_stats(IFN, type = "common"))    #### problème avec titre col ####

# Test kruskal
res.kruskal <- data %>% kruskal_test(formula(paste0( 'IFN'," ~ new_time_point")))
print(res.kruskal)

# Taille de l’effet
print(data %>% kruskal_effsize(formula(paste0('IFN'," ~ new_time_point"))))

# Comparaisons par paires test de Dunn
pwc <- data %>% 
  dunn_test(formula(paste0('IFN'," ~ new_time_point")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "new_time_point")
print(ggboxplot(data, x = "new_time_point", y = 'IFN', color = "Reponse", palette = c("cornflowerblue","brown3")) +
        # geom_point( color = "grey70") +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(subtitle = get_test_label(res.kruskal, detailed = TRUE),
             caption = get_pwc_label(pwc))) +
  labs(title = "Titre en IFN en fonction des groupes",
       y = "Titre en IFN") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=20)) +
  theme(legend.position="right")




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
                title = "Day of sampling after symptom in V1",
                x = "Reponse", 
                y = paste0("IFN"), 
                ylab = "IFN",
                xlab = "Groups", 
                add = "point") + 
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
                title = "Day of sampling after symptom in V2",
                x = "Reponse", 
                y = paste0("IFN"), 
                ylab = "IFN",
                xlab = "Groups", 
                add = "point") + 
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
                title = "Day of sampling after symptom in V3",
                x = "Reponse", 
                y = paste0("IFN"), 
                ylab = "IFN",
                xlab = "Groups", 
                add = "point") + 
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
                title = "Day of sampling after symptom in V4",
                x = "Reponse", 
                y = paste0("IFN"), 
                ylab = "IFN",
                xlab = "Groups", 
                add = "point") + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))






