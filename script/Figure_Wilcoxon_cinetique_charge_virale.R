Cinétique de la charge virale 

library(readxl)
library(ggplot2)
library(ggrepel)
library(rstatix)
library(ggpubr)
library(dplyr)

# Import data-------------------------------------------------------------------
rm(list = ls())
charge_vir <- read_xlsx("data/titre_virale_IFI27_ac.xlsx", sheet = 3)

charge_vir_filter <- charge_vir[charge_vir$Reponse %in% c("R", "RP"),]

charge_vir_filter <- na.omit(charge_vir_filter)

count(charge_vir_filter, new_time_point, Reponse)


## Visualisation cinétique------------------------------------------------------
ggplot(data = charge_vir_filter, 
       aes(x = new_time_point, 
           y = charge_virale_log10, 
           color = Reponse))+ 
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  labs(title = "Wilcoxon test for viral load", 
       x = "Time point", 
       y = "Viral load (log10)")+
  theme_classic() + 
  theme(legend.position = c(0.9,0.8))

# Test stat time point R/RP-----------------------------------------------------

### charge_virale_log10 V1 ###
V1 <- charge_vir_filter[charge_vir_filter$new_time_point == "VT1",]

V1 <- V1%>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V1 %>%
  wilcox_test(formula(paste0("charge_virale_log10"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V1, 
                title = "Wilcoxon test for viral load in V1",
                x = "Reponse", 
                y = paste0("charge_virale_log10"), 
                ylab = "Viral load (log10)",
                xlab = "Groups", 
                add = "point",
                color = "Reponse",
                palette = c("gray50", "#CB2027"),
                legend = c(0.9,0.95)) + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))

### charge_virale_log10 V2 ###
V2 <- charge_vir_filter[charge_vir_filter$new_time_point == "VT2",]
V2 <- V2 %>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V2 %>%
  wilcox_test(formula(paste0("charge_virale_log10"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V2, 
                title = "Wilcoxon test for viral load in V2",
                x = "Reponse", 
                y = paste0("charge_virale_log10"), 
                ylab = "Viral load (log10)",
                xlab = "Groups", 
                add = "point",
                color = "Reponse",
                palette = c("gray50", "#CB2027"),
                legend = c(0.9,0.95)) + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))


### charge_virale_log10 V3 ###
V3 <- charge_vir_filter[charge_vir_filter$new_time_point == "VT3",]
V3 <- V3 %>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V3 %>%
  wilcox_test(formula(paste0("charge_virale_log10"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V3, 
                title = "Wilcoxon test for viral load in V3",
                x = "Reponse", 
                y = paste0("charge_virale_log10"), 
                ylab = "Viral load (log10)",
                xlab = "Groups", 
                add = "point",
                color = "Reponse",
                palette = c("gray50", "#CB2027"),
                legend = c(0.9,0.95)) + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))

### charge_virale_log10 V4 ###
V4 <- charge_vir_filter[charge_vir_filter$new_time_point == "VT4",]
V4 <- V4 %>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V4 %>%
  wilcox_test(formula(paste0("charge_virale_log10"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V4, 
                title = "Wilcoxon test for viral load in V4",
                x = "Reponse", 
                y = paste0("charge_virale_log10"), 
                ylab = "Viral load (log10)",
                xlab = "Groups", 
                add = "point",
                color = "Reponse",
                palette = c("gray50", "#CB2027"),
                legend = c(0.9,0.95)) + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))
