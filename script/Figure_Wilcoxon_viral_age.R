Test de wilcoxon sur la charge virale à V1

library(readxl)
library(ggplot2)
library(ggrepel)
library(rstatix)
library(ggpubr)

# Import data-------------------------------------------------------------------
rm(list = ls())
charge_vir <- read_xlsx("data/unique_patient.xlsx")
charge_vir <- charge_vir[charge_vir$Reponse != "NC",]
charge_vir <- charge_vir[charge_vir$Reponse != "NR",]
charge_vir <- charge_vir[charge_vir$Reponse != "NR-",]
charge_vir <- charge_vir[charge_vir$Reponse != "RP-",]
charge_vir <- charge_vir[charge_vir$Reponse != "A",]

charge_vir<- charge_vir%>% 
  reorder_levels(Reponse, order = c( "R", "RP")) 

### Charge virale ###

# Test stat
res.wilcox <- charge_vir %>%
  wilcox_test(formula(paste0("charge_virale_log10"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(charge_vir, 
                title = "Viral load according to groups",
                x = "Reponse", 
                y = paste0("charge_virale_log10"), 
                ylab = "log10 /1 000 000 cells",
                xlab = "Groups",
                add = "point",
                color = "Reponse",
                palette = c("gray60","#CB2027"),
                legend = c(0.9,0.95)) +  
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))


### Age ###

# Test stat
res.wilcox <- charge_vir %>%
  wilcox_test(formula(paste0("Age"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(charge_vir, 
                title = "Age according to groups",
                x = "Reponse", 
                y = paste0("Age"), 
                ylab = "Age",
                xlab = "Groups", 
                add = "point",
                color = "Reponse",
                palette = c("gray60","#CB2027"),
                legend = c(0.9,0.95)) +  
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))
