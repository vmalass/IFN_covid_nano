library(readxl)
library(ggplot2)
library(ggrepel)
library(rstatix)
library(ggpubr)
library(dplyr)

# Import data-------------------------------------------------------------------
rm(list = ls())
charge_vir <- read_xlsx("data/titre_virale_IFI27_ac.xlsx", sheet = 3)

# all data----------------------------------------------------------------------
charge_vir <- na.omit(charge_vir)

ggplot(charge_vir, aes(x = jours_prelevement, 
                       y = `charge_virale_log10`, 
                       color = Reponse, 
                       label = numero_patient_victor))+
  geom_point(size = 3) +
  geom_text_repel(show.legend = F,
                  max.overlaps  = Inf) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A", "NC"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB", "black")) +
  ggtitle("Charge virale de tous les patients")+
  theme_bw()

# Uniquement R & RP-------------------------------------------------------------
charge_vir_filter <- charge_vir[charge_vir$Reponse %in% c("R", "RP"),]

ggplot(charge_vir_filter, aes(x = jours_prelevement, 
                       y = `charge_virale_log10`, 
                       color = Reponse, 
                       label = numero_patient_victor))+
  geom_point(size = 3) +
  geom_text_repel(show.legend = F,
                  max.overlaps  = Inf) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A", "NC"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB", "black")) +
  ggtitle("Charge virale de tous les patients R et RP")

# Identification des patients qui on VT1 & VT2----------------------------------
VT1 <- dplyr::filter(charge_vir, new_time_point == "VT1")
VT2 <- dplyr::filter(charge_vir, new_time_point == "VT2")

a <- intersect(VT1$numero_patient_victor, VT2$numero_patient_victor)

charge_vir_filter_2 <- charge_vir[which(charge_vir$numero_patient_victor %in% a),]
b <- unique(charge_vir_filter_2$numero_patient_victor)

ggplot(charge_vir_filter_2, aes(x = jours_prelevement, 
                              y = `charge_virale_log10`, 
                              color = Reponse, 
                              label = numero_patient_victor))+
  geom_point(size = 3) +
  geom_text_repel(show.legend = F,
                  max.overlaps  = Inf) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A", "NC"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB", "black")) +
  ggtitle("Charge virale de tous les patients ayant un prélèvement à VT1 et VT2")

# Identification des patients R et RP qui on VT1 & VT2--------------------------
charge_vir_filter_3 <- charge_vir_filter_2[charge_vir_filter_2$Reponse %in% c("R", "RP"),]

ggplot(charge_vir_filter_3, aes(x = jours_prelevement, 
                                y = `charge_virale_log10`, 
                                color = Reponse, 
                                label = numero_patient_victor))+
  geom_point(size = 3) +
  geom_text_repel(show.legend = F,
                  max.overlaps  = Inf) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A", "NC"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB", "black"))  +
  ggtitle("Charge virale de tous les patients R et RP ayant un prélèvement à VT1 et VT2")

# Identification des patients avec 2 times points-------------------------------
VT3 <- dplyr::filter(charge_vir, new_time_point == "VT3")
VT4 <- dplyr::filter(charge_vir, new_time_point == "VT4")

a <- intersect(VT1$numero_patient_victor, VT2$numero_patient_victor)
b <- intersect(VT1$numero_patient_victor, VT3$numero_patient_victor)
c <- intersect(VT1$numero_patient_victor, VT4$numero_patient_victor)

d <- union(a, b)
d <- union(d, c)

charge_vir_filter_4 <- charge_vir[which(charge_vir$numero_patient_victor %in% d),]

ggplot(charge_vir_filter_4, aes(x = jours_prelevement, 
                                y = `charge_virale_log10`, 
                                color = Reponse, 
                                label = numero_patient_victor))+
  geom_point(size = 3) +
  geom_text_repel(show.legend = F,
                  max.overlaps  = Inf) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A", "NC"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB", "black")) +
  ggtitle("Charge virale de tous les patients ayant deux prélèvements minimum")

# Identification des patients R et RP qui on deux VT mini-----------------------
charge_vir_filter_5 <- charge_vir_filter_4[charge_vir_filter_4$Reponse %in% c("R", "RP"),]

ggplot(charge_vir_filter_5, aes(x = jours_prelevement, 
                                y = `charge_virale_log10`, 
                                color = Reponse, 
                                label = numero_patient_victor))+
  geom_point(size = 3) +
  geom_text_repel(show.legend = F,
                  max.overlaps  = Inf) +
  scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A", "NC"),
                     values = c("darkorange", "#DDCC77","cornflowerblue",
                                "#882255" ,"brown3", "chartreuse4", "#BBBBBB", "black")) +
  ggtitle("Charge virale de tous les patients R et RP ayant deux prélèvements minimum")

## Cinétique--------------------------------------------------------------------
charge_vir_filter_5 <- charge_vir_filter_5 %>% 
  reorder_levels(new_time_point, order = c( "VT1", "VT2", "VT3", "VT4"))

ggplot(data = charge_vir_filter_5, 
       aes(x = jours_prelevement, 
           y = `charge_virale_log10`, 
           color = Reponse, 
           group = numero_patient_victor))+ 
  geom_line() +
  geom_point() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  ggtitle("Cinétique charge virale de tous les patients R et RP ayant deux prélèvements minimum ")+
  theme_bw()


## Test stat---------------------------------------------------------------------
summary(charge_vir_filter_5)

print(charge_vir_filter_5 %>% 
        group_by(new_time_point) %>%
        get_summary_stats(charge_virale_log10, type = "common"))    #### problème avec titre col ####

# Test kruskal
res.kruskal <- charge_vir_filter_5 %>% kruskal_test(formula(paste0( 'charge_virale_log10'," ~ new_time_point")))
print(res.kruskal)

# Taille de l’effet
print(charge_vir_filter_5 %>% kruskal_effsize(formula(paste0('charge_virale_log10'," ~ new_time_point"))))

# Comparaisons par paires test de Dunn
pwc <- charge_vir_filter_5 %>% 
  dunn_test(formula(paste0('charge_virale_log10'," ~ new_time_point")), p.adjust.method = "bonferroni") 
print(pwc)

ggplot(data = charge_vir_filter_5, 
       aes(x = new_time_point, 
           y = `charge_virale_log10`, 
           color = Reponse))+ 
  geom_boxplot() + 
  geom_dotplot(binaxis='y', 
               stackdir='center',
               position=position_dodge(0.75), 
               binwidth = 0.15,
               fill = "white")+
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  labs(title = "Viral load", 
       subtitle = "Kruuskal test and Dunn test",
       x = "Time point", 
       y = "Viral load")+
  theme_bw() +
  theme(axis.text = element_text(size = 20), 
        axis.title = element_text(size = 20), 
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 14))


## Test stat R vs RP à VT1/2/3/4------------------------------------------------
### Viral load V1 ###
V1 <- charge_vir_filter_5[charge_vir_filter_5$new_time_point == "VT1",]

V1 <- V1%>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V1 %>%
  wilcox_test(formula(paste0("charge_virale_log10"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V1, 
                title = "Day of sampling after symptom in V1",
                x = "Reponse", 
                y = paste0("charge_virale_log10"), 
                ylab = "charge_virale_log10",
                xlab = "Groups", 
                add = "point") + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))

### Viral load V2 ###
V2 <- charge_vir_filter_5[charge_vir_filter_5$new_time_point == "VT2",]
V2 <- V2 %>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V2 %>%
  wilcox_test(formula(paste0("charge_virale_log10"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V2, 
                title = "Day of sampling after symptom in V2",
                x = "Reponse", 
                y = paste0("charge_virale_log10"), 
                ylab = "charge_virale_log10",
                xlab = "Groups", 
                add = "point") + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))


### Viral load V3 ###
V3 <- charge_vir_filter_5[charge_vir_filter_5$new_time_point == "VT3",]
V3 <- V3 %>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V3 %>%
  wilcox_test(formula(paste0("charge_virale_log10"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V3, 
                title = "Day of sampling after symptom in V3",
                x = "Reponse", 
                y = paste0("charge_virale_log10"), 
                ylab = "charge_virale_log10",
                xlab = "Groups", 
                add = "point") + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))

### Viral load V4 ###
V4 <- charge_vir_filter_5[charge_vir_filter_5$new_time_point == "VT4",]
V4 <- V4 %>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",

# Test stat
res.wilcox <- V4 %>%
  wilcox_test(formula(paste0("charge_virale_log10"," ~ Reponse"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
print(ggboxplot(V4, 
                title = "Day of sampling after symptom in V4",
                x = "Reponse", 
                y = paste0("charge_virale_log10"), 
                ylab = "charge_virale_log10",
                xlab = "Groups", 
                add = "point") + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))


