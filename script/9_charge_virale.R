library(readxl)
library(ggplot2)
library(ggrepel)
library(rstatix)

# Import data-------------------------------------------------------------------
rm(list = ls())
charge_vir <- read_xlsx("data/unique_patient.xlsx")
charge_vir <- charge_vir[charge_vir$Reponse != "NC",]

# Visualisation & cor test------------------------------------------------------
test_cor <- cor.test(charge_vir$jours_prelevement, charge_vir$`charge_virale_log10`, method= "pearson")

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
  annotate(geom="text", 
           x=5, 
           y=1.5, 
           label= paste0("cor : ",round(test_cor$estimate, 2), " pvalue : ", signif(test_cor$p.value, 3)),
           color="red",
           size = 8) +
  ggtitle("Charge virale en fonction du jours de prélèvement") +
  theme_bw() +
  theme(plot.title = element_text(size=20),
        axis.title = element_text(size = 13), 
        axis.text = element_text(size = 13),
        legend.title = element_text(size=16),
        legend.text = element_text(size=15)) 

# Test stat---------------------------------------------------------------------
charge_vir %>% sample_n_by(Reponse, size = 1)
data <- charge_vir %>%
  reorder_levels(Reponse, order = c("NR", "NR-", "R", "RP-", "RP"))
data <- na.omit(data)
summary(data)
print(data %>% 
        group_by(Reponse) %>%
        get_summary_stats(charge_virale_log10, type = "common"))    #### problème avec titre col ####

# Test kruskal
res.kruskal <- data %>% kruskal_test(formula(paste0( 'charge_virale_log10'," ~ Reponse")))
print(res.kruskal)

# Taille de l’effet
print(data %>% kruskal_effsize(formula(paste0('charge_virale_log10'," ~ Reponse"))))

# Comparaisons par paires test de Dunn
pwc <- data %>% 
  dunn_test(formula(paste0('charge_virale_log10'," ~ Reponse")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "Reponse")
print(ggboxplot(data, x = "Reponse", y = 'charge_virale_log10') +
        geom_point( color = "grey70") +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(subtitle = get_test_label(res.kruskal, detailed = TRUE),
             caption = get_pwc_label(pwc))) +
  labs(title = "Charge virale en fonction des groupes",
       y = "log10copies /1 000 000 cellules") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=20))








