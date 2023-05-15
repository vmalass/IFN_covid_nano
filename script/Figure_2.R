library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(rstatix)
library(cowplot)

# import datat------------------------------------------------------------------
rm(list = ls())
data <- read.xlsx("data/titre_virale_IFI27_Ac.xlsx", sheet = 6)
unique_pat <- read.xlsx("data/unique_patient.xlsx")
V1 <- read.xlsx("data/Frequence_cell_V1_V4.xlsx", sheet = 1)
V4 <- read.xlsx("data/Frequence_cell_V1_V4.xlsx", sheet = 2)


data$numero_patient<- as.factor(data$numero_patient)
data$real_time_point<- as.factor(data$real_time_point)

# Elimination VT8 S1 S2---------------------------------------------------------
data_filtre <- data[data$real_time_point != "VT8",]
data_filtre <- data_filtre[data_filtre$real_time_point != "S1",]
data_filtre <- data_filtre[data_filtre$real_time_point != "S2",]

#
ggplot(data = V1, 
       aes(x = jours_prelevement, 
           y = MONO, 
           color = response, 
           group = numero_patient_victor))+ 
  geom_point(size = 4) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  ggtitle("Monocyte")+
  theme_bw() +

ggplot(data = V4, 
       aes(x = jours_prelevement, 
           y = `T_cell_IFNg+`, 
           color = Responder, 
           group = numero_patient_victor))+ 
  geom_point(size = 4) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  ggtitle("T cell IFNg")+
  theme_bw() +

ggplot(data = V4, 
       aes(x = jours_prelevement, 
           y = `CD8+_IFNg+`, 
           color = Responder, 
           group = numero_patient_victor))+ 
  geom_point(size = 4) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  ggtitle("CD8+ IFNg+")+
  theme_bw()

# Visualisation cinétique boxplot-----------------------------------------------
data_filtre <- data_filtre %>% 
  reorder_levels(real_time_point, order = c( "VT1", "VT2", "VT3", "VT4", "VT5", "VT6", "VT7"))

ggplot(data = data_filtre, 
       aes(x = real_time_point, 
           y = Vidas, 
           color = Groupe, 
           group = numero_patient))+ 
  geom_line() +
  geom_point() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  ggtitle("Titre en anti-corps")+
  theme_bw()


ggplot(data = data_filtre, 
       aes(x = real_time_point, 
           y = Vidas, 
           color = Groupe))+ 
  geom_boxplot() +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("cornflowerblue","brown3")) +
  ggtitle("Titre en anti-corps")+
  theme_bw()


# Test stat---------------------------------------------------------------------
data_filtre <- na.omit(data_filtre)
summary(data_filtre)
print(data_filtre %>% 
        group_by(real_time_point) %>%
        get_summary_stats(Vidas, type = "common"))    #### problème avec titre col ####

# Test kruskal
res.kruskal <- data_filtre %>% kruskal_test(formula(paste0( 'Vidas'," ~ real_time_point")))
print(res.kruskal)

# Taille de l’effet
print(data_filtre %>% kruskal_effsize(formula(paste0('Vidas'," ~ real_time_point"))))

# Comparaisons par paires test de Dunn
pwc <- data_filtre %>% 
  dunn_test(formula(paste0('Vidas'," ~ real_time_point")), p.adjust.method = "bonferroni") 
print(pwc)

# Visualisation : Boxplots avec p-values
pwc <- pwc %>% add_xy_position(x = "real_time_point")
print(ggboxplot(data_filtre, x = "real_time_point", y = 'Vidas', color = "Groupe", palette = c("cornflowerblue","brown3")) +
        # geom_point( color = "grey70") +
        stat_pvalue_manual(pwc, hide.ns = TRUE) +
        labs(subtitle = get_test_label(res.kruskal, detailed = TRUE),
             caption = get_pwc_label(pwc))) +
  labs(title = "Titre en Ac en fonction des groupes",
       y = "Titre en Anticorps") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size=20)) +
  theme(legend.position="right")








data_filtre <- data[data$time_point == "day_V1",]

data_filtre<- data_filtre%>% 
  reorder_levels(Groupe, order = c( "R", "RP"))  ## "T",

### Viral load ###

# Test stat
res.wilcox <- data_filtre %>%
  wilcox_test(formula(paste0("day"," ~ Groupe"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_filtre, 
                title = "Day of sampling after symptom in V1",
                x = "Groupe", 
                y = paste0("day"), 
                ylab = "day",
                xlab = "Groups", 
                add = "point") + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))



data_filtre <- data[data$time_point == "day_V4",]

data_filtre<- data_filtre%>% 
  reorder_levels(Groupe, order = c( "R", "RP"))  ## "T",

### Viral load ###

# Test stat
res.wilcox <- data_filtre %>%
  wilcox_test(formula(paste0("day"," ~ Groupe"))) %>%
  add_significance()

# Visualization
res.wilcox <- res.wilcox %>% add_xy_position(x = "Groupe")
print(ggboxplot(data_filtre, 
                title = "Day of sampling after symptom in V4",
                x = "Groupe", 
                y = paste0("day"), 
                ylab = "day",
                xlab = "Groups", 
                add = "point") + 
        stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
        labs(subtitle = get_test_label(res.wilcox, detailed = T)))


# Sexe ratio R et RP------------------------------------------------------------
data_fi <- unique_pat[unique_pat$Reponse == "R",]
data_fi <- rbind(data_fi, unique_pat[unique_pat$Reponse == "RP",])
data_fi <- data_fi[data_fi$numero_patient_victor != "54",]

data_fi<- data_fi%>% 
  reorder_levels(Reponse, order = c( "R", "RP"))  ## "T",


dat <- data.frame("M" = c(4,1),
                  "F" = c(18,7),
                  row.names = c("R", "RP"),
                  stringsAsFactors = F)

mosaicplot(dat, 
           main = "Répartion des patients en fonction de leur groupe et de leur type de symptôme",
           color = T)

test <- fisher.test(dat)
test
test$p.value

# Visualisation
x <- c()
for (row in rownames(dat)) {
  for (col in colnames(dat)) {
    x <- rbind(x, matrix(rep(c(row, col), dat[row, col]), ncol = 2, byrow = TRUE))
  }
}
df <- as.data.frame(x)
colnames(df) <- c("Group", "Sex")
df

library(ggstatsplot)
ggbarstats(df, 
           Group, 
           Sex,
           results.subtitle = FALSE,
           subtitle = paste0("Fisher's exact test", ", p-value = ",ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))),
           title = "Sex ratio in groups"
)



