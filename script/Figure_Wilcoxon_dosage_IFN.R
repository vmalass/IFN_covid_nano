Dosage IFN R et RP 

library(readxl)
library(ggplot2)
library(ggrepel)
library(rstatix)
library(ggpubr)
library(dplyr)

# Import data-------------------------------------------------------------------
rm(list = ls())

data <- read_xlsx("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/Dosage_IFN.xlsx", sheet = 2)
data_2 <- read_xlsx("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/Dosage_IFN.xlsx", sheet = 3) # avec IFN alpha

T_F <- data$Reponse %in% c("R", "RP")
data_fi <- data[T_F == T,]

T_F <- data_2$Reponse %in% c("R", "RP")
data_fi_2 <- data_2[T_F == T,]

data_fi[,1:3] <- lapply(data_fi[,1:3], as.numeric)

tab <- data_fi %>%
  dplyr::count(Reponse, new_time_point)  #numero_patient_victor

tab <- tab %>% 
  dplyr::rename("nombre echantillon" = "n")

tab_2 <- data_fi_2 %>%
  dplyr::count(Reponse, new_time_point)  #numero_patient_victor

tab_2 <- tab_2 %>% 
  dplyr::rename("nombre echantillon" = "n")

# Visualisation IFN-------------------------------------------------------------
ggplot(data = data_fi, 
       aes(x = new_time_point, 
           y = `IFN_l1_IL_29`, 
           color = Reponse))+ 
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(y_position = c(57, 46, 33, 20), 
              xmin = c(0.8, 1.8, 2.8, 3.8), 
              xmax = c(1.2, 2.2, 3.2, 4.2),
              annotation = c("0.76", "0.4", "0.52", "0.74")) + 
  labs(title = "Wilcoxon test for title IFN-lambda", 
       subtitle = "VT1 6R 4RP / VT2 4R 3RP / VT3 7R 3RP / VT4 6R 2RP",
       x = "", 
       y = "IFN-lambda (pg/mL)")+
  theme_classic() + 
  theme(legend.position = c(0.9,0.8))

ggplot(data = data_fi, 
       aes(x = new_time_point, 
           y = `IFN_gamma`, 
           color = Reponse))+ 
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25)) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  geom_signif(y_position = c(1200, 450, 150, 150), 
              xmin = c(0.8, 1.8, 2.8, 3.8), 
              xmax = c(1.2, 2.2, 3.2, 4.2),
              annotation = c("0.35", "0.057", "0.18", "0.29")) + 
  labs(title = "Wilcoxon test for title IFN gamma", 
       subtitle = "VT1 6R 4RP / VT2 4R 3RP / VT3 7R 3RP / VT4 6R 2RP",
       x = "", 
       y = "IIFN gamma (pg/mL)")+
  theme_classic() + 
  theme(legend.position = c(0.9,0.8))

dplyr::count(data_fi, new_time_point, Reponse)

ggplot(data = data_fi, 
       aes(x = new_time_point, 
           y = `IFN_beta`, 
           color = Reponse))+ 
  geom_boxplot() +
  geom_point(position = position_dodge(width=0.75)) +
  scale_color_manual(breaks = c("R", "RP"),
                     values = c("gray50","#CB2027")) +
  labs(title = "Wilcoxon test for title IFN beta", 
       x = "Time point", 
       y = "IIFN beta (pg/mL)")+
  theme_classic() + 
  theme(legend.position = c(0.9,0.8))


# Test de Wilcoxon--------------------------------------------------------------
vector <- c("VT1", "VT2", "VT3", "VT4")
for (time_point in vector) {
  ### Table ###
  table <- data_fi[data_fi$new_time_point == time_point,]
  table <- table %>% 
    reorder_levels(Reponse, order = c( "R", "RP"))
  ### test stat ###
  res.wilcox <- table %>%
    wilcox_test(formula(paste0("IFN_l1_IL_29"," ~ Reponse"))) %>%
    add_significance()

  ### Visualization ###
  res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
  print(ggboxplot(table, 
                  title = paste0("Wilcoxon test for title IFN I1 IL29 in ", time_point),
                  x = "Reponse", 
                  y = paste0("IFN_l1_IL_29"), 
                  ylab = "quantity of IFN (pg/mL)",
                  xlab = "Groups", 
                  add = "point",
                  color = "Reponse",
                  palette = c("gray50", "#CB2027"),
                  legend = c(0.9,0.95)) + 
          stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
          labs(subtitle = get_test_label(res.wilcox, detailed = T)))
}

for (time_point in vector) {
  ### Table ###
  table <- data_fi[data_fi$new_time_point == time_point,]
  table <- table %>% 
    reorder_levels(Reponse, order = c( "R", "RP"))
  ### test stat ###
  res.wilcox <- table %>%
    wilcox_test(formula(paste0("IFN_gamma"," ~ Reponse"))) %>%
    add_significance()
  
  ### Visualization ###
  res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
  print(ggboxplot(table, 
                  title = paste0("Wilcoxon test for title IFN gamma in ", time_point),
                  x = "Reponse", 
                  y = paste0("IFN_gamma"), 
                  ylab = "quantity of IFN (pg/mL)",
                  xlab = "Groups", 
                  add = "point",
                  color = "Reponse",
                  palette = c("gray50", "#CB2027"),
                  legend = c(0.9,0.3)) + 
          stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
          labs(subtitle = get_test_label(res.wilcox, detailed = T)))
}

for (time_point in vector) {
  ### Table ###
  table <- data_fi[data_fi$new_time_point == time_point,]
  table <- table %>% 
    reorder_levels(Reponse, order = c( "R", "RP"))
  ### test stat ###
  res.wilcox <- table %>%
    wilcox_test(formula(paste0("IFN_beta"," ~ Reponse"))) %>%
    add_significance()
  
  ### Visualization ###
  res.wilcox <- res.wilcox %>% add_xy_position(x = "Reponse")
  print(ggboxplot(table, 
                  title = paste0("Wilcoxon test for title IFN Beta in ", time_point),
                  x = "Reponse", 
                  y = paste0("IFN_beta"), 
                  ylab = "quantity of IFN (pg/mL)",
                  xlab = "Groups", 
                  add = "point",
                  color = "Reponse",
                  palette = c("gray50", "#CB2027"),
                  legend = c(0.9,0.95)) + 
          stat_pvalue_manual(res.wilcox, tip.length = 0.03, hide.ns = T) +
          labs(subtitle = get_test_label(res.wilcox, detailed = T)))
}
