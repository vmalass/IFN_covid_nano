# 1 Library---------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)

# 2 Import data-----------------------------------------------------------------
rm(list = ls())
PC1_VT1_T <- read.table("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/PC1_VT1.txt")
PC1_VT2_T <- read.table("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/PC1_VT2.txt")
data_V1 <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",6)
data_6m <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",8)

# 3 Frequence cell V1-----------------------------------------------------------
## 3.1 VT1 soit separation NR---------------------------------------------------
data_fi <- inner_join(data_V1, PC1_VT1_T, by = c("numero_patient_nanostring" = "numero_patient",
                                                 "REPONSE"="REPONSE")) %>%
  filter(REPONSE != "T") %>%
  rename(numero_patient = "numero_patient_nanostring") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

namevars_to_plot <- names(data_fi[8:50])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  p_i <- ggplot(
    data_fi, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
    geom_point()+
    ggrepel::geom_text_repel(nudge_x=0.6,
                             # nudge_y=0.15,
                             show.legend = F,
                             max.overlaps  = Inf)+
    scale_color_manual(breaks = c("NR","R","RP"),
                       values = c("darkorange","cornflowerblue","brown3"))+
    labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i
  
}

# different moyen d'afficher tes graphes avec 2 libraries : cowplot et patchwork

for (p in list_plot) {
  print(p)
}


cowplot::plot_grid(plotlist = list_plot[1:6], ncol = 2)

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/main_subset_VT1.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

matcor <- cbind(data_fi$`Age à S1`, data_fi[11:50])
mcor <- t(as.data.frame(cor(data_fi$SelectPCA.PC1, matcor)))
write.table(mcor, file = "result/correlation_main_subset_VT1.txt")


# library(patchwork)
# p_all2 <- wrap_plots(list_plot) + plot_layout(guides = "collect" ) & theme(legend.position = 'top')
# ggsave(filename = "all_plots2.pdf", plot = p_all2, width = 20, height = 30, limitsize = F )

## 3.2 VT2 soit separation RP---------------------------------------------------
data_fi <- inner_join(data_V1, PC1_VT2_T, by = c("numero_patient_nanostring" = "numero_patient",
                                                 "REPONSE"="REPONSE")) %>%
  filter(REPONSE != "T") %>%
  rename(numero_patient = "numero_patient_nanostring") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

namevars_to_plot <- names(data_fi[8:50])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  p_i <- ggplot(
    data_fi, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
    geom_point()+
    ggrepel::geom_text_repel(nudge_x=0.6,
                             # nudge_y=0.15,
                             show.legend = F,
                             max.overlaps  = Inf)+
    scale_color_manual(breaks = c("NR","R","RP"),
                       values = c("darkorange","cornflowerblue","brown3"))+
    labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i
  
}

# different moyen d'afficher tes graphes avec 2 libraries : cowplot et patchwork

for (p in list_plot) {
  print(p)
}


cowplot::plot_grid(plotlist = list_plot[1:6], ncol = 2)

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/main_subset_VT2.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

matcor <- cbind(data_fi$`Age à S1`, data_fi[11:50])
mcor <- t(as.data.frame(cor(data_fi$SelectPCA.PC1, matcor)))
write.table(mcor, file = "result/correlation_main_subset_VT2.txt")

# 4 Frequence cell 6 mois-------------------------------------------------------
## 4.1 VT1 soit separation NR---------------------------------------------------

data_fi <- inner_join(data_6m, PC1_VT1_T, by = c("numero_patient_nanostring" = "numero_patient",
                                                 "REPONSE"="REPONSE")) %>%
  filter(REPONSE != "T") %>%
  rename(numero_patient = "numero_patient_nanostring") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

namevars_to_plot <- names(data_fi[14:99])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  p_i <- ggplot(
    data_fi, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
    geom_point()+
    ggrepel::geom_text_repel(nudge_x=0.6,
                             # nudge_y=0.15,
                             show.legend = F,
                             max.overlaps  = Inf)+
    scale_color_manual(breaks = c("NR","R","RP"),
                       values = c("darkorange","cornflowerblue","brown3"))+
    labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i
  
}

# different moyen d'afficher tes graphes avec 2 libraries : cowplot et patchwork

for (p in list_plot) {
  print(p)
}


cowplot::plot_grid(plotlist = list_plot[1:6], ncol = 2)

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/6_months_VT1.pdf", plot = p_all, width = 20, height = 100, limitsize = F )

mcor <- t(as.data.frame(cor(data_fi$SelectPCA.PC1, data_fi[14:99])))
write.table(mcor, file = "result/correlation_6_months_VT1.txt")

## 4.2 VT2 soit separation RP---------------------------------------------------

data_fi <- inner_join(data_6m, PC1_VT2_T, by = c("numero_patient_nanostring" = "numero_patient",
                                                 "REPONSE"="REPONSE")) %>%
  filter(REPONSE != "T") %>%
  rename(numero_patient = "numero_patient_nanostring") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

namevars_to_plot <- names(data_fi[14:99])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  p_i <- ggplot(
    data_fi, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
    geom_point()+
    ggrepel::geom_text_repel(nudge_x=0.6,
                             # nudge_y=0.15,
                             show.legend = F,
                             max.overlaps  = Inf)+
    scale_color_manual(breaks = c("NR","R","RP"),
                       values = c("darkorange","cornflowerblue","brown3"))+
    labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i
  
}

# different moyen d'afficher tes graphes avec 2 libraries : cowplot et patchwork

for (p in list_plot) {
  print(p)
}


cowplot::plot_grid(plotlist = list_plot[1:6], ncol = 2)

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/6_months_VT2.pdf", plot = p_all, width = 20, height = 100, limitsize = F )

mcor <- t(as.data.frame(cor(data_fi$SelectPCA.PC1, data_fi[14:99])))
write.table(mcor, file = "result/correlation_6_months_VT2.txt")


















