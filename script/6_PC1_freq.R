# 1 Library---------------------------------------------------------------------
library(tidyverse)
library(ggpubr)
library(rstatix)
library(readxl)
library(ggplot2)
library(dplyr)
library(stringr)
library(cowplot)
library(openxlsx)

# 2 Import data-----------------------------------------------------------------
rm(list = ls())
PC1_VT1_T <- read.table("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/PC1_VT1.txt")
PC1_VT2_T <- read.table("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/PC1_VT2.txt")
data_V1 <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",6)
data_6m <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",8)
data_V4_CD4_8_Agspe <- read_xlsx("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/nano_freq_V4_T_cell_omiq.xlsx",2)
data_V4_T_CD4_8<- read_xlsx("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/nano_freq_V4_T_cell_omiq.xlsx",1)

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

### visalisation & save ###

for (p in list_plot) {
  print(p)
} # affiche les plots


cowplot::plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/main_subset_VT1.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

### 3.1.1 VT1 correlation-------------------------------------------------------
### Calcule des correlations & save ###
matcor <- cbind(data_fi$`Age à S1`, data_fi[11:50])
colnames(matcor)[1] <- "age"
mcor <- t(as.data.frame(cor(data_fi$SelectPCA.PC1, matcor, method= "pearson")))
write.table(mcor, file = "result/correlation_main_subset_VT1.txt")

### correlation with p-value ###
namecor <- names(matcor)
list_plot <- list()
for (i in namecor) {
  pi<-cor.test(data_fi$SelectPCA.PC1, matcor[[i]], method= "pearson")
  list_plot[[i]] <- pi
}
for (p in list_plot) {
  print(p)
}
dat <-  NULL
for (i in list_plot) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- colnames(matcor)
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_main_subset_VT1.xlsx", rowNames = T)



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
colnames(matcor)[1] <- "age"
mcor <- t(as.data.frame(cor(data_fi$SelectPCA.PC1, matcor)))
write.table(mcor, file = "result/correlation_main_subset_VT2.txt")

### 3.1.2 VT2 correlation-------------------------------------------------------
### correlation with p-value ###
namecor <- names(matcor)
list_plot <- list()
for (i in namecor) {
  pi<-cor.test(data_fi$SelectPCA.PC1, matcor[[i]], method= "pearson")
  list_plot[[i]] <- pi
}
for (p in list_plot) {
  print(p)
}
dat <-  NULL
for (i in list_plot) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- colnames(matcor)
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_main_subset_VT2.xlsx", rowNames = T)


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
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$REPONSE <- data_fi$REPONSE
  data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  # plot 
  p_i <- ggplot(
    data_clean, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
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

### 4.1.1 VT1 correlation-------------------------------------------------------
### Calcule correlation ###

list_cor <- NULL
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- t(as.data.frame(cor(data_clean$SelectPCA.PC1, data_clean[1])))
  list_cor[namevar_i] <-p_i
}

list_cor <- as.data.frame(list_cor)
write.table(list_cor, file = "result/correlation_6_months_VT1.txt")


### correlation with p-value ###

list_cor_pval <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- cor.test(data_clean$SelectPCA.PC1, data_clean[[namevar_i]], method= "pearson")
  list_cor_pval[[namevar_i]] <- p_i
}

for (p in list_cor_pval) {
  print(p)
}

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_6_months_VT1.xlsx", rowNames = T)

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
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$REPONSE <- data_fi$REPONSE
  data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  # plot 
  p_i <- ggplot(
    data_clean, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
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

### 4.1.2 VT2 correlation-------------------------------------------------------
### Calcule correlation ###

list_cor <- NULL
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- t(as.data.frame(cor(data_clean$SelectPCA.PC1, data_clean[1])))
  list_cor[namevar_i] <-p_i
}

list_cor <- as.data.frame(list_cor)
write.table(list_cor, file = "result/correlation_6_months_VT2.txt")

### correlation with p-value ###

list_cor_pval <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- cor.test(data_clean$SelectPCA.PC1, data_clean[[namevar_i]], method= "pearson")
  list_cor_pval[[namevar_i]] <- p_i
}

for (p in list_cor_pval) {
  print(p)
}

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_6_months_VT2.xlsx", rowNames = T)

# 5 Frequence cell V4-----------------------------------------------------------
## 5.1 V4 CD4 CD8 Ag spe--------------------------------------------------------
### 5.1.1 VT1 soit separation NR------------------------------------------------

data_fi <- inner_join(data_V4_CD4_8_Agspe, PC1_VT1_T, by = c("numero_patient" = "numero_patient",
                                                 "Responder"="REPONSE")) %>%
  # filter(REPONSE != "T") %>%
  rename(REPONSE = "Responder") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

namevars_to_plot <- names(data_fi[9:40])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$REPONSE <- data_fi$REPONSE
  data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  # plot 
  p_i <- ggplot(
    data_clean, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
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
ggsave(filename = "result/V4_CD4_CD8_Ag_spe_VT1.pdf", plot = p_all, width = 20, height = 35, limitsize = F )

#### 5.1.1.1 VT1 correlation----------------------------------------------------
### Calcule correlation ###

list_cor <- NULL
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- t(as.data.frame(cor(data_clean$SelectPCA.PC1, data_clean[1])))
  list_cor[namevar_i] <-p_i
}

list_cor <- as.data.frame(list_cor)
write.table(list_cor, file = "result/correlation_V4_CD4_CD8_Ag_spe_VT1.txt")


### correlation with p-value ###

list_cor_pval <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- cor.test(data_clean$SelectPCA.PC1, data_clean[[namevar_i]], method= "pearson")
  list_cor_pval[[namevar_i]] <- p_i
}

for (p in list_cor_pval) {
  print(p)
}

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_V4_CD4_CD8_Ag_spe_VT1.xlsx", rowNames = T)

### 5.1.2 VT2 soit separation RP------------------------------------------------

data_fi <- inner_join(data_V4_CD4_8_Agspe, PC1_VT2_T, by = c("numero_patient" = "numero_patient",
                                                            "Responder"="REPONSE")) %>%
  # filter(REPONSE != "T") %>%
  rename(REPONSE = "Responder") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

namevars_to_plot <- names(data_fi[9:40])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$REPONSE <- data_fi$REPONSE
  data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  # plot 
  p_i <- ggplot(
    data_clean, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
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
ggsave(filename = "result/V4_CD4_CD8_Ag_spe_VT2.pdf", plot = p_all, width = 20, height = 35, limitsize = F )

#### 5.1.2.1 VT2 correlation----------------------------------------------------
### Calcule correlation ###

list_cor <- NULL
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- t(as.data.frame(cor(data_clean$SelectPCA.PC1, data_clean[1])))
  list_cor[namevar_i] <-p_i
}

list_cor <- as.data.frame(list_cor)
write.table(list_cor, file = "result/correlation_V4_CD4_CD8_Ag_spe_VT2.txt")


### correlation with p-value ###

list_cor_pval <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- cor.test(data_clean$SelectPCA.PC1, data_clean[[namevar_i]], method= "pearson")
  list_cor_pval[[namevar_i]] <- p_i
}

for (p in list_cor_pval) {
  print(p)
}

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_V4_CD4_CD8_Ag_spe_VT2.xlsx", rowNames = T)

## 5.2 V4 T CD4 CD8 totaux------------------------------------------------------
data_V4_T_CD4_8[,9:62]<- lapply(data_V4_T_CD4_8[,9:62], as.numeric) 
data_True_false <- data_V4_T_CD4_8$Stimulation %in% "no stim"
data_V4_T_CD4_8_no_stim <- data_V4_T_CD4_8[data_True_false==T,]
data_V4_T_CD4_8_stim <- data_V4_T_CD4_8[data_True_false==F,]

### 5.2.1 V4 T CD4 CD8 avec STIM------------------------------------------------
#### 5.2.1.1 VT1 soit separation NR---------------------------------------------

data_fi <- inner_join(data_V4_T_CD4_8_stim, PC1_VT1_T, by = c("numero_patient" = "numero_patient",
                                                             "Responder"="REPONSE")) %>%
  # filter(REPONSE != "T") %>%
  rename(REPONSE = "Responder") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

namevars_to_plot <- names(data_fi[9:62])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$REPONSE <- data_fi$REPONSE
  data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  # plot 
  p_i <- ggplot(
    data_clean, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
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
ggsave(filename = "result/V4_T_CD4_CD8_stim_totaux_VT1.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

##### 5.2.1.1.1 VT1 correlat----------------------------------------------------
### Calcule correlation ###

list_cor <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- t(as.data.frame(cor(data_clean$SelectPCA.PC1, data_clean[1])))
  list_cor[namevar_i] <-p_i
}

list_cor <- as.data.frame(list_cor)
write.table(list_cor, file = "result/correlation_V4_T_CD4_CD8_stim_totaux_VT1.txt")


### correlation with p-value ###

list_cor_pval <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- cor.test(data_clean$SelectPCA.PC1, data_clean[[namevar_i]], method= "pearson")
  list_cor_pval[[namevar_i]] <- p_i
}

for (p in list_cor_pval) {
  print(p)
}

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_V4_T_CD4_CD8_stim_totaux_VT1.xlsx", rowNames = T)

#### 5.2.1.2 VT2 soit separation RP------------------------------------------------

data_fi <- inner_join(data_V4_T_CD4_8_stim, PC1_VT2_T, by = c("numero_patient" = "numero_patient",
                                                             "Responder"="REPONSE")) %>%
  # filter(REPONSE != "T") %>%
  rename(REPONSE = "Responder") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

namevars_to_plot <- names(data_fi[9:62])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$REPONSE <- data_fi$REPONSE
  data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  # plot 
  p_i <- ggplot(
    data_clean, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
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
ggsave(filename = "result/V4_T_CD4_CD8_stim_totaux_VT2.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

##### 5.2.2.1 VT2 correlation----------------------------------------------------
### Calcule correlation ###

list_cor <- NULL
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- t(as.data.frame(cor(data_clean$SelectPCA.PC1, data_clean[1])))
  list_cor[namevar_i] <-p_i
}

list_cor <- as.data.frame(list_cor)
write.table(list_cor, file = "result/correlation_V4_T_CD4_CD8_stim_totaux_VT2.txt")


### correlation with p-value ###

list_cor_pval <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- cor.test(data_clean$SelectPCA.PC1, data_clean[[namevar_i]], method= "pearson")
  list_cor_pval[[namevar_i]] <- p_i
}

for (p in list_cor_pval) {
  print(p)
}

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_V4_T_CD4_CD8_stim_totaux_VT2.xlsx", rowNames = T)

### 5.2.2 V4 T CD4 CD8 avec NO STIM------------------------------------------------
#### 5.2.2.1 VT1 soit separation NR---------------------------------------------

data_fi <- inner_join(data_V4_T_CD4_8_no_stim, PC1_VT1_T, by = c("numero_patient" = "numero_patient",
                                                              "Responder"="REPONSE")) %>%
  # filter(REPONSE != "T") %>%
  rename(REPONSE = "Responder") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

namevars_to_plot <- names(data_fi[9:62])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$REPONSE <- data_fi$REPONSE
  data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  # plot 
  p_i <- ggplot(
    data_clean, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
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
ggsave(filename = "result/V4_T_CD4_CD8_no_stim_totaux_VT1.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

##### 5.2.2.1.1 VT1 correlat----------------------------------------------------
### Calcule correlation ###

list_cor <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- t(as.data.frame(cor(data_clean$SelectPCA.PC1, data_clean[1])))
  list_cor[namevar_i] <-p_i
}

list_cor <- as.data.frame(list_cor)
write.table(list_cor, file = "result/correlation_V4_T_CD4_CD8_no_stim_totaux_VT1.txt")


### correlation with p-value ###

list_cor_pval <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- cor.test(data_clean$SelectPCA.PC1, data_clean[[namevar_i]], method= "pearson")
  list_cor_pval[[namevar_i]] <- p_i
}

for (p in list_cor_pval) {
  print(p)
}

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_V4_T_CD4_CD8_no_stim_totaux_VT1.xlsx", rowNames = T)

#### 5.2.2.2 VT2 soit separation RP------------------------------------------------

data_fi <- inner_join(data_V4_T_CD4_8_no_stim, PC1_VT2_T, by = c("numero_patient" = "numero_patient",
                                                              "Responder"="REPONSE")) %>%
  # filter(REPONSE != "T") %>%
  rename(REPONSE = "Responder") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

namevars_to_plot <- names(data_fi[9:62])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$REPONSE <- data_fi$REPONSE
  data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  # plot 
  p_i <- ggplot(
    data_clean, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
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
ggsave(filename = "result/V4_T_CD4_CD8_no_stim_totaux_VT2.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

##### 5.2.2.2.1 VT2 correlation----------------------------------------------------
### Calcule correlation ###

list_cor <- NULL
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- t(as.data.frame(cor(data_clean$SelectPCA.PC1, data_clean[1])))
  list_cor[namevar_i] <-p_i
}

list_cor <- as.data.frame(list_cor)
write.table(list_cor, file = "result/correlation_V4_T_CD4_CD8_no_stim_totaux_VT2.txt")

### correlation with p-value ###

list_cor_pval <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  # data_clean$REPONSE <- data_fi$REPONSE
  # data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  
  p_i <- cor.test(data_clean$SelectPCA.PC1, data_clean[[namevar_i]], method= "pearson")
  list_cor_pval[[namevar_i]] <- p_i
}

for (p in list_cor_pval) {
  print(p)
}

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_V4_T_CD4_CD8_no_stim_totaux_VT2.xlsx", rowNames = T)
