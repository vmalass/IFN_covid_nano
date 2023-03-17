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
data_V1 <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",6)
data_V1_non_sup <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",7)
# data_V4_CD4_8_Agspe <- read_xlsx("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/nano_freq_V4_T_cell_omiq.xlsx",2)
# data_V4_T_CD4_8<- read_xlsx("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/nano_freq_V4_T_cell_omiq.xlsx",1)
# data_6m <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",8)
PC1_V1 <- read.table("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/PC1_V1_gene_DE_VT1vsT_mat1.3.txt")
load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg

# 3 Frequence cell V1 supervisé-------------------------------------------------
## 3.1 Freq vs jour de prelevement----------------------------------------------
ma <- mat_pat_clean$time_point %in% "V1"
mat <- mat_pat_clean[ma == T,]
mat[,1:736] <- NULL
mat$numero_patient <- as.character(mat$numero_patient)
mat$numero_patient <- as.numeric(mat$numero_patient)

data_fi <- inner_join(data_V1, mat, by = c("numero_patient_nanostring" = "numero_patient",
                                               "REPONSE"="REPONSE")) %>%
  # filter(REPONSE != "T") %>%
  rename(numero_patient = "numero_patient_nanostring") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

ma <- data_V1$REPONSE %in% "T"
mat <- data_V1[ma == T,]
mat$jours_prelevement <- 0
mat$jours_prelevement <- as.factor(mat$jours_prelevement)
mat$label_patient <- "T"
data_fi <- bind_rows(data_fi, mat )

data_fi$jours_prelevement <- as.character(data_fi$jours_prelevement)
data_fi$jours_prelevement <- as.numeric(data_fi$jours_prelevement)

setdiff(data_V1$numero_patient_nanostring, data_fi$numero_patient)
## ne manque pas de patient

namevars_to_plot <- names(data_fi[8:56])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  p_i <- ggplot(
    data_fi, aes(x = jours_prelevement , y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
    geom_point()+
    ggrepel::geom_text_repel(# nudge_x=0.6,
                             # nudge_y=0.15,
                             show.legend = F,
                             max.overlaps  = Inf)+
    # scale_color_manual(breaks = c("NR","R","RP","T"),
    #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
    scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                       values = c("darkorange", "#DDCC77","cornflowerblue",
                                  "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) 
  labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i
  
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/main_subset_V1_sup_jour_prelevement_T.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

### 3.1.1 VT1 correlation-------------------------------------------------------
### Calcule des correlations & save ###
ma <- data_fi$REPONSE %in% "T"
data_fi <- data_fi[ma == F,]
matcor <- cbind(data_fi$`Age_à_S1`, data_fi[8:56])
colnames(matcor)[1] <- "age"

### correlation with p-value ###
namecor <- names(matcor)
list_plot <- list()
for (i in namecor) {
  pi<-cor.test(data_fi$jours_prelevement, matcor[[i]], method= "pearson")
  list_plot[[i]] <- pi
}


dat <-  NULL
for (i in list_plot) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- colnames(matcor)
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_main_subset_V1_sup_jour_prelevement.xlsx", rowNames = T)

## 3.2 Freq vs PCA--------------------------------------------------------------
data_fi <- inner_join(data_V1, PC1_V1, by = c("numero_patient_nanostring" = "numero_patient",
                                           "REPONSE"="REPONSE")) %>%
  filter(REPONSE != "T") %>%
  rename(numero_patient = "numero_patient_nanostring") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

setdiff(data_V1$numero_patient_nanostring, data_fi$numero_patient)
## ne manque pas de patient

namevars_to_plot <- names(data_fi[8:56])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  p_i <- ggplot(
    data_fi, aes(x = SelectPCA.PC1 , y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
    geom_point()+
    ggrepel::geom_text_repel(# nudge_x=0.6,
      # nudge_y=0.15,
      show.legend = F,
      max.overlaps  = Inf)+
    # scale_color_manual(breaks = c("NR","R","RP","T"),
    #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
    scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                       values = c("darkorange", "#DDCC77","cornflowerblue",
                                  "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) 
  labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i
  
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/main_subset_V1_sup_PC1.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

### 3.2.1 VT1 correlation-------------------------------------------------------
### Calcule des correlations & save ###
matcor <- cbind(data_fi$`Age_à_S1`, data_fi[8:56])
colnames(matcor)[1] <- "age"

### correlation with p-value ###
namecor <- names(matcor)
list_plot <- list()
for (i in namecor) {
  pi<-cor.test(data_fi$SelectPCA.PC1, matcor[[i]], method= "pearson")
  list_plot[[i]] <- pi
}


dat <-  NULL
for (i in list_plot) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- colnames(matcor)
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_main_subset_V1_sup_PC1.xlsx", rowNames = T)





# 3' Frequence cell V1 non supervisé-------------------------------------------------
## 3'.1 Freq vs jour de prelevement----------------------------------------------
ma <- mat_pat_clean$time_point %in% "V1"
mat <- mat_pat_clean[ma == T,]
mat[,1:736] <- NULL
mat$numero_patient <- as.character(mat$numero_patient)
mat$numero_patient <- as.numeric(mat$numero_patient)

data_fi <- inner_join(data_V1_non_sup, mat, by = c("numero_patient_nanostring" = "numero_patient",
                                                   "REPONSE"="REPONSE")) %>%
  # filter(REPONSE != "T") %>%
  rename(numero_patient = "numero_patient_nanostring") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

ma <- data_V1_non_sup$REPONSE %in% "T"
mat <- data_V1_non_sup[ma == T,]
mat$jours_prelevement <- 0
mat$jours_prelevement <- as.factor(mat$jours_prelevement)
mat$label_patient <- "T"
data_fi <- bind_rows(data_fi, mat )

data_fi$jours_prelevement <- as.character(data_fi$jours_prelevement)
data_fi$jours_prelevement <- as.numeric(data_fi$jours_prelevement)

setdiff(data_V1_non_sup$numero_patient_nanostring, data_fi$numero_patient)
## ne manque pas de patient

namevars_to_plot <- names(data_fi[9:41])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  p_i <- ggplot(
    data_fi, aes(x = jours_prelevement , y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
    geom_point()+
    ggrepel::geom_text_repel(# nudge_x=0.6,
      # nudge_y=0.15,
      show.legend = F,
      max.overlaps  = Inf)+
    # scale_color_manual(breaks = c("NR","R","RP","T"),
    #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
    scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                       values = c("darkorange", "#DDCC77","cornflowerblue",
                                  "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) 
  labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i
  
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/main_subset_V1_non_sup_jour_prelevement_T.pdf", plot = p_all, width = 20, height = 35, limitsize = F )

### 3'.1.1 VT1 correlation-------------------------------------------------------
### correlation with p-value ###
ma <- data_fi$REPONSE %in% "T"
data_fi <- data_fi[ma == F,]

list_plot <- list()
for (i in namevars_to_plot) {
  pi<-cor.test(data_fi$jours_prelevement, data_fi[[i]], method= "pearson")
  list_plot[[i]] <- pi
}


dat <-  NULL
for (i in list_plot) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_main_subset_V1_non_sup_jour_prelevement.xlsx", rowNames = T)

## 3'.2 Freq vs PCA--------------------------------------------------------------
data_fi <- inner_join(data_V1_non_sup, PC1_V1, by = c("numero_patient_nanostring" = "numero_patient",
                                                      "REPONSE"="REPONSE")) %>%
  filter(REPONSE != "T") %>%
  rename(numero_patient = "numero_patient_nanostring") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

setdiff(data_V1_non_sup$numero_patient_nanostring, data_fi$numero_patient)
## ne manque pas de patient

namevars_to_plot <- names(data_fi[9:41])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  p_i <- ggplot(
    data_fi, aes(x = SelectPCA.PC1 , y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
    geom_point()+
    ggrepel::geom_text_repel(# nudge_x=0.6,
      # nudge_y=0.15,
      show.legend = F,
      max.overlaps  = Inf)+
    # scale_color_manual(breaks = c("NR","R","RP","T"),
    #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
    scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                       values = c("darkorange", "#DDCC77","cornflowerblue",
                                  "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) 
  labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i
  
}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/main_subset_V1_non_sup_PC1.pdf", plot = p_all, width = 20, height = 35, limitsize = F )

### 3'.2.1 VT1 correlation-------------------------------------------------------
### correlation with p-value ###
list_plot <- list()
for (i in namevars_to_plot) {
  pi<-cor.test(data_fi$SelectPCA.PC1, data_fi[[i]], method= "pearson")
  list_plot[[i]] <- pi
}


dat <-  NULL
for (i in list_plot) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_main_subset_V1_non_sup_PC1.xlsx", rowNames = T)










# .o
# .o
# .o
# .o



























# 4 Frequence cell V4 CD4-8 Ag spe----------------------------------------------
## 4.2 Freq vs PCA--------------------------------------------------------------
data_V4_CD4_8_Agspe <- data_V4_CD4_8_Agspe[data_V4_CD4_8_Agspe$numero_patient != 50,] ## supression du patient 50

data_fi <- inner_join(data_V4_CD4_8_Agspe, PC1_V1, by = c("numero_patient" = "numero_patient",
                                              "Responder"="REPONSE")) %>%
  filter(Responder != "T") %>%
  rename(numero_patient = "numero_patient") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

setdiff(data_V4_CD4_8_Agspe$numero_patient, data_fi$numero_patient)
##  manque 56 64 66 normale car un seul prelevement

namevars_to_plot <- names(data_fi[9:40])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  p_i <- ggplot(
    data_fi, aes(x = SelectPCA.PC1 , y = .data[[namevar_i]], color = Responder, label = label_patient)) +
    geom_point()+
    ggrepel::geom_text_repel(show.legend = F,
                             max.overlaps  = Inf)+
    scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                       values = c("darkorange", "#DDCC77","cornflowerblue",
                                  "#882255" ,"brown3", "chartreuse4", "#BBBBBB"))
  labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i

}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/V4_CD4_8_Agspe_non_supervise_sans50.pdf", plot = p_all, width = 20, height = 35, limitsize = F )

### 4.2.1 V4 correlation-------------------------------------------------------
### Calcule des correlations & save ###
list_cor_pval <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)

  p_i <- cor.test(data_clean$SelectPCA.PC1, data_clean[[namevar_i]], method= "pearson")
  list_cor_pval[[namevar_i]] <- p_i
}

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_V4_CD4_8_Agspe_non_supervise_sans50.xlsx", rowNames = T)


# 5 Frequence cell V4 CD4-8 ----------------------------------------------------
## 5.1 V4 T CD4 CD8 totaux------------------------------------------------------
data_V4_T_CD4_8[,9:62]<- lapply(data_V4_T_CD4_8[,9:62], as.numeric)

data_V4_T_CD4_8 <- data_V4_T_CD4_8[data_V4_T_CD4_8$numero_patient != 50,] ## supprime le patient 50

data_True_false <- data_V4_T_CD4_8$Stimulation %in% "no stim"
data_V4_T_CD4_8_no_stim <- data_V4_T_CD4_8[data_True_false==T,]
data_V4_T_CD4_8_stim <- data_V4_T_CD4_8[data_True_false==F,]

## 5.2 Freq vs PCA non stim-----------------------------------------------------
data_fi <- inner_join(data_V4_T_CD4_8_no_stim, PC1_V1, by = c("numero_patient" = "numero_patient",
                                                          "Responder"="REPONSE")) %>%
  filter(Responder != "T") %>%
  rename(numero_patient = "numero_patient") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

setdiff(data_V4_T_CD4_8_no_stim$numero_patient, data_fi$numero_patient)
##  manque 56 64 66 normale car un seul prelevement

namevars_to_plot <- names(data_fi[9:62])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$Responder <- data_fi$Responder
  data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  # plot
  p_i <- ggplot(
    data_clean, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = Responder, label = label_patient)) +
    geom_point()+
    ggrepel::geom_text_repel(show.legend = F,
                             max.overlaps  = Inf)+
    scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                       values = c("darkorange", "#DDCC77","cornflowerblue",
                                  "#882255" ,"brown3", "chartreuse4", "#BBBBBB"))
  labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i

}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/V4_CD4_8_non_supervise_no_stim_sans50.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

### 5.2.1 V4 non stim correlation-----------------------------------------------
### Calcule des correlations & save ###
list_cor_pval <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)

  p_i <- cor.test(data_clean$SelectPCA.PC1, data_clean[[namevar_i]], method= "pearson")
  list_cor_pval[[namevar_i]] <- p_i
}

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_V4_CD4_8_non_supervise_no_stim_sans50.xlsx", rowNames = T)


## 5.2 Freq vs PCA stim---------------------------------------------------------
data_fi <- inner_join(data_V4_T_CD4_8_stim, PC1_V1, by = c("numero_patient" = "numero_patient",
                                                              "Responder"="REPONSE")) %>%
  filter(Responder != "T") %>%
  rename(numero_patient = "numero_patient") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

setdiff(data_V4_T_CD4_8_stim$numero_patient, data_fi$numero_patient)
##  manque 56 64 66 normale car un seul prelevement

namevars_to_plot <- names(data_fi[9:62])

list_plot <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  # suppression des NA dans la variable
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$Responder <- data_fi$Responder
  data_clean$label_patient <- data_fi$label_patient
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)
  # plot
  p_i <- ggplot(
    data_clean, aes(x = SelectPCA.PC1, y = .data[[namevar_i]], color = Responder, label = label_patient)) +
    geom_point()+
    ggrepel::geom_text_repel(show.legend = F,
                             max.overlaps  = Inf)+
    scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                       values = c("darkorange", "#DDCC77","cornflowerblue",
                                  "#882255" ,"brown3", "chartreuse4", "#BBBBBB"))
  labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i

}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/V4_CD4_8_non_supervise_stim_sans50.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

### 5.2.1 V4 stim correlation---------------------------------------------------
### Calcule des correlations & save ###
list_cor_pval <- list()
for (namevar_i in namevars_to_plot) {
  print(namevar_i)
  data_clean <- as.data.frame(data_fi[namevar_i])
  data_clean$SelectPCA.PC1 <- data_fi$SelectPCA.PC1
  data_clean <- na.omit(data_clean)

  p_i <- cor.test(data_clean$SelectPCA.PC1, data_clean[[namevar_i]], method= "pearson")
  list_cor_pval[[namevar_i]] <- p_i
}

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_V4_CD4_8_non_supervise_stim_sans50.xlsx", rowNames = T)


# 6 Frequence cell 6 months-----------------------------------------------------
## 6.1 Freq vs jour de prelevement----------------------------------------------
# ma <- mat_pat_clean$time_point %in% "V1"
# mat <- mat_pat_clean[ma == T,]
# mat[,1:736] <- NULL
# mat$numero_patient <- as.character(mat$numero_patient)
# mat$numero_patient <- as.numeric(mat$numero_patient)
#
# data_fi <- inner_join(data_6m, mat, by = c("numero_patient_nanostring" = "numero_patient",
#                                            "REPONSE"="REPONSE")) %>%
#   # filter(REPONSE != "T") %>%
#   rename(numero_patient = "numero_patient_nanostring") %>%
#   mutate(label_patient = as.character(numero_patient)) %>%
#   arrange(numero_patient)
#
# ma <- data_6m$REPONSE %in% "T"
# mat <- data_6m[ma == T,]
# mat$jours_prelevement <- 0
# mat$jours_prelevement <- as.factor(mat$jours_prelevement)
# mat$label_patient <- "T"
# data_fi <- bind_rows(data_fi, mat )
#
# data_fi$jours_prelevement <- as.character(data_fi$jours_prelevement)
# data_fi$jours_prelevement <- as.numeric(data_fi$jours_prelevement)
# data_fi <- data_fi[data_fi$label_patient !=50,] ## Suprime le patient 50
#
# setdiff(data_6m$numero_patient_nanostring, data_fi$numero_patient)
# ## ne manque pas de patient
#
# namevars_to_plot <- names(data_fi[14:99])
#
# list_plot <- list()
# for (namevar_i in namevars_to_plot) {
#   print(namevar_i)
#   # suppression des NA dans la variable
#   data_clean <- as.data.frame(data_fi[namevar_i])
#   data_clean$REPONSE <- data_fi$REPONSE
#   data_clean$label_patient <- data_fi$label_patient
#   data_clean$jours_prelevement <- data_fi$jours_prelevement
#   data_clean <- na.omit(data_clean)
#   # plot
#   p_i <- ggplot(
#     data_clean, aes(x = jours_prelevement, y = .data[[namevar_i]], color = REPONSE, label = label_patient)) +
#     geom_point()+
#     ggrepel::geom_text_repel(# nudge_x=0.6,
#       # nudge_y=0.15,
#       show.legend = F,
#       max.overlaps  = Inf)+
#     # scale_color_manual(breaks = c("NR","R","RP","T"),
#     #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
#     scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
#                        values = c("darkorange", "#DDCC77","cornflowerblue",
#                                   "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
#     labs(title=namevar_i) +
#     theme_bw()
#   list_plot[[namevar_i]] <- p_i
#
# }
#
# plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page
#
# p_all = plot_grid(plotlist = list_plot, ncol = 3)
# ggsave(filename = "result/6_months_V1_jour_prelevement_T_sans50.pdf", plot = p_all, width = 20, height = 100, limitsize = F )
#
# ### 6.1.1 VT1 correlation-------------------------------------------------------
# ### correlation with p-value ###
# list_cor_pval <- list()
# for (namevar_i in namevars_to_plot) {
#   print(namevar_i)
#   # suppression des NA dans la variable
#   data_clean <- as.data.frame(data_fi[namevar_i])
#   # data_clean$REPONSE <- data_fi$REPONSE
#   # data_clean$label_patient <- data_fi$label_patient
#   data_clean$jours_prelevement <- data_fi$jours_prelevement
#   data_clean <- na.omit(data_clean)
#
#   p_i <- cor.test(data_clean$jours_prelevement, data_clean[[namevar_i]], method= "pearson")
#   list_cor_pval[[namevar_i]] <- p_i
# }
#
# dat <-  NULL
# for (i in list_cor_pval) {
#   test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
#   dat <- rbind(dat, test)
# }
# rownames(dat) <- namevars_to_plot
# colnames(dat) <- c("p-val", "cor")
# write.xlsx(dat, file = "result/correlation_pvalue_6_months_V1_jour_prelevement_T.xlsx", rowNames = T)

## 6.2 Freq vs PCA--------------------------------------------------------------
data_fi <- inner_join(data_6m, PC1_V1, by = c("numero_patient_nanostring" = "numero_patient",
                                              "REPONSE"="REPONSE")) %>%
  filter(REPONSE != "T") %>%
  rename(numero_patient = "numero_patient_nanostring") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

data_fi <- data_fi[data_fi$label_patient !=50,] ## Suprime le patient 50

setdiff(data_6m$numero_patient_nanostring, data_fi$numero_patient)
## ne manque pas de patient

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
    ggrepel::geom_text_repel(# nudge_x=0.6,
      # nudge_y=0.15,
      show.legend = F,
      max.overlaps  = Inf)+
    # scale_color_manual(breaks = c("NR","R","RP","T"),
    #                    values = c("darkorange","cornflowerblue","brown3","chartreuse4"))+
    scale_color_manual(breaks = c("NR", "NR-", "R", "RP-", "RP", "T", "A"),
                       values = c("darkorange", "#DDCC77","cornflowerblue",
                                  "#882255" ,"brown3", "chartreuse4", "#BBBBBB")) +
    labs(title=namevar_i) +
    theme_bw()
  list_plot[[namevar_i]] <- p_i

}

plot_grid(plotlist = list_plot[1:6], ncol = 2) # pour orga ne nombre de plot sur la page

p_all = plot_grid(plotlist = list_plot, ncol = 3)
ggsave(filename = "result/6_months_V1_APC_sans50.pdf", plot = p_all, width = 20, height = 100, limitsize = F )

### 6.2.1 6 months correlation-------------------------------------------------------
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

dat <-  NULL
for (i in list_cor_pval) {
  test <- cbind(as.data.frame(i$p.value),as.data.frame(i$estimate))
  dat <- rbind(dat, test)
}
rownames(dat) <- namevars_to_plot
colnames(dat) <- c("p-val", "cor")
write.xlsx(dat, file = "result/correlation_pvalue_6_months_V1_APC_sans50.xlsx", rowNames = T)

