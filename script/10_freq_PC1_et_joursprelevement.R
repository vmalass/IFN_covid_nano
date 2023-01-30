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
data_6m <- read_xlsx("/Users/victor/Documents/JM/NanoString/Data_brut/Covid_IFN-avec_data_V5.xlsx",8)
PC1_V1 <- read.table("/Users/victor/Documents/JM/NanoString/IFN_covid_nano/data/PC1_V1_gene_DE_NRvsR.txt")
load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg

# 3 Frequence cell V1-----------------------------------------------------------
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

namevars_to_plot <- names(data_fi[8:50])

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
ggsave(filename = "result/main_subset_V1_jour_prelevement_T.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

### 3.1.1 VT1 correlation-------------------------------------------------------
### Calcule des correlations & save ###
matcor <- cbind(data_fi$`Age à S1`, data_fi[11:50])
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
write.xlsx(dat, file = "result/correlation_pvalue_main_subset_V1_jour_prelevement_T.xlsx", rowNames = T)

## 3.2 Freq vs PCA--------------------------------------------------------------
data_fi <- inner_join(data_V1, PC1_V1, by = c("numero_patient_nanostring" = "numero_patient",
                                           "REPONSE"="REPONSE")) %>%
  filter(REPONSE != "T") %>%
  rename(numero_patient = "numero_patient_nanostring") %>%
  mutate(label_patient = as.character(numero_patient)) %>%
  arrange(numero_patient)

setdiff(data_V1$numero_patient_nanostring, data_fi$numero_patient)
## ne manque pas de patient

namevars_to_plot <- names(data_fi[8:50])

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
ggsave(filename = "result/main_subset_V1.pdf", plot = p_all, width = 20, height = 50, limitsize = F )

### 3.2.1 VT1 correlation-------------------------------------------------------
### Calcule des correlations & save ###
matcor <- cbind(data_fi$`Age à S1`, data_fi[11:50])
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
write.xlsx(dat, file = "result/correlation_pvalue_main_subset_V1.xlsx", rowNames = T)

# 4 Frequence cell 6 months-----------------------------------------------------
## 4.1 Freq vs jour de prelevement----------------------------------------------
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
# ### 4.1.1 VT1 correlation-------------------------------------------------------
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

## 4.2 Freq vs PCA--------------------------------------------------------------
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

### 4.2.1 VT1 correlation-------------------------------------------------------
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

