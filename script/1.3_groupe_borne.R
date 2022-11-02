# 1-library -----
library(openxlsx)
library(umap)
library(factoextra)
library(S4Vectors)
library(magrittr)
library(DESeq2)
library(ggplot2)
library(scater)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(clustree)
library(fpc)
library(scran)
library(pheatmap)
library(RColorBrewer)
library(gplots)
library(dplyr)
library(forcats)
library(scales)
library(pheatmap)
library(apeglm)

# 2-ouverture des fichier----
rm(list = ls())
load("data/1.1_mat_pat_clean.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]

# 3-Bornage VT1----
## 3.1-bornage VT1 / RSAD2 / real_time_point -----
vec_VT1 <- mat_pat_clean$real_time_point %in% "VT1" #vec logique avec True = VT1
borne_VT1<- mat_pat_clean[vec_VT1==T,] #nouvel obj avec que les ligne qui ont VT1 en real time point
borne_VT1<- arrange(borne_VT1,numero_patient) #ordo par n° pat

round(median(borne_VT1$RSAD2),5) # calcul de la median poàur le gene RSAD2
# median = 26527
round(mad(borne_VT1$RSAD2, low = F , high = F),5) # calcul du MAD (median absolut deviation) (= sd pour la moyene)
# MAD = 16978.74
borne_VT1$RSAD2 < median(borne_VT1$RSAD2)-1.2*mad(borne_VT1$RSAD2, low = F , high = F) 

groupe1_J1 <- subset(borne_VT1, borne_VT1$RSAD2 < median(borne_VT1$RSAD2)-1.2*mad(borne_VT1$RSAD2, low = F , high = F)) 

num_pat_J1 <- groupe1_J1$numero_patient # 2  13 49 50 62

### 3.1.1- RSAD2 / VT1----
new <- mat_pat_clean$numero_patient %in% num_pat_J1
groupe_RP_autre <- mat_pat_clean[new==F,]

numgp_RP <- groupe_RP_autre$numero_patient
GP3 <- mat_pat_clean_sans_R_T$numero_patient %in% numgp_RP
GP3 <- mat_pat_clean_sans_R_T[GP3==F,]

ggplot(data = GP3, 
       aes_string(x="jours_prelevement", 
                              y= "RSAD2",
                              color = "numero_patient ",
                              group = "numero_patient")) + 
  geom_line() + 
  geom_point() + 
  ggtitle("Patients NR avec borne J1 à J5 sur le gène RSAD2") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

YY <- as.character(num_pat_J1)

GP5 <- GP3 %>% mutate(numero_patient  = as.character(numero_patient))
GP5$numero_patient
couleur_patient <- mat_pat_clean_sans_R_T$numero_patient %in% YY
couleur_patient
GP5<- cbind(mat_pat_clean_sans_R_T, couleur_patient)

ggplot(data = GP5, 
       aes_string(x="jours_prelevement", 
                  y= "RSAD2", 
                  color = "couleur_patient ", 
                  group = "numero_patient")) + 
  geom_line() + 
  geom_point() + 
  ggtitle("Patients NR avec borne J1 à J5 sur le gène RSAD2") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))


## 3.2-bornage VT1 / CXCL10 / real_time_point -----
round(median(borne_VT1$CXCL10),5)
round(mad(borne_VT1$CXCL10, low = F , high = F),5)

borne_VT1$CXCL10 < median(borne_VT1$CXCL10)-0.8*mad(borne_VT1$CXCL10, low = F , high = F)

groupe1_J1 <- subset(borne_VT1, borne_VT1$CXCL10 < median(borne_VT1$CXCL10)-0.8*mad(borne_VT1$CXCL10, low = F , high = F))  
num_pat_J1 <- groupe1_J1$numero_patient  #13 49 50 51 62

### 3.2.1-CXCL10 / VT1-----
new <- mat_pat_clean$numero_patient %in% num_pat_J1
groupe_RP_autre <- mat_pat_clean[new==F,]

numgp_RP <- groupe_RP_autre$numero_patient
GP3 <- mat_pat_clean_sans_R_T$numero_patient %in% numgp_RP
GP3 <- mat_pat_clean_sans_R_T[GP3==F,]

ggplot(data = GP3, aes_string(x="jours_prelevement", y= "CXCL10", color = "numero_patient ", group = "numero_patient"))+ geom_line() + geom_point() + ggtitle("Patients NR avec borne J1 à J5 sur le gène CXCL10")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

YY <- as.character(num_pat_J1)

GP5 <- GP3 %>% mutate(numero_patient  = as.character(numero_patient))
GP5$numero_patient
couleur_patient <- mat_pat_clean_sans_R_T$numero_patient %in% YY
couleur_patient
GP5<- cbind(mat_pat_clean_sans_R_T, couleur_patient)

ggplot(data = GP5, aes_string(x="jours_prelevement", y= "CXCL10", color = "couleur_patient ", group = "numero_patient"))+ geom_line() + geom_point() + ggtitle("Patients NR avec borne J1 à J5 sur le gène CXCL10")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

## 3.3-bornage VT1 / IFI27 / real_time_point -----
round(median(borne_VT1$IFI27),5)
round(mad(borne_VT1$IFI27, low = F , high = F),5)

borne_VT1$IFI27 < median(borne_VT1$IFI27)-0.8*mad(borne_VT1$IFI27, low = F , high = F)

groupe1_J1 <- subset(borne_VT1, borne_VT1$IFI27 < median(borne_VT1$IFI27)-0.9*mad(borne_VT1$IFI27, low = F , high = F)) 
num_pat_J1 <- groupe1_J1$numero_patient #10 49


### 3.3.1-IFI27 / VT1-----
new <- mat_pat_clean$numero_patient %in% num_pat_J1
groupe_RP_autre <- mat_pat_clean[new==F,]

numgp_RP <- groupe_RP_autre$numero_patient
GP3 <- mat_pat_clean_sans_R_T$numero_patient %in% numgp_RP
GP3 <- mat_pat_clean_sans_R_T[GP3==F,]

ggplot(data = GP3, aes_string(x="jours_prelevement", y= "IFI27", color = "numero_patient ", group = "numero_patient"))+ geom_line() + geom_point() + ggtitle("Patients NR avec borne J1 à J6 sur le gène IFI27")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

YY <- as.character(num_pat_J1)

GP5 <- GP3 %>% mutate(numero_patient  = as.character(numero_patient))
GP5$numero_patient
couleur_patient <- mat_pat_clean_sans_R_T$numero_patient %in% YY
couleur_patient
GP5<- cbind(mat_pat_clean_sans_R_T, couleur_patient)

ggplot(data = GP5, aes_string(x="jours_prelevement", y= "IFI27", color = "couleur_patient ", group = "numero_patient"))+ geom_line() + geom_point() + ggtitle("Patients NR avec borne J1 à J6 sur le gène IFI27")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

# 4-Bornage VT2 -----

## 4.1-bornage VT2 / RSAD2 / real_time_point -----
vec_VT2 <- mat_pat_clean$real_time_point %in% "VT2" #vec logique avec True = VT2
borne_VT2<- mat_pat_clean[vec_VT2==T,] #nouvel obj avec que les ligne qui ont VT2 en real time point
borne_VT2<- arrange(borne_VT2,numero_patient) #ordo par n° pat

round(median(borne_VT2$RSAD2),5)
round(mad(borne_VT2$RSAD2, low = F , high = F),5)

borne_VT2$RSAD2 > median(borne_VT2$RSAD2)+3*mad(borne_VT2$RSAD2, low = F , high = F)

groupe1_J7 <- subset(borne_VT2, borne_VT2$RSAD2 > median(borne_VT2$RSAD2)+3*mad(borne_VT2$RSAD2, low = F , high = F))
num_pat_J7 <- groupe1_J7$numero_patient # 4  6  10 15 21 23 26 29 48 60 61

### 4.1.1-RSAD2 / VT2-----
new <- mat_pat_clean$numero_patient %in% num_pat_J7
groupe_RP_autre <- mat_pat_clean[new==F,]

numgp_RP <- groupe_RP_autre$numero_patient
GP3 <- mat_pat_clean_sans_R_T$numero_patient %in% numgp_RP
GP3 <- mat_pat_clean_sans_R_T[GP3==F,]

ggplot(data = GP3, aes_string(x="jours_prelevement", y= "RSAD2", color = "numero_patient ", group = "numero_patient"))+ geom_line() + geom_point() + ggtitle("Patients RP avec borne J7 à J13 sur le gène RSAD2")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

YY <- as.character(num_pat_J7)

GP5 <- GP3 %>% mutate(numero_patient  = as.character(numero_patient))
GP5$numero_patient
couleur_patient <- mat_pat_clean_sans_R_T$numero_patient %in% YY
couleur_patient
GP5<- cbind(mat_pat_clean_sans_R_T, couleur_patient)

ggplot(data = GP5, aes_string(x="jours_prelevement", y= "RSAD2", color = "couleur_patient ", group = "numero_patient"))+ geom_line() + geom_point() + ggtitle("Patients RP avec borne J7 à J13 sur le gène RSAD2")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

## 4.2- bornage VT2 / CXCL10 / real_time_point-----
round(median(borne_VT2$CXCL10),5)
round(mad(borne_VT2$CXCL10, low = F , high = F),5)

borne_VT2$CXCL10 > median(borne_VT2$CXCL10)+3*mad(borne_VT2$CXCL10, low = F , high = F)

groupe1_J7 <- subset(borne_VT2, borne_VT2$CXCL10 > median(borne_VT2$CXCL10)+3*mad(borne_VT2$CXCL10, low = F , high = F))  
num_pat_J7 <- groupe1_J7$numero_patient #4  6  10 15 21 26 29 48 61

### 4.2.1- CXCL10 / VT2-----
new <- mat_pat_clean$numero_patient %in% num_pat_J7
groupe_RP_autre <- mat_pat_clean[new==F,]

numgp_RP <- groupe_RP_autre$numero_patient
GP3 <- mat_pat_clean_sans_R_T$numero_patient %in% numgp_RP
GP3 <- mat_pat_clean_sans_R_T[GP3==F,]

ggplot(data = GP3, aes_string(x="jours_prelevement", y= "CXCL10", color = "numero_patient ", group = "numero_patient"))+ geom_line() + geom_point() + ggtitle("Patients RP avec borne J7 à J14 sur le gène CXCL10")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

YY <- as.character(num_pat_J7)

GP5 <- GP3 %>% mutate(numero_patient  = as.character(numero_patient))
GP5$numero_patient
couleur_patient <- mat_pat_clean_sans_R_T$numero_patient %in% YY
couleur_patient
GP5<- cbind(mat_pat_clean_sans_R_T, couleur_patient)

ggplot(data = GP5, aes_string(x="jours_prelevement", y= "CXCL10", color = "couleur_patient ", group = "numero_patient"))+ geom_line() + geom_point() + ggtitle("Patients RP avec borne J7 à J14 sur le gène CXCL10")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

## 4.3- bornage VT2 / IFI27 / real_time_point-----
round(median(borne_VT2$IFI27),5)
round(mad(borne_VT2$IFI27, low = F , high = F),5)

borne_VT2$IFI27 > median(borne_VT2$IFI27)+1.8*mad(borne_VT2$IFI27, low = F , high = F)

groupe1_J7 <- subset(borne_VT2, borne_VT2$IFI27 > median(borne_VT2$IFI27)+1.8*mad(borne_VT2$IFI27, low = F , high = F))  
num_pat_J7 <- groupe1_J7$numero_patient # 1  10 15 48 59

### 4.3.1-IFI27 / VT2-----
new <- mat_pat_clean$numero_patient %in% num_pat_J7
groupe_RP_autre <- mat_pat_clean[new==F,]

numgp_RP <- groupe_RP_autre$numero_patient
GP3 <- mat_pat_clean_sans_R_T$numero_patient %in% numgp_RP
GP3 <- mat_pat_clean_sans_R_T[GP3==F,]

ggplot(data = GP3, aes_string(x="jours_prelevement", y= "IFI27", color = "numero_patient ", group = "numero_patient"))+ geom_line() + geom_point() + ggtitle("Patients RP avec borne J7 à J13 sur le gène IFI27 1.8MAD")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))

YY <- as.character(num_pat_J7)

GP5 <- GP3 %>% mutate(numero_patient  = as.character(numero_patient))
GP5$numero_patient
couleur_patient <- mat_pat_clean_sans_R_T$numero_patient %in% YY
couleur_patient
GP5<- cbind(mat_pat_clean_sans_R_T, couleur_patient)

ggplot(data = GP5, aes_string(x="jours_prelevement", y= "IFI27", color = "couleur_patient ", group = "numero_patient"))+ geom_line() + geom_point() + ggtitle("Patients RP avec borne J7 à J13 sur le gène IFI27 1.8MAD")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))


# Pour VT1 et VT2 on selectionne tout les patients qui ressorte au moins dans deux gènes pour les attribuer aux groupes NR & RP
# NR 13 / 49 / 50 / 62
# RP 4 / 6 / 10 / 15 / 21 / 26 / 29 / 48 / 61

