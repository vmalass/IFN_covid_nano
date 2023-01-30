# 1-library --------------------------------------------------------------------
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

# 2-ouverture des fichier ------------------------------------------------------
rm(list = ls())
load("data/1.1_mat_pat_clean.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]


# 3-Création des groupes -------------------------------------------------------
## 3.1-REA ---------------------------------------------------------------------
mat_REA <- mat_pat_clean$real_time_point %in% "REA" #Selection des prelevement REA creation d'un vecteur logique TRUE = REA
count_matrix <- mat_pat_clean[mat_REA==T,]  # on conserve les lignes qui on REA c-a-d TRUE
count_matrix<-arrange(count_matrix,numero_patient)
count_matrix <- t(count_matrix[,1:736])

coldata <- mat_pat_clean[,737:742]
mat_REA <- coldata$real_time_point %in% "REA" #Selection des prelevement REA creation d'un vecteur logique TRUE = REA
coldata <- coldata[mat_REA==T,]  # on conserve les lignes qui on REA c-a-d TRUE

coldata<-arrange(coldata,numero_patient)

all(row.names(coldata) == colnames(count_matrix))  # vérif que les patients sont dans le meme ordre pour les deux tables

coldata_num_REA<- coldata$numero_patient
REA <- data.frame(as.character(coldata_num_REA))
REA$classe = rep("REA",12)
colnames(REA) = c("numero_patient", "classe")

## 3.2-Temoins -----------------------------------------------------------------
mat_T <- mat_pat_clean$real_time_point %in% "T" #Selection des prelevement T creation d'un vecteur logique TRUE = T
count_matrix <- mat_pat_clean[mat_T==T,]  # on conserve les lignes qui on T c-a-d TRUE
count_matrix<-arrange(count_matrix,numero_patient)
count_matrix <- t(count_matrix[,1:736])

coldata <- mat_pat_clean[,737:742]
mat_T <- coldata$real_time_point %in% "T" #Selection des prelevement T creation d'un vecteur logique TRUE = T
coldata <- coldata[mat_T==T,]  # on conserve les lignes qui on T c-a-d TRUE

coldata<-arrange(coldata,numero_patient)

all(row.names(coldata) == colnames(count_matrix))  # vérif que les patients sont dans le meme ordre pour les deux tables

coldata_num_T<- coldata$numero_patient

TEMOIN= data.frame( as.character(coldata_num_T))
TEMOIN$classe = rep("T",7)
colnames(TEMOIN) = c("numero_patient", "classe")

classe = rbind(REA,TEMOIN)

## 3.3-VT1 ---------------------------------------------------------------------
mat_V1 <- mat_pat_clean$real_time_point %in% "VT1" #Selection des prelevement V1 creation d'un vecteur logique TRUE = V1
count_matrix <- mat_pat_clean[mat_V1==T,]  # on conserve les lignes qui on V1 c-a-d TRUE
count_matrix<-arrange(count_matrix,numero_patient)
count_matrix <- t(count_matrix[,1:736])

coldata <- mat_pat_clean[,737:742]
mat_V1 <- coldata$real_time_point %in% "VT1" #Selection des prelevement V1 creation d'un vecteur logique TRUE = V1
coldata <- coldata[mat_V1==T,]  # on conserve les lignes qui on V1 c-a-d TRUE

coldata<-arrange(coldata,numero_patient)

all(row.names(coldata) == colnames(count_matrix))  # vérif que les patients sont dans le meme ordre pour les deux tables

coldata_num_V1<- coldata$numero_patient

NR= data.frame( "49","50")
NR <- as_tibble(t(NR))
NR$classe = rep("NR",2)
colnames(NR) = c("numero_patient", "classe")

classe = rbind(classe,NR)

NRm= data.frame( "2", "13", "62")
NRm <- as_tibble(t(NRm))
NRm$classe = rep("NR-",3)
colnames(NRm) = c("numero_patient", "classe")

classe = rbind(classe,NRm)

R = data.frame("1", "3", "5" , "6", "8", "9", "11", "12", "14", "16", "17", "20",
               "22","25", "27", "30", "31", "47", "51", "54", "58", "63")
R<- as_tibble(t(R))
R$classe <- rep("R", 22)
colnames(R)<- c("numero_patient", "classe")

classe = rbind(classe, R)

RPm = data.frame(  "19", "23", "52", "59")
RPm <- as_tibble(t(RPm))
RPm$classe = rep("RP-",4)
colnames(RPm) = c("numero_patient", "classe")

classe = rbind(classe, RPm)

RP = data.frame(  "10", "21", "26", "29", "48", "61")
RP <- as_tibble(t(RP))
RP$classe = rep("RP",6)
colnames(RP) = c("numero_patient", "classe")

classe = rbind(classe, RP)

A = data.frame("60")
A<- as_tibble(t(A))
A$classe <- rep("A", 1)
colnames(A)<- c("numero_patient", "classe")

classe = rbind(classe, A)


## 3.4-VT2 ---------------------------------------------------------------------
mat_V2 <- mat_pat_clean$real_time_point %in% "VT2" #Selection des prelevement V1 creation d'un vecteur logique TRUE = V1
count_matrix <- mat_pat_clean[mat_V2==T,]  # on conserve les lignes qui on V1 c-a-d TRUE
count_matrix<-arrange(count_matrix,numero_patient)
count_matrix <- t(count_matrix[,1:736])

coldata <- mat_pat_clean[,737:742]
mat_V2 <- coldata$real_time_point %in% "VT2" #Selection des prelevement V1 creation d'un vecteur logique TRUE = V1
coldata <- coldata[mat_V2==T,]  # on conserve les lignes qui on V1 c-a-d TRUE

coldata<-arrange(coldata,numero_patient)

all(row.names(coldata) == colnames(count_matrix))  # vérif que les patients sont dans le meme ordre pour les deux tables
coldata_num_V2<- coldata$numero_patient

NR= data.frame( "49","50")
NR <- as_tibble(t(NR))
NR$classe = rep("NR",2)
colnames(NR) = c("numero_patient", "classe")

classe = rbind(classe,NR)

NRm= data.frame( "2", "13")
NRm <- as_tibble(t(NRm))
NRm$classe = rep("NR-",2)
colnames(NRm) = c("numero_patient", "classe")

classe = rbind(classe,NRm)

R = data.frame("1", "3", "5" , "6", "7", "8", "9", "11", "12", "14", "16", "17",
               "18", "20", "22", "24", "24", "25", "27", "30", "31", "47", "51", 
               "54", "58", "63")
R<- as_tibble(t(R))
R$classe <- rep("R", 26)
colnames(R)<- c("numero_patient", "classe")

classe = rbind(classe, R)

RPm = data.frame("19", "23", "52", "59")
RPm <- as_tibble(t(RPm))
RPm$classe = rep("RP-",4)
colnames(RPm) = c("numero_patient", "classe")

classe = rbind(classe, RPm)

RP = data.frame( "4", "10", "15", "21", "26", "29", "48", "61")
RP <- as_tibble(t(RP))
RP$classe = rep("RP",8)
colnames(RP) = c("numero_patient", "classe")

classe = rbind(classe, RP)

A = data.frame("60")
A<- as_tibble(t(A))
A$classe <- rep("A", 1)
colnames(A)<- c("numero_patient", "classe")

classe = rbind(classe, A)

## 3.5-VT3 ---------------------------------------------------------------------
mat_V3 <- mat_pat_clean$real_time_point %in% "VT3" #Selection des prelevement V1 creation d'un vecteur logique TRUE = V1
count_matrix <- mat_pat_clean[mat_V3==T,]  # on conserve les lignes qui on V1 c-a-d TRUE
count_matrix<-arrange(count_matrix,numero_patient)
count_matrix <- t(count_matrix[,1:736])

coldata <- mat_pat_clean[,737:742]
mat_V3 <- coldata$real_time_point %in% "VT3" #Selection des prelevement V1 creation d'un vecteur logique TRUE = V1
coldata <- coldata[mat_V3==T,]  # on conserve les lignes qui on V1 c-a-d TRUE

coldata<-arrange(coldata,numero_patient)

all(row.names(coldata) == colnames(count_matrix))  # vérif que les patients sont dans le meme ordre pour les deux tables
coldata_num_V3<- coldata$numero_patient

NRm= data.frame( "2", "13", "62")
NRm <- as_tibble(t(NRm))
NRm$classe = rep("NR-",3)
colnames(NRm) = c("numero_patient", "classe")

classe = rbind(classe,NRm)

R = data.frame("1", "3", "5" , "6", "7", "8", "9", "11", "12", "14", "16", "17",
               "18", "20", "22", "24", "25", "27", "30", "31")
R<- as_tibble(t(R))
R$classe <- rep("R", 20)
colnames(R)<- c("numero_patient", "classe")

classe = rbind(classe, R)

RPm = data.frame(  "19", "23")
RPm <- as_tibble(t(RPm))
RPm$classe = rep("RP-",2)
colnames(RPm) = c("numero_patient", "classe")

classe = rbind(classe, RPm)

RP = data.frame( "4", "10", "15", "21", "26", "29")
RP <- as_tibble(t(RP))
RP$classe = rep("RP",6)
colnames(RP) = c("numero_patient", "classe")

classe = rbind(classe, RP)


## 3.6-VT4 ---------------------------------------------------------------------
mat_V4 <- mat_pat_clean$real_time_point %in% "VT4" #Selection des prelevement V1 creation d'un vecteur logique TRUE = V1
count_matrix <- mat_pat_clean[mat_V4==T,]  # on conserve les lignes qui on V1 c-a-d TRUE
count_matrix<-arrange(count_matrix,numero_patient)
count_matrix <- t(count_matrix[,1:736])

coldata <- mat_pat_clean[,737:742]
mat_V4 <- coldata$real_time_point %in% "VT4" #Selection des prelevement V1 creation d'un vecteur logique TRUE = V1
coldata <- coldata[mat_V4==T,]  # on conserve les lignes qui on V1 c-a-d TRUE

coldata<-arrange(coldata,numero_patient)

all(row.names(coldata) == colnames(count_matrix))  # vérif que les patients sont dans le meme ordre pour les deux tables
coldata_num_V4<- coldata$numero_patient

NRm= data.frame( "2", "13")
NRm <- as_tibble(t(NRm))
NRm$classe = rep("NR-",2)
colnames(NRm) = c("numero_patient", "classe")

classe = rbind(classe,NRm)

R = data.frame("1", "3", "5" , "6", "7", "7", "8", "9", "11", "12", "14", "16",
               "17", "18", "18", "20", "22", "25", "27", "30", "31")
R<- as_tibble(t(R))
R$classe <- rep("R", 21)
colnames(R)<- c("numero_patient", "classe")

classe = rbind(classe, R)

RPm = data.frame(  "19", "23")
RPm <- as_tibble(t(RPm))
RPm$classe = rep("RP-",2)
colnames(RPm) = c("numero_patient", "classe")

classe = rbind(classe, RPm)

RP = data.frame( "4", "4", "10", "15", "15", "21", "26")
RP <- as_tibble(t(RP))
RP$classe = rep("RP",7)
colnames(RP) = c("numero_patient", "classe")

classe = rbind(classe, RP)

# 4-matrice finale -------------------------------------------------------------
classe$numero_patient<- as.numeric(as.character(classe$numero_patient)) #pour ordo le vecteur lvl
classe <- arrange(classe, numero_patient) #ordo par num pat
mat_pat_clean<-arrange(mat_pat_clean, numero_patient)
all(mat_pat_clean$numero_patient == classe$numero_patient) #verif que les num patient sont dans le meme ordre
mat_pat_clean$REPONSE <- classe$classe #ajout de la colonne REPONSE (NR_R_RP_REA_T)
mat_pat_clean<-arrange(mat_pat_clean, real_time_point, REPONSE) #ordo real_time_point en premiere puis REPONSE

# 5-Svg fichier ----------------------------------------------------------------
save(mat_pat_clean, file = "data/1.3_mat_pat_clean_final.rds") #svg du data
load("data/1.3_mat_pat_clean_final.rds") #ouverture de la svg
