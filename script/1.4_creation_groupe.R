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


# 3-Création des groupes----

## 3.1-REA----
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

## 3.2-Temoins----
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

## 3.3-VT1----
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

NR= data.frame(  "13", "49","50", "62")
NR <- as_tibble(t(NR))
NR$classe = rep("NR",4)
colnames(NR) = c("numero_patient", "classe")

classe = rbind(classe,NR)

RP = data.frame(  "10", "21", "26", "29", "48", "61")
RP <- as_tibble(t(RP))
RP$classe = rep("RP",6)
colnames(RP) = c("numero_patient", "classe")

classe = rbind(classe, RP)

R = data.frame("1","2","3"  ,  "5" ,  "8"  ,"9",  "11","12","14", "16","17",  "19", "20" , "22" ,"23" ,"25",  "27" , "30" ,"31", "47" ,  "51",   "52", "54", "58","59","60",  "63")
R<- as_tibble(t(R))
R$classe <- rep("R", 27)
colnames(R)<- c("numero_patient", "classe")

classe = rbind(classe, R)

## 3.4-VT2----
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

NR= data.frame("13",  "49","50")
NR <- as_tibble(t(NR))
NR$classe = rep("NR",3)
colnames(NR) = c("numero_patient", "classe")
classe = rbind(NR, classe)

RP = data.frame( "4",  "6", "10", "15", "21", "26", "29","48","61")
RP <- as_tibble(t(RP))
RP$classe = rep("RP",9)
colnames(RP) = c("numero_patient", "classe")

classe = rbind(classe, RP)

R = data.frame("1","2" ,"3"  ,  "5"  ,"7",  "8"  , "9",  "11","12",   "14",  "16","17", "18", "19", "20" , "22" ,"23" ,"25",  "27" , "30" ,"31", "47" ,"51",  "52", "54", "58","59", "60",  "63")
R<- as_tibble(t(R))
R$classe <- rep("R", 29)
colnames(R)<- c("numero_patient", "classe")

classe = rbind(classe, R)

## 3.5-VT3----
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

NR= data.frame( "13","62")
NR <- as_tibble(t(NR))
NR$classe = rep("NR",2)
colnames(NR) = c("numero_patient", "classe")
classe = rbind(classe, NR)

RP = data.frame( "4", "6", "10","15", "21", "26", "29")
RP <- as_tibble(t(RP))
RP$classe = rep("RP",7)
colnames(RP) = c("numero_patient", "classe")

classe = rbind(classe, RP)

R = data.frame("1","2","3"  ,  "5",   "7"  ,  "8"  , "9",  "11","12",   "14",  "16", "17","18",  "19", "20" , "22" ,"23" ,"25",  "27" , "30" ,"31")
R<- as_tibble(t(R))
R$classe <- rep("R", 21)
colnames(R)<- c("numero_patient", "classe")

classe = rbind(classe, R)

## 3.6-VT4----
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

NR= data.frame(  "13")
NR <- as_tibble(t(NR))
NR$classe = rep("NR",1)
colnames(NR) = c("numero_patient", "classe")
classe = rbind(classe, NR)

RP = data.frame( "4","4",  "6",  "6", "10", "15","15", "21", "26")
RP <- as_tibble(t(RP))
RP$classe = rep("RP",9)
colnames(RP) = c("numero_patient", "classe")

classe = rbind(classe, RP)

R = data.frame("1","2","3"  ,  "5"  ,"7","7",  "8"  , "9",  "11","12",   "14",  "16","17", "18", "18",  "19", "20" , "22" ,"23","25",  "27" , "30" ,"31")
R<- as_tibble(t(R))
R$classe <- rep("R", 23)
colnames(R)<- c("numero_patient", "classe")

classe = rbind(classe, R)

# 4-matrice finale----

classe$numero_patient<- as.numeric(as.character(classe$numero_patient)) #pour ordo le vecteur lvl
classe <- arrange(classe, numero_patient) #ordo par num pat
mat_pat_clean<-arrange(mat_pat_clean, numero_patient)
all(mat_pat_clean$numero_patient == classe$numero_patient) #verif que les num patient sont dans le meme ordre
mat_pat_clean$REPONSE <- classe$classe #ajout de la colonne REPONSE (NR_R_RP_REA_T)
mat_pat_clean<-arrange(mat_pat_clean, real_time_point, REPONSE) #ordo real_time_point en premiere puis REPONSE

# 5-Svg fichier ----
save(mat_pat_clean, file = "data/1.2_mat_pat_clean.rds") #svg du data
load("data/1.2_mat_pat_clean.rds") #ouverture de la svg











coldata<-mat_pat_clean[,737:743] # creation du coldata pour les heatmap
coldata_num<- coldata$numero_patient  # ordo
coldata[,1:6]<-NULL # sup des data non utile
mat<-mat_pat_clean #pour la mat de count
mat[,737:743]<-NULL # sup des data non utile
mat<-t(mat)  # translation de la matrice
all(rownames(coldata) == colnames(mat)) #verif col et row identique
all(rownames(coldata) %in% colnames(mat)) #attention de verif pas l'ordre mais seulement présence
mat_log<-log1p(mat) # passage en log 1 p pour les heatmap
annC <- data.frame(condition= coldata) # pour les annotation en col pour la heatmap
rownames(annC) <- colnames(mat_log)  # pour les annotation en col pour la heatmap
all(rownames(annC) == colnames(mat)) #vérif du même ordre / position
all(rownames(annC) == rownames(mat_pat_clean))#vérif du même ordre / position
annC$real_time_point <- mat_pat_clean$real_time_point  # ajout real_time_point

my_palette = colorRampPalette(c("royalblue4", "lightskyblue3", "white", "lightsalmon3","darkred"))(n = 256)

pheatmap(mat_log, scale="row", fontsize_row=1, fontsize_col = 5,annotation_col = annC, color = my_palette, cutree_rows = 1, main = "Heatmap all prelevement all time point all gene", labels_col = coldata_num,cluster_cols = T,cluster_rows = T,)

pheatmap(mat_log, scale="row", fontsize_row=1, fontsize_col = 5,annotation_col = annC, color = my_palette, cutree_rows = 1, main = "Heatmap all prelevement all time point all gene", labels_col = coldata_num,cluster_cols = F,cluster_rows = T,)


classe_SANS_REA_T<-classe[classe$classe !="REA",]
classe_SANS_REA_T<-classe[classe$classe !="T",]
