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
setwd("~/Documents/JM/NanoString/NanoString_Covid/nanostring_covid/data") #folder data
rm(list = ls())
load("1.2_mat_pat_clean.rds") #ouverture de la svg
mat_pat_clean_sans_R_T<-mat_pat_clean[20:160,]
load("HVG_scran.rds") #ouverture de la svg

my_palette = colorRampPalette(c("royalblue4", "lightskyblue3", "white", "lightsalmon3","darkred"))(n = 256)

mat_REA <- mat_pat_clean$real_time_point %in% "REA" #Selection des prelevement REA creation d'un vecteur logique TRUE = REA
mat_REA <- mat_pat_clean[mat_REA==T,]  # on conserve les lignes qui on REA c-a-d TRUE
mat_T <- mat_pat_clean$real_time_point %in% "T" #Selection des prelevement T creation d'un vecteur logique TRUE = T
mat_T <- mat_pat_clean[mat_T==T,]  # on conserve les lignes qui on T c-a-d TRUE
mat_REA_T<-rbind(mat_REA,mat_T)

# 3-VT1----

## 3.1-heatmap all gene / all groups----

mat_VT1<-mat_pat_clean$real_time_point %in% "VT1" #
mat_VT1<-mat_pat_clean[mat_VT1==T,]
mat_VT1<-arrange(mat_VT1, REPONSE)
mat_VT1<-rbind(mat_VT1,mat_REA_T)
coldata_VT1<-as.data.frame(mat_VT1$REPONSE)
row.names(coldata_VT1)<-row.names(mat_VT1)

coldata_num_VT1<- mat_VT1$numero_patient
mat_VT1[,737:743]<-NULL
all(row.names(coldata_VT1) == row.names(mat_VT1))


mat_log_VT1<-log1p(t(mat_VT1)) # passage en log 1 p pour les heatmap
annC_VT1 <- data.frame(condition= coldata_VT1) # pour les annotation en col pour la heatmap
all(rownames(annC_VT1) == colnames(mat_log_VT1))  # pour les annotation en col pour la heatmap
annC_VT1<-rename(annC_VT1,Groupe=mat_VT1.REPONSE)

pheatmap(mat_log_VT1, 
         scale="row", 
         fontsize_row=1,
         fontsize_col = 5,
         annotation_colors = list(Groupe = c(NR = "#7570BE", 
                                                      R="#F15854", 
                                                      RP="#882255",
                                                      REA = "#0000FF", 
                                                      "T"= "#117733")),
         annotation_col = annC_VT1, color = my_palette, 
         cutree_rows = 1,
         main = "Heatmap au prélèvement VT1 sur tout les gènes avec clusterisation des patients", # autre titre : Heatmap au prélèvement VT1 sur tout les gènes
         labels_col = coldata_num_VT1,
         cluster_cols = T, # ou F
         cluster_rows = T,)

## 3.2-heatmap HVG / all groups----
top.hvgs<-as.data.frame(top.hvgs)    ##creation d'un DF
test <- rownames(mat_log_VT1) %in% top.hvgs$V1  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_top_hvgs_VT1<- mat_log_VT1[test==T,]   #ne garde en colone que les True

pheatmap(mat_top_hvgs_VT1,
         scale="row", 
         fontsize_row=1, 
         fontsize_col = 5,
         annotation_colors = list(Groupe = c(NR = "#7570BE",
                                                      R="#F15854", 
                                                      RP="#882255",
                                                      REA = "#0000FF",
                                                      "T"= "#117733")),
         annotation_col = annC_VT1, 
         color = my_palette,
         cutree_rows = 1, 
         main = "Heatmap au prélèvement VT1 sur les HVG avec clusterisation des patients", # autre titre si clust col F = Heatmap au prélèvement VT1 sur les HVG
         labels_col = coldata_num_VT1,
         cluster_cols = T, # ou F
         cluster_rows = T,)

# 4-VT2----

## 4.1-heatmap all gene / all groups----

mat_VT2<-mat_pat_clean$real_time_point %in% "VT2"
mat_VT2<-mat_pat_clean[mat_VT2==T,]
mat_VT2<-arrange(mat_VT2, REPONSE)
mat_VT2<-rbind(mat_VT2,mat_REA_T)
coldata_VT2<-as.data.frame(mat_VT2$REPONSE)
row.names(coldata_VT2)<-row.names(mat_VT2)

coldata_num_VT2<- mat_VT2$numero_patient
mat_VT2[,737:743]<-NULL
all(row.names(coldata_VT2) == row.names(mat_VT2))


mat_log_VT2<-log1p(t(mat_VT2)) # passage en log 1 p pour les heatmap
annC_VT2 <- data.frame(condition= coldata_VT2) # pour les annotation en col pour la heatmap
all(rownames(annC_VT2) == colnames(mat_log_VT2))  # pour les annotation en col pour la heatmap
annC_VT2<-rename(annC_VT2,Groupe=mat_VT2.REPONSE)

pheatmap(mat_log_VT2,
         scale="row",
         fontsize_row=1,
         fontsize_col = 5,
         annotation_colors = list(Groupe = c(NR = "#7570BE",
                                                      R="#F15854",
                                                      RP="#882255",
                                                      REA = "#0000FF",
                                                      "T"= "#117733")),
         annotation_col = annC_VT2,
         color = my_palette, 
         cutree_rows = 1,
         main = "Heatmap au prélèvement VT2 sur tout les gènes avec clusterisation des patients", #Heatmap au prélèvement VT2 sur tout les gènes 
         labels_col = coldata_num_VT2,
         cluster_cols = T, #F
         cluster_rows = T,)


## 4.2-heatmap HVG / all groups----
top.hvgs<-as.data.frame(top.hvgs)    ##creation d'un DF
test <- rownames(mat_log_VT2) %in% top.hvgs$V1  #creation d'un vecteur logique T F, les T correspondent au nom identique dans les deux
mat_top_hvgs_VT2<- mat_log_VT2[test==T,]   #ne garde en colone que les True

pheatmap(mat_top_hvgs_VT2, 
         scale="row",   
         fontsize_row=1, 
         fontsize_col = 5,
         annotation_colors = list(Groupe = c(NR = "#7570BE", 
                                                      R="#F15854", 
                                                      RP="#882255",
                                                      REA = "#0000FF",
                                                      "T"= "#117733")),
         annotation_col = annC_VT2, 
         color = my_palette,
         cutree_rows = 1, 
         main = "Heatmap au prélèvement VT2 sur les HVG avec clusterisation des patients",   #Heatmap au prélèvement VT2 sur les HVG
         labels_col = coldata_num_VT2,
         cluster_cols = T,   #F
         cluster_rows = T,)
